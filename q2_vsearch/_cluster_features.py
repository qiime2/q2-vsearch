# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import subprocess
import sqlite3

import biom
import skbio
import pandas as pd
from qiime2 import Metadata
from q2_types.feature_data import DNAFASTAFormat


class VSearchError(Exception):
    pass


def run_command(cmd, verbose=True):
    print("Running external command line application. This may print "
          "messages to stdout and/or stderr.")
    print("The command being run is below. This command cannot "
          "be manually re-run as it will depend on temporary files that "
          "no longer exist.")
    print("\nCommand:", end=' ')
    print(" ".join(cmd), end='\n\n')
    subprocess.run(cmd, check=True)


def _uc_to_sqlite(uc):
    '''Parse uc-style file into a SQLite in-memory database.

    This populates an in-memory database with the following schema (displayed
    below with dummy data):

        feature_id | cluster_id | count
        -----------|------------|-------
        feature1   | r1         | 204
        feature2   | r2         | 4
        feature3   | r1         | 15
        feature4   | r2         | 24
        feature5   | r2         | 16
    '''
    conn = sqlite3.connect(':memory:')
    c = conn.cursor()
    # The PK constraint ensures that there are no duplicate Feature IDs
    c.execute('CREATE TABLE feature_cluster_map (feature_id TEXT PRIMARY KEY,'
              'cluster_id TEXT NOT NULL, count INTEGER);')
    c.execute('CREATE INDEX idx1 ON '
              'feature_cluster_map(feature_id, cluster_id);')
    conn.commit()
    insert_stmt = 'INSERT INTO feature_cluster_map VALUES (?, ?, ?);'

    for line in uc:
        line = line.strip()
        if len(line) == 0 or line.startswith(b'#'):
            continue
        else:
            fields = line.split(b'\t')
            if fields[0] == b'S':
                sequence_id = fields[8].decode('utf-8').split(';')[0]
                c.execute(insert_stmt, (sequence_id, sequence_id, None))
            elif fields[0] == b'H':
                centroid_id = fields[9].decode('utf-8').split(';')[0]
                sequence_id = fields[8].decode('utf-8').split(';size=')
                sequence_id, count = sequence_id[0], sequence_id[1]
                c.execute(insert_stmt, (sequence_id, centroid_id, count))
            else:
                pass
    conn.commit()
    return conn


def _collapse_f_from_sqlite(conn):
    c = conn.cursor()
    # This query produces the following results (displayed below with dummy
    # data):
    # feature_id | cluster_id
    # -----------|------------
    # feature1   | r1
    # feature2   | r2
    # feature3   | r1
    # feature4   | r2
    # feature4   | r2
    c.execute('SELECT feature_id, cluster_id FROM feature_cluster_map;')
    id_to_centroid = dict(c.fetchall())

    if len(id_to_centroid) == 0:
        raise ValueError("No sequence matches were identified by vsearch.")

    def collapse_f(id_, x):
        return id_to_centroid[id_]

    return collapse_f


def _fasta_from_sqlite(conn, input_fasta_fp, output_fasta_fp):
    input_seqs = skbio.read(input_fasta_fp, format='fasta',
                            constructor=skbio.DNA)
    c = conn.cursor()
    # Create a second in-memory table with the following schema (displayed
    # below with dummy data):
    # feature_id | sequence_string
    # -----------|------------------
    # feature1   | ACGTACGTACGTACGT
    # feature2   | GGGGAAAACCCCTTTT
    # feature3   | TCAGAAAATTTTTCAG
    # feature4   | AAAAAAAAAAAAAAAA
    # feature5   | GGGGGGGGGGGGGGGG
    c.execute('CREATE TABLE rep_seqs (feature_id TEXT PRIMARY KEY, '
              'sequence_string TEXT NOT NULL);')
    c.executemany(
        'INSERT INTO rep_seqs VALUES (?, ?);',
        [(seq.metadata['id'], str(seq)) for seq in input_seqs]
    )
    conn.commit()
    # Preemptively sort the table to deal with tie-breaking, later.
    # This is a table, not a view, because we want/need sqlite's rowid.
    c.execute('CREATE TABLE sorted_feature_cluster_map AS '
              'SELECT * FROM feature_cluster_map ORDER BY cluster_id ASC,'
              'feature_id ASC;')
    c.execute('CREATE INDEX idx2 ON '
              'sorted_feature_cluster_map(cluster_id, count);')
    conn.commit()
    # The results from this query should look like the following (displayed
    # below with dummy data):
    # cluster_id | sequence_string
    # -----------|------------------
    # r1         | ACGTACGTACGTACGT
    # r2         | AAAAAAAAAAAAAAAA
    c.execute('''SELECT fcm.cluster_id, rs.sequence_string, MAX(fcm.count)
                   FROM sorted_feature_cluster_map fcm
             INNER JOIN rep_seqs rs ON rs.feature_id = fcm.feature_id
               GROUP BY fcm.cluster_id
               ORDER BY fcm.cluster_id ASC;
    ''')
    with open(output_fasta_fp, 'w') as output_seqs:
        while True:
            partial_results = c.fetchmany(size=100)
            if partial_results:
                output_seqs.writelines(
                    ['>%s\n%s\n' % (i, s) for (i, s, _) in partial_results]
                )
            else:
                break


def _fasta_with_sizes(input_fasta_fp, output_fasta_fp, table):
    table_ids = table.ids(axis='observation')
    sizes = {id_: size for id_, size in zip(table_ids,
                                            table.sum(axis='observation'))}
    output_fasta_f = open(output_fasta_fp, 'w')
    sequence_ids = set()
    for e in skbio.io.read(input_fasta_fp, constructor=skbio.DNA,
                           format='fasta'):
        feature_id = e.metadata['id']
        feature_seq = str(e)
        sequence_ids.add(feature_id)
        try:
            feature_size = sizes[feature_id]
        except KeyError:
            raise ValueError('Feature %s is present in sequences, but not '
                             'in table. The set of features in sequences must '
                             'be identical to the set of features in table.'
                             % feature_id)
        output_fasta_f.write('>%s;size=%d\n%s\n' %
                             (feature_id, feature_size, feature_seq))
    output_fasta_f.close()

    _error_on_nonoverlapping_ids(set(table_ids), sequence_ids,
                                 check_extra_table_ids=True,
                                 check_extra_sequence_ids=False)


def cluster_features_de_novo(sequences: DNAFASTAFormat, table: biom.Table,
                             perc_identity: float, threads: int = 1
                             ) -> (biom.Table, DNAFASTAFormat):
    clustered_sequences = DNAFASTAFormat()
    with tempfile.NamedTemporaryFile() as fasta_with_sizes:
        with tempfile.NamedTemporaryFile() as out_uc:
            _fasta_with_sizes(str(sequences), fasta_with_sizes.name, table)
            cmd = ['vsearch',
                   '--cluster_size', fasta_with_sizes.name,
                   '--id', str(perc_identity),
                   '--centroids', str(clustered_sequences),
                   '--uc', out_uc.name,
                   '--qmask', 'none',  # ensures no lowercase DNA chars
                   '--xsize',
                   '--threads', str(threads),
                   '--minseqlength', '1',
                   '--fasta_width', '0']
            run_command(cmd)
            out_uc.seek(0)

            conn = _uc_to_sqlite(out_uc)
            collapse_f = _collapse_f_from_sqlite(conn)

    table = table.collapse(collapse_f, norm=False, min_group_size=1,
                           axis='observation',
                           include_collapsed_metadata=False)

    return table, clustered_sequences


def _error_on_nonoverlapping_ids(table_ids, sequence_ids,
                                 check_extra_table_ids=True,
                                 check_extra_sequence_ids=True):
    if check_extra_table_ids:
        extra_table_ids = table_ids - sequence_ids
        if len(extra_table_ids):
            raise ValueError('Some feature ids are present in table, but not '
                             'in sequences. The set of features in sequences '
                             'must be identical to the set of features in '
                             'table. Feature ids present in table but not '
                             'sequences are: %s' % ', '.join(extra_table_ids))

    if check_extra_sequence_ids:
        extra_sequence_ids = sequence_ids - table_ids
        if len(extra_sequence_ids):
            raise ValueError('Some feature ids are present in sequences, but '
                             'not in table. The set of features in sequences '
                             'must be identical to the set of features in '
                             'table. Feature ids present in sequences but not '
                             'table are: %s' % ', '.join(extra_sequence_ids))


def cluster_features_closed_reference(sequences: DNAFASTAFormat,
                                      table: biom.Table,
                                      reference_sequences: DNAFASTAFormat,
                                      perc_identity: float,
                                      strand: str = 'plus',
                                      threads: int = 1
                                      ) -> (biom.Table, DNAFASTAFormat,
                                            DNAFASTAFormat):

    table_ids = set(table.ids(axis='observation'))
    sequence_ids = {e.metadata['id'] for e in skbio.io.read(
                    str(sequences), constructor=skbio.DNA, format='fasta')}
    _error_on_nonoverlapping_ids(table_ids, sequence_ids)
    matched_seqs, unmatched_seqs = DNAFASTAFormat(), DNAFASTAFormat()

    with tempfile.NamedTemporaryFile() as fasta_with_sizes, \
            tempfile.NamedTemporaryFile() as out_uc, \
            tempfile.NamedTemporaryFile() as tmp_unmatched_seqs:
        _fasta_with_sizes(str(sequences), fasta_with_sizes.name, table)
        cmd = ['vsearch',
               '--usearch_global', fasta_with_sizes.name,
               '--id', str(perc_identity),
               '--db', str(reference_sequences),
               '--uc', out_uc.name,
               '--strand', str(strand),
               '--qmask', 'none',  # ensures no lowercase DNA chars
               '--notmatched', tmp_unmatched_seqs.name,
               '--threads', str(threads),
               '--minseqlength', '1',
               '--fasta_width', '0']
        run_command(cmd)
        out_uc.seek(0)

        # It is possible for there to be no unmatched sequences --- if that
        # is the case, skip thie following clean-up.
        if os.path.getsize(tmp_unmatched_seqs.name) > 0:
            # We don't really need to sort the matched sequences, this
            # is just to let us use --xsize, which strips the counts from
            # the Feature ID. It would be more ideal if --usearch_global,
            # above let us pass in --xsize, but unfortunately it isn't
            # supported.
            cmd = ['vsearch',
                   '--sortbysize', tmp_unmatched_seqs.name,
                   '--xsize',
                   '--output', str(unmatched_seqs),
                   '--minseqlength', '1',
                   '--fasta_width', '0']
            run_command(cmd)

        try:
            conn = _uc_to_sqlite(out_uc)
            collapse_f = _collapse_f_from_sqlite(conn)
            _fasta_from_sqlite(conn, str(sequences), str(matched_seqs))
        except ValueError:
            raise VSearchError('No matches were identified to '
                               'reference_sequences. This can happen if '
                               'sequences are not homologous to '
                               'reference_sequences, or if sequences are '
                               'not in the same orientation as reference_'
                               'sequences (i.e., if sequences are reverse '
                               'complemented with respect to reference '
                               'sequences). Sequence orientation can be '
                               'adjusted with the strand parameter.')

        unmatched_ids = [e.metadata['id']
                         for e in skbio.io.read(open(str(unmatched_seqs)),
                                                constructor=skbio.DNA,
                                                format='fasta')]
    table.filter(ids_to_keep=unmatched_ids, invert=True, axis='observation',
                 inplace=True)
    table = table.collapse(collapse_f, norm=False, min_group_size=1,
                           axis='observation',
                           include_collapsed_metadata=False)

    return table, matched_seqs, unmatched_seqs


def cluster_features_open_reference(ctx, sequences, table, reference_sequences,
                                    perc_identity, strand='plus', threads=1):

    cluster_features_closed_reference = ctx.get_action(
        'vsearch', 'cluster_features_closed_reference')
    filter_features = ctx.get_action('feature_table', 'filter_features')
    cluster_features_de_novo = ctx.get_action(
        'vsearch', 'cluster_features_de_novo')
    merge = ctx.get_action('feature_table', 'merge')
    merge_seqs = ctx.get_action('feature_table', 'merge_seqs')

    skipped_closed_ref = True
    try:
        closed_ref_table, rep_seqs, unmatched_seqs = \
            cluster_features_closed_reference(
                sequences=sequences, table=table,
                reference_sequences=reference_sequences,
                perc_identity=perc_identity,
                strand=strand, threads=threads)
        skipped_closed_ref = False
    except VSearchError:  # No matches
        pass

    # If cluster_features_closed_reference fails to match, we need to
    # pass the source data into cluster_features_de_novo wholesale.
    if skipped_closed_ref:
        unmatched_seqs, closed_ref_table = sequences, table

    # It is possible that all of the sequences matched the reference database,
    # if that is the case, don't worry about running cluster_features_de_novo.
    if unmatched_seqs.view(pd.Series).size > 0:
        unmatched_seqs_md = unmatched_seqs.view(Metadata)
        unmatched_table, = filter_features(table=table,
                                           metadata=unmatched_seqs_md)

        de_novo_table, de_novo_seqs = cluster_features_de_novo(
            sequences=unmatched_seqs, table=unmatched_table,
            perc_identity=perc_identity, threads=threads)

        if skipped_closed_ref:
            merged_reference_seqs, = merge_seqs(data=[reference_sequences,
                                                      de_novo_seqs])
            outputs = (de_novo_table, de_novo_seqs, merged_reference_seqs)
        else:
            merged_table, = merge(
                tables=[closed_ref_table, de_novo_table],
                overlap_method='error_on_overlapping_feature')

            merged_rep_seqs, = merge_seqs(data=[rep_seqs, de_novo_seqs])

            merged_reference_seqs, = merge_seqs(data=[reference_sequences,
                                                      de_novo_seqs])
            outputs = (merged_table, merged_rep_seqs, merged_reference_seqs)
    else:  # skipped de novo
        outputs = (closed_ref_table, rep_seqs, reference_sequences)
    return outputs
