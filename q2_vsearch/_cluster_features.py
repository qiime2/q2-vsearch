# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
import subprocess

import biom
import skbio
from q2_types.feature_data import DNAFASTAFormat


def run_command(cmd, verbose=True):
    print("Running external command line application. This may print "
          "messages to stdout and/or stderr.")
    print("The command being run is below. This command cannot "
          "be manually re-run as it will depend on temporary files that "
          "no longer exist.")
    print("\nCommand:", end=' ')
    print(" ".join(cmd), end='\n\n')
    subprocess.run(cmd, check=True)


def _collapse_f_from_uc(uc):
    id_to_centroid = {}
    for line in uc:
        line = line.strip()
        if len(line) == 0 or line.startswith(b'#'):
            continue
        else:
            fields = line.split(b'\t')
            if fields[0] == b'S':
                sequence_id = fields[8].decode('utf-8').split(';')[0]
                centroid_id = sequence_id
                id_to_centroid[sequence_id] = centroid_id
            elif fields[0] == b'H':
                sequence_id = fields[8].decode('utf-8').split(';')[0]
                centroid_id = fields[9].decode('utf-8').split(';')[0]
                id_to_centroid[sequence_id] = centroid_id
            else:
                pass

    if len(id_to_centroid) == 0:
        raise ValueError("No sequence matches were identified by vsearch.")

    def collapse_f(id_, x):
        return id_to_centroid[id_]

    return collapse_f


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

    table_ids = set(table_ids)
    non_overlapping_ids = table_ids - sequence_ids
    if len(non_overlapping_ids) != 0:
        raise ValueError('Some feature ids are present in table, but not in '
                         'sequences. The set of features in sequences must be '
                         'identical to the set of features in table. Feature '
                         'ids present in table but not sequences are: %s'
                         % ', '.join(non_overlapping_ids))


def cluster_features_de_novo(sequences: DNAFASTAFormat, table: biom.Table,
                             perc_identity: float
                             )-> (biom.Table, DNAFASTAFormat):
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
                   '--xsize']
            run_command(cmd)
            out_uc.seek(0)
            collapse_f = _collapse_f_from_uc(out_uc)

    table = table.collapse(collapse_f, norm=False, min_group_size=1,
                           axis='observation',
                           include_collapsed_metadata=False)

    return table, clustered_sequences

def cluster_features_closed_reference(sequences: DNAFASTAFormat,
                                      table: biom.Table,
                                      reference_sequences: DNAFASTAFormat,
                                      perc_identity: float,
                                      strand:str ='plus',
                                      # cores?
                                      )-> biom.Table:
    with tempfile.NamedTemporaryFile() as out_uc:
        with tempfile.NamedTemporaryFile() as notmatched:
            cmd = ['vsearch',
                   '--usearch_global', str(sequences),
                   '--id', str(perc_identity),
                   '--db', str(reference_sequences),
                   '--uc', out_uc.name,
                   '--strand', str(strand),
                   '--qmask', 'none',  # ensures no lowercase DNA chars
                   '--notmatched', notmatched.name]
            run_command(cmd)

            out_uc.seek(0)
            try:
                collapse_f = _collapse_f_from_uc(out_uc)
            except ValueError:
                raise ValueError('No matches were identified to '
                                 'reference_sequences. This can happen if '
                                 'sequences are not homologous to '
                                 'reference_sequences, or if sequences are '
                                 'not in the same orientation as reference_'
                                 'sequences (i.e., if sequences are reverse '
                                 'complemented with respect to reference '
                                 'sequences). Sequence orientation can be '
                                 'adjusted with the strand parameter.')

            notmatched.seek(0)
            notmatched_ids =  [e.metadata['id']
                               for e in skbio.io.read(notmatched.name,
                                                      constructor=skbio.DNA,
                                                      format='fasta')]
    table.filter(ids_to_keep=notmatched_ids, invert=True, axis='observation',
                 inplace=True)
    table = table.collapse(collapse_f, norm=False, min_group_size=1,
                           axis='observation',
                           include_collapsed_metadata=False)

    return table
