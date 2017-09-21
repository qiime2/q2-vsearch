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
                id_to_centroid[fields[8].decode('utf-8')] = \
                    fields[8].decode('utf-8')
            elif fields[0] == b'H':
                id_to_centroid[fields[8].decode('utf-8')] = \
                    fields[9].decode('utf-8')
            else:
                pass

    def collapse_f(id_, x):
        return id_to_centroid[id_]

    return collapse_f


def cluster_features_denovo(sequences: DNAFASTAFormat, table: biom.Table,
                            perc_identity: float
                            )-> (biom.Table, DNAFASTAFormat):
    sequences_fp = str(sequences)

    feature_ids_seqs = {e.metadata['id'] for e in
                        skbio.io.read(sequences_fp, constructor=skbio.DNA,
                                      format='fasta')}
    feature_ids_table = set(table.ids(axis='observation'))
    non_overlapping_ids = feature_ids_seqs ^ feature_ids_table

    if len(non_overlapping_ids) > 0:
        raise ValueError('All feature ids must be present in table and '
                         'sequences, but some are not. Feature ids not '
                         'present in table and sequences are: '
                         '%s' % ', '.join(non_overlapping_ids))

    clustered_sequences = DNAFASTAFormat()
    with tempfile.NamedTemporaryFile() as out_uc:
        cmd = ['vsearch', '--cluster_fast', sequences_fp, '--id',
               str(perc_identity), '--centroids', str(clustered_sequences),
               '--uc', out_uc.name, '--qmask', 'none']
        run_command(cmd)
        out_uc.seek(0)
        collapse_f = _collapse_f_from_uc(out_uc)

    table = table.collapse(collapse_f, norm=False, min_group_size=1,
                           axis='observation',
                           include_collapsed_metadata=False)

    return table, clustered_sequences
