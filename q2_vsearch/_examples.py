# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugins import ArtifactAPIUsage

rep_seqs_url = \
    ('https://s3-us-west-2.amazonaws.com/qiime2-data/2024.2/'
     'tutorials/chimera/atacama-rep-seqs.qza')

table_url = \
    ('https://s3-us-west-2.amazonaws.com/qiime2-data/2024.2/'
     'tutorials/chimera/atacama-table.qza')

use = ArtifactAPIUsage()


def cluster_features_de_novo(use):
    rep_seqs = use.init_artifact_from_url('seqs1', rep_seqs_url)
    table = use.init_artifact_from_url('table1', table_url)

    perc_identity = 0.97
    strand = 'plus'
    threads = 1

    clustered_table, clustered_sequences = use.action(
        use.UsageAction('vsearch', 'cluster_features_de_novo'),
        use.UsageInputs(sequences=rep_seqs, table=table,
                        perc_identity=perc_identity, strand=strand,
                        threads=threads),
        use.UsageOutputNames(clustered_table='clustered_table',
                             clustered_sequences='clustered_sequences')
    )

    clustered_table.assert_output_type('FeatureTable[Frequency]')
    clustered_sequences.assert_output_type('FeatureData[Sequence]')
