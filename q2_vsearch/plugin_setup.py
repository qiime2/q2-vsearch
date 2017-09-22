# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin

import q2_vsearch
import q2_vsearch._cluster_features
from q2_types.feature_data import FeatureData, Sequence
from q2_types.feature_table import FeatureTable, Frequency

plugin = qiime2.plugin.Plugin(
    name='vsearch',
    version=q2_vsearch.__version__,
    website='https://github.com/qiime2/q2-vsearch',
    package='q2_vsearch',
    user_support_text=None,
    citation_text=("Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) "
                   "VSEARCH: a versatile open source tool for metagenomics. "
                   "PeerJ 4:e2584. doi: 10.7717/peerj.2584")
)

plugin.methods.register_function(
    function=q2_vsearch._cluster_features.cluster_features_de_novo,
    inputs={
        'table': FeatureTable[Frequency],
        'sequences': FeatureData[Sequence]},
    parameters={
        'perc_identity': qiime2.plugin.Float % qiime2.plugin.Range(
                          0, 1, inclusive_start=False, inclusive_end=True)},
    outputs=[
        ('clustered_table', FeatureTable[Frequency]),
        ('clustered_sequences', FeatureData[Sequence]),
    ],
    input_descriptions={
        'table': 'The feature table to be clustered.',
        'sequences': 'The sequences corresponding to the features in table.',
    },
    parameter_descriptions={
        'perc_identity': ('The percent identity at which clustering should be '
                          'performed. This parameter maps to vsearch\'s --id '
                          'parameter.'),
    },
    output_descriptions={
        'clustered_table': 'The table following clustering of features.',
        'clustered_sequences': 'Sequences representing clustered features.',
    },
    name='De novo clustering of features.',
    description=('Given a feature table and the associated feature '
                 'sequences, cluster the features based on user-specified '
                 'percent identity threshold of their sequences. This is not '
                 'a general-purpose de novo clustering method, but rather is '
                 'intended to be used for clustering the results of '
                 'quality-filtering/dereplication methods, such as DADA2, or '
                 'for re-clustering a FeatureTable at a lower percent '
                 'identity than it was originally clustered at. When a group '
                 'of features in the input table are clustered into a single '
                 'feature, the frequency of that single feature in a given '
                 'sample is the sum of the frequencies of the features that '
                 'were clustered in that sample. Feature identifiers and '
                 'sequences will be inherited from the centroid feature '
                 'of each cluster. See the vsearch documentation for details '
                 'on how sequence clustering is performed.')
)

plugin.methods.register_function(
    function=q2_vsearch._cluster_features.cluster_features_closed_reference,
    inputs={
        'table': FeatureTable[Frequency],
        'sequences': FeatureData[Sequence],
        'reference_sequences': FeatureData[Sequence]
    },
    parameters={
        'perc_identity': qiime2.plugin.Float % qiime2.plugin.Range(
                          0, 1, inclusive_start=False, inclusive_end=True),
        'strand': qiime2.plugin.Str % qiime2.plugin.Choices(['plus', 'both']),
    },
    outputs=[
        ('clustered_table', FeatureTable[Frequency]),
    ],
    input_descriptions={
        'table': 'The feature table to be clustered.',
        'sequences': 'The sequences corresponding to the features in table.',
        'reference_sequences': 'The sequences to use as cluster centroids.',
    },
    parameter_descriptions={
        'perc_identity': ('The percent identity at which clustering should be '
                          'performed. This parameter maps to vsearch\'s --id '
                          'parameter.'),
        'strand': ('Search plus (i.e., forward) or both (i.e., forward and '
                   'reverse complement) strands.')
    },
    output_descriptions={
        'clustered_table': 'The table following clustering of features.',
    },
    name='Closed-reference clustering of features.',
    description=('Given a feature table and the associated feature '
                 'sequences, cluster the features against a reference '
                 'database based on user-specified '
                 'percent identity threshold of their sequences. This is not '
                 'a general-purpose closed-reference clustering method, but rather is '
                 'intended to be used for clustering the results of '
                 'quality-filtering/dereplication methods, such as DADA2, or '
                 'for re-clustering a FeatureTable at a lower percent '
                 'identity than it was originally clustered at. When a group '
                 'of features in the input table are clustered into a single '
                 'feature, the frequency of that single feature in a given '
                 'sample is the sum of the frequencies of the features that '
                 'were clustered in that sample. Feature identifiers '
                 'will be inherited from the centroid feature '
                 'of each cluster. See the vsearch documentation for details '
                 'on how sequence clustering is performed.')
)
