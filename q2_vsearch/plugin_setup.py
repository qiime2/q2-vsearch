# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

import qiime2.plugin
import q2_vsearch._cluster_features
import q2_vsearch._cluster_sequences
import q2_vsearch._join_pairs
import q2_vsearch._chimera
import q2_vsearch._stats

from q2_vsearch._type import UchimeStats
from q2_vsearch._format import UchimeStatsFmt, UchimeStatsDirFmt
from q2_types.feature_data import FeatureData, Sequence
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (
    Sequences, SequencesWithQuality, PairedEndSequencesWithQuality,
    JoinedSequencesWithQuality)

citations = qiime2.plugin.Citations.load('citations.bib', package='q2_vsearch')
plugin = qiime2.plugin.Plugin(
    name='vsearch',
    version=q2_vsearch.__version__,
    website='https://github.com/qiime2/q2-vsearch',
    package='q2_vsearch',
    user_support_text=None,
    short_description='Plugin for clustering and dereplicating with vsearch.',
    description=('This plugin wraps the vsearch application, and provides '
                 'methods for clustering and dereplicating features and '
                 'sequences.'),
    citations=[citations['rognes2016vsearch']]
)

plugin.register_formats(UchimeStatsFmt, UchimeStatsDirFmt)

plugin.register_semantic_types(UchimeStats)
plugin.register_semantic_type_to_format(
    UchimeStats,
    artifact_format=UchimeStatsDirFmt)

plugin.methods.register_function(
    function=q2_vsearch._cluster_features.cluster_features_de_novo,
    inputs={
        'table': FeatureTable[Frequency],
        'sequences': FeatureData[Sequence]},
    parameters={
        'perc_identity': qiime2.plugin.Float % qiime2.plugin.Range(
                          0, 1, inclusive_start=False, inclusive_end=True),
        'threads': qiime2.plugin.Int % qiime2.plugin.Range(
                          0, 256, inclusive_start=True, inclusive_end=True)
    },
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
        'threads': ('The number of threads to use for computation. Passing 0 '
                    'will launch one thread per CPU core.')
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
        'threads': qiime2.plugin.Int % qiime2.plugin.Range(
                0, 256, inclusive_start=True, inclusive_end=True)
    },
    outputs=[
        ('clustered_table', FeatureTable[Frequency]),
        ('clustered_sequences', FeatureData[Sequence]),
        ('unmatched_sequences', FeatureData[Sequence]),
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
                   'reverse complement) strands.'),
        'threads': ('The number of threads to use for computation. Passing 0 '
                    'will launch one thread per CPU core.')
    },
    output_descriptions={
        'clustered_table': 'The table following clustering of features.',
        'clustered_sequences': 'The sequences representing clustered '
                               'features, relabeled by the reference IDs.',
        'unmatched_sequences': 'The sequences which failed to match any '
                               'reference sequences. This output maps to '
                               'vsearch\'s --notmatched parameter.'
    },
    name='Closed-reference clustering of features.',
    description=('Given a feature table and the associated feature '
                 'sequences, cluster the features against a reference '
                 'database based on user-specified '
                 'percent identity threshold of their sequences. This is not '
                 'a general-purpose closed-reference clustering method, but '
                 'rather is intended to be used for clustering the results of '
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

plugin.pipelines.register_function(
    function=q2_vsearch._cluster_features.cluster_features_open_reference,
    inputs={
        'table': FeatureTable[Frequency],
        'sequences': FeatureData[Sequence],
        'reference_sequences': FeatureData[Sequence]
    },
    parameters={
        'perc_identity': qiime2.plugin.Float % qiime2.plugin.Range(
                          0, 1, inclusive_start=False, inclusive_end=True),
        'strand': qiime2.plugin.Str % qiime2.plugin.Choices(['plus', 'both']),
        'threads': qiime2.plugin.Int % qiime2.plugin.Range(
                0, 256, inclusive_start=True, inclusive_end=True)
    },
    outputs=[
        ('clustered_table', FeatureTable[Frequency]),
        ('clustered_sequences', FeatureData[Sequence]),
        ('new_reference_sequences', FeatureData[Sequence]),
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
                   'reverse complement) strands.'),
        'threads': ('The number of threads to use for computation. Passing 0 '
                    'will launch one thread per CPU core.')
    },
    output_descriptions={
        'clustered_table': 'The table following clustering of features.',
        'clustered_sequences': 'Sequences representing clustered features.',
        'new_reference_sequences': 'The new reference sequences. This can be '
                                   'used for subsequent runs of '
                                   'open-reference clustering for consistent '
                                   'definitions of features across '
                                   'open-reference feature tables.',
    },
    name='Open-reference clustering of features.',
    description='Given a feature table and the associated feature sequences, '
                'cluster the features against a reference database based on '
                'user-specified percent identity threshold of their sequences.'
                ' Any sequences that don\'t match are then clustered de novo. '
                'This is not a general-purpose clustering method, but rather '
                'is intended to be used for clustering the results of '
                'quality-filtering/dereplication methods, such as DADA2, or '
                'for re-clustering a FeatureTable at a lower percent identity '
                'than it was originally clustered at. When a group of '
                'features in the input table are clustered into a single '
                'feature, the frequency of that single feature in a given '
                'sample is the sum of the frequencies of the features that '
                'were clustered in that sample. Feature identifiers will be '
                'inherited from the centroid feature of each cluster. For '
                'features that match a reference sequence, the centroid '
                'feature is that reference sequence, so its identifier will '
                'become the feature identifier. The clustered_sequences '
                'result will contain feature representative sequences that '
                'are derived from the sequences input for all features in '
                'clustered_table. This will always be the most abundant '
                'sequence in the cluster. The new_reference_sequences result '
                'will contain the entire reference database, plus feature '
                'representative sequences for any de novo features. This is '
                'intended to be used as a reference database in subsequent '
                'iterations of cluster_features_open_reference, if '
                'applicable. See the vsearch documentation for details on how '
                'sequence clustering is performed.',
    citations=[citations['rideout2014subsampled']]
)

plugin.methods.register_function(
    function=q2_vsearch._cluster_sequences.dereplicate_sequences,
    inputs={
        'sequences': (SampleData[Sequences] |
                      SampleData[SequencesWithQuality] |
                      SampleData[JoinedSequencesWithQuality])
    },
    parameters={
        'derep_prefix': qiime2.plugin.Bool,
    },
    outputs=[
        ('dereplicated_table', FeatureTable[Frequency]),
        ('dereplicated_sequences', FeatureData[Sequence]),
    ],
    input_descriptions={
        'sequences': 'The sequences to be dereplicated.',
    },
    parameter_descriptions={
        'derep_prefix': ('Merge sequences with identical prefixes. If a '
                         'sequence is identical to the prefix of two or more '
                         'longer sequences, it is clustered with the shortest '
                         'of them. If they are equally long, it is clustered '
                         'with the most abundant.'),
    },
    output_descriptions={
        'dereplicated_table': 'The table of dereplicated sequences.',
        'dereplicated_sequences': 'The dereplicated sequences.',
    },
    name='Dereplicate sequences.',
    description=('Dereplicate sequence data and create a feature table and '
                 'feature representative sequences. Feature identifiers '
                 'in the resulting artifacts will be the sha1 hash '
                 'of the sequence defining each feature. If clustering of '
                 'features into OTUs is desired, the resulting artifacts '
                 'can be passed to the cluster_features_* methods in this '
                 'plugin.')
)

plugin.methods.register_function(
    function=q2_vsearch._join_pairs.join_pairs,
    inputs={
        'demultiplexed_seqs': SampleData[PairedEndSequencesWithQuality]
    },
    parameters={
        'truncqual': qiime2.plugin.Int % qiime2.plugin.Range(0, None),
        'minlen': qiime2.plugin.Int % qiime2.plugin.Range(0, None),
        'maxns': qiime2.plugin.Int % qiime2.plugin.Range(0, None),
        'allowmergestagger': qiime2.plugin.Bool,
        'minovlen': qiime2.plugin.Int % qiime2.plugin.Range(0, None),
        'maxdiffs': qiime2.plugin.Int % qiime2.plugin.Range(0, None),
        'minmergelen': qiime2.plugin.Int % qiime2.plugin.Range(0, None),
        'maxmergelen': qiime2.plugin.Int % qiime2.plugin.Range(0, None),
        'maxee': qiime2.plugin.Float % qiime2.plugin.Range(0., None),
        'qmin': qiime2.plugin.Int % qiime2.plugin.Range(
            -5, 2, inclusive_start=True, inclusive_end=True),
        'qminout': qiime2.plugin.Int % qiime2.plugin.Range(
            -5, 2, inclusive_start=True, inclusive_end=True),
        'qmax': qiime2.plugin.Int % qiime2.plugin.Range(
            40, 41, inclusive_start=True, inclusive_end=True),
        'qmaxout': qiime2.plugin.Int % qiime2.plugin.Range(
            40, 41, inclusive_start=True, inclusive_end=True),
        'threads': qiime2.plugin.Int % qiime2.plugin.Range(
            0, 8, inclusive_start=True, inclusive_end=True)
    },
    outputs=[
        ('joined_sequences', SampleData[JoinedSequencesWithQuality])
    ],
    input_descriptions={
        'demultiplexed_seqs': ('The demultiplexed paired-end sequences to '
                               'be joined.'),
    },
    parameter_descriptions={
        'truncqual': ('Truncate sequences at the first base with the '
                      'specified quality score value or lower.'),
        'minlen': ('Sequences shorter than minlen after truncation are '
                   'discarded.'),
        'maxns': ('Sequences with more than maxns N characters are '
                  'discarded.'),
        'allowmergestagger': ('Allow joining of staggered read pairs.'),
        'minovlen': ('Minimum overlap length of forward and reverse reads '
                     'for joining.'),
        'maxdiffs': ('Maximum number of mismatches in the forward/reverse '
                     'read overlap for joining.'),
        'minmergelen': ('Minimum length of the joined read to be retained.'),
        'maxmergelen': ('Maximum length of the joined read to be retained.'),
        'maxee': ('Maximum number of expected errors in the joined read '
                  'to be retained.'),
        'qmin': ('The minimum allowed quality score in the input.'),
        'qminout': ('The minimum allowed quality score to use in output.'),
        'qmax': ('The maximum allowed quality score in the input.'),
        'qmaxout': ('The maximum allowed quality score to use in output.'),
        'threads': ('The number of threads to use for computation. Does '
                    'not scale much past 4 threads.')
    },
    output_descriptions={
        'joined_sequences': ('The joined sequences.'),
    },
    name='Join paired-end reads.',
    description=('Join paired-end sequence reads using vsearch\'s '
                 'merge_pairs function. The qmin, qminout, qmax, and qmaxout '
                 'parameters should only need to be modified when working '
                 'with older fastq sequence data. See the vsearch '
                 'documentation for details on how paired-end joining is '
                 'performed, and for more information on the parameters to '
                 'this method.')
)

plugin.methods.register_function(
    function=q2_vsearch._chimera.uchime_ref,
    inputs={
        'sequences': FeatureData[Sequence],
        'table': FeatureTable[Frequency],
        'reference_sequences': FeatureData[Sequence]},
    parameters={
        'dn': qiime2.plugin.Float % qiime2.plugin.Range(0., None),
        'mindiffs': qiime2.plugin.Int % qiime2.plugin.Range(1, None),
        'mindiv': qiime2.plugin.Float % qiime2.plugin.Range(0., None),
        'minh': qiime2.plugin.Float % qiime2.plugin.Range(
                          0., 1.0, inclusive_end=True),
        'xn': qiime2.plugin.Float % qiime2.plugin.Range(
                          1., None, inclusive_start=False),
        'threads': qiime2.plugin.Int % qiime2.plugin.Range(
                          0, 256, inclusive_start=True, inclusive_end=True)
    },
    outputs=[
        ('chimeras', FeatureData[Sequence]),
        ('nonchimeras', FeatureData[Sequence]),
        ('stats', UchimeStats)
    ],
    input_descriptions={
        'sequences': 'The feature sequences to be chimera-checked.',
        'table': ('Feature table (used for computing total feature '
                  'abundances).'),
        'reference_sequences': 'The non-chimeric reference sequences.'
    },
    parameter_descriptions={
        'dn': ('No vote pseudo-count, corresponding to the parameter n in '
               'the chimera scoring function.'),
        'mindiffs': 'Minimum number of differences per segment.',
        'mindiv': 'Minimum divergence from closest parent.',
        'minh': ('Minimum score (h). Increasing this value tends to reduce '
                 'the number of false positives and to decrease sensitivity.'),
        'xn': ('No vote weight, corresponding to the parameter beta in the '
               'scoring function.'),
        'threads': ('The number of threads to use for computation. Passing 0 '
                    'will launch one thread per CPU core.')
    },
    output_descriptions={
        'chimeras': 'The chimeric sequences.',
        'nonchimeras': 'The non-chimeric sequences.',
        'stats': 'Summary statistics from chimera checking.'
    },
    name='Reference-based chimera filtering with vsearch.',
    description=('Apply the vsearch uchime_ref method to identify chimeric '
                 'feature sequences. The results of this method can be used '
                 'to filter chimeric features from the corresponding feature '
                 'table. For additional details, please refer to the vsearch '
                 'documentation.')
)

plugin.methods.register_function(
    function=q2_vsearch._chimera.uchime_denovo,
    inputs={
        'sequences': FeatureData[Sequence],
        'table': FeatureTable[Frequency]},
    parameters={
        'dn': qiime2.plugin.Float % qiime2.plugin.Range(0., None),
        'mindiffs': qiime2.plugin.Int % qiime2.plugin.Range(1, None),
        'mindiv': qiime2.plugin.Float % qiime2.plugin.Range(0., None),
        'minh': qiime2.plugin.Float % qiime2.plugin.Range(
                          0., 1.0, inclusive_end=True),
        'xn': qiime2.plugin.Float % qiime2.plugin.Range(
                          1., None, inclusive_start=False)
    },
    outputs=[
        ('chimeras', FeatureData[Sequence]),
        ('nonchimeras', FeatureData[Sequence]),
        ('stats', UchimeStats)
    ],
    input_descriptions={
        'sequences': 'The feature sequences to be chimera-checked.',
        'table': ('Feature table (used for computing total feature '
                  'abundances).'),
    },
    parameter_descriptions={
        'dn': ('No vote pseudo-count, corresponding to the parameter n in '
               'the chimera scoring function.'),
        'mindiffs': 'Minimum number of differences per segment.',
        'mindiv': 'Minimum divergence from closest parent.',
        'minh': ('Minimum score (h). Increasing this value tends to reduce '
                 'the number of false positives and to decrease sensitivity.'),
        'xn': ('No vote weight, corresponding to the parameter beta in the '
               'scoring function.'),
    },
    output_descriptions={
        'chimeras': 'The chimeric sequences.',
        'nonchimeras': 'The non-chimeric sequences.',
        'stats': 'Summary statistics from chimera checking.'
    },
    name='De novo chimera filtering with vsearch.',
    description=('Apply the vsearch uchime_denovo method to identify chimeric '
                 'feature sequences. The results of this method can be used '
                 'to filter chimeric features from the corresponding feature '
                 'table. For additional details, please refer to the vsearch '
                 'documentation.')
)


plugin.visualizers.register_function(
    function=q2_vsearch._stats.fastq_stats,
    inputs={
        'sequences': SampleData[
            SequencesWithQuality | PairedEndSequencesWithQuality],
    },
    parameters={
        'threads': qiime2.plugin.Int % qiime2.plugin.Range(
            1, None) | qiime2.plugin.Str % qiime2.plugin.Choices(['auto'])
    },
    input_descriptions={
        'sequences': 'Fastq sequences'
    },
    parameter_descriptions={
        'threads': 'The number of threads used for computation.',
    },
    name='Fastq stats with vsearch.',
    description='A fastq overview via vsearch\'s fastq_stats, fastq_eestats '
                'and fastq_eestats2 utilities. Please see '
                'https://github.com/torognes/vsearch for detailed '
                'documentation of these tools.',
)

importlib.import_module('q2_vsearch._transformer')
