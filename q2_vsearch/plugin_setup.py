import qiime2.plugin

import q2_vsearch
import q2_vsearch._cluster_features
from q2_types.feature_data import DNAFASTAFormat, FeatureData, Sequence
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
    function=q2_vsearch._cluster_features.cluster_features,
    inputs={
        'table': FeatureTable[Frequency],
        'sequences': FeatureData[Sequence]},
    parameters={
        'id': qiime2.plugin.Float % qiime2.plugin.Range(
            0, 1, inclusive_start=False, inclusive_end=True)},
    outputs=[
        ('clustered_table', FeatureTable[Frequency]),
        ('clustered_sequences', FeatureData[Sequence]),
    ],
    input_descriptions = {
        'table': 'The feature table to be clustered.',
        'sequences': 'The sequences corresponding to the features in table.',
    },
    parameter_descriptions = {
        'id': 'The percent identity at which clustering should be performed.',
    },
    output_descriptions = {
        'clustered_table': 'The table following clustering of features.',
        'clustered_sequences': 'Sequences representing clustered features.',
    },
    name='Cluster features at user-specified percent identity.',
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
                 'on how sequence clustering is performed.'
                )
)
