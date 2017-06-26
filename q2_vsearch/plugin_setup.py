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
        'represenative_seqs': FeatureData[Sequence]},
    parameters={
        'id': qiime2.plugin.Float % qiime2.plugin.Range(0, 1, inclusive_start=False, inclusive_end=True)},
    outputs=[
        ('clustered_table', FeatureTable[Frequency]),
        ('clustered_represenative_seqs', FeatureData[Sequence])
    ],
    name='Clusters features at user-specified percent identity.',
    description=('Given a feature table and the associated representative '
                 'sequences, cluster the features based on user-specified '
                 'percent identity threshold of their sequences. This is '
                 'intended to be used for exploring the impact of clustering '
                 'on the results of quality-filtering/dereplication methods, '
                 'such as DADA2. It is not a general-purpose de novo OTU '
                 'clustering method. The output feature ids will be the ids '
                 'of the input features that become cluster centroids.')
)
