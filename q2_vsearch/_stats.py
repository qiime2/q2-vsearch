import os
# import tempfile
# import hashlib
import subprocess

# import biom
# import skbio
# import qiime2.util
import pandas as pd

# from q2_types.feature_data import DNAIterator
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)

from ._cluster_features import run_command


def _build_commands(output_dir: str, filelist, direction='forward'):
    cmd = ['zcat'] + filelist

    results = os.path.join(
                output_dir, 'fastq_stats_{0}.txt'.format(direction))
    stats = ['vsearch', '--fastq_stats', '-', '--log', results]

    results = os.path.join(
                output_dir, 'fastq_eestats_{0}.txt'.format(direction))
    eestats = ['vsearch', '--fastq_eestats', '-', '--output', results]

    results = os.path.join(
                output_dir, 'fastq_eestats2_{0}.txt'.format(direction))
    eestats2 = ['vsearch', '--fastq_eestats2', '-', '--output', results]

    for cmd2 in [stats, eestats, eestats2]:
        # pipe it!
        # manual console 4MB, here 160MB RAM!!!
        process1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        process2 = subprocess.run(cmd2, stdin=process1.stdout)


# for pairedEnd only
def fastq_stats(output_dir: str,
                sequences: SingleLanePerSamplePairedEndFastqDirFmt) -> None:
    # SingleLanePerSampleSingleEndFastqDirFmt
    # read manifest and add complete path
    manifest = pd.read_csv(os.path.join(str(sequences),
                           sequences.manifest.pathspec),
                           header=0, comment='#')
    manifest.filename = manifest.filename.apply(
        lambda x: os.path.join(str(sequences), x))

    # filter read direction
    r_fwd = manifest[manifest.direction == 'forward'].filename.values.tolist()
    r_rev = manifest[manifest.direction == 'reverse'].filename.values.tolist()

    _build_commands(output_dir, r_fwd)
    _build_commands(output_dir, r_rev, 'reverse')
