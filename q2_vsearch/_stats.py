import os
import pkg_resources
import subprocess
import pandas as pd

from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)
import q2templates

TEMPLATES = pkg_resources.resource_filename('q2_vsearch', 'assets')


def _build_commands(output_dir: str, filelist, direction='forward'):
    url = {direction: {}}
    cmd = ['zcat'] + filelist

    results = os.path.join(
                output_dir, 'fastq_stats_{0}.txt'.format(direction))
    stats = ['vsearch', '--fastq_stats', '-', '--log', results]
    url[direction]['stats'] = os.path.basename(results)

    results = os.path.join(
                output_dir, 'fastq_eestats_{0}.txt'.format(direction))
    eestats = ['vsearch', '--fastq_eestats', '-', '--output', results]
    url[direction]['eestats'] = os.path.basename(results)

    results = os.path.join(
                output_dir, 'fastq_eestats2_{0}.txt'.format(direction))
    eestats2 = ['vsearch', '--fastq_eestats2', '-', '--output', results]
    url[direction]['eestats2'] = os.path.basename(results)

    for cmd2 in [stats, eestats, eestats2]:
        # pipe
        process1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        process2 = subprocess.run(cmd2, stdin=process1.stdout)
    return(url)


def _get_html(datafiles):
    html = {}
    for direction in datafiles:
        html[direction] = {}
        for stats_type in datafiles[direction]:
            filename = datafiles[direction][stats_type]
            data_df = pd.read_csv(filename, sep='\t')
            html[direction][stats_type] = q2templates.df_to_html(data_df,
                                                                 index=False)
    return(html)


def _fastq_stats(output_dir: str, sequences, paired=False):
    url = {}
    # read manifest and add complete path
    manifest = pd.read_csv(os.path.join(str(sequences),
                           sequences.manifest.pathspec),
                           header=0, comment='#')
    manifest.filename = manifest.filename.apply(
        lambda x: os.path.join(str(sequences), x))

    # filter read direction
    r_fwd = manifest[manifest.direction == 'forward'].filename.values.tolist()
    datafiles = _build_commands(output_dir, r_fwd)
    if (paired):
        r_rev = manifest[
            manifest.direction == 'reverse'].filename.values.tolist()
        datafiles.update(_build_commands(output_dir, r_rev, 'reverse'))

    html = _get_html(datafiles)

    index = os.path.join(TEMPLATES, 'index.html')
    stats_template = os.path.join(TEMPLATES, 'fastq_stats.html')
    eestats_template = os.path.join(TEMPLATES, 'fastq_eestats.html')
    eestats2_template = os.path.join(TEMPLATES, 'fastq_eestats2.html')
    context = {
        'paired': paired,
        'datafiles': datafiles,
        'html': html,
        'tabs': [{'title': 'fastq_stats',
                  'url': 'fastq_stats.html'},
                 {'title': 'fastq_eestats',
                  'url': 'fastq_eestats.html'},
                 {'title': 'fastq_eestats2',
                  'url': 'fastq_eestats2.html'}],
    }
    templates = [index, stats_template, eestats_template, eestats2_template]
    q2templates.render(templates, output_dir, context=context)


def fastq_stats_paired(output_dir: str,
                       sequences: SingleLanePerSamplePairedEndFastqDirFmt
                       ) -> None:
    _fastq_stats(output_dir, sequences, True)


def fastq_stats_single(output_dir: str,
                       sequences: SingleLanePerSampleSingleEndFastqDirFmt
                       ) -> None:
    _fastq_stats(output_dir, sequences)
