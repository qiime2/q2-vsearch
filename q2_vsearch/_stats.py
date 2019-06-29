import os
import pkg_resources
import subprocess
import pandas as pd
from multiprocessing import Pool, cpu_count

from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)
import q2templates

TEMPLATES = pkg_resources.resource_filename('q2_vsearch', 'assets')


def _get_stats(cmds) -> None:
    cmd, cmd2 = cmds
    process1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    process2 = subprocess.Popen(cmd2, stdin=process1.stdout)

    while True:
        process1.poll()
        process2.poll()
        if process1.returncode is not None and process2.returncode is not None:
            process1.communicate(timeout=15)
            break


def _build_cmds(output_dir: str, filelist, direction='forward'):
    datafiles = {direction: {}}
    cmd = ['zcat'] + filelist

    results = os.path.join(
                output_dir, 'fastq_stats_{0}.txt'.format(direction))
    stats = ['vsearch', '--quiet', '--fastq_stats', '-', '--log', results]
    datafiles[direction]['stats'] = os.path.basename(results)

    results = os.path.join(
                output_dir, 'fastq_eestats_{0}.txt'.format(direction))
    eestats = ['vsearch', '--quiet', '--fastq_eestats',
               '-', '--output', results]
    datafiles[direction]['eestats'] = os.path.basename(results)

    results = os.path.join(
                output_dir, 'fastq_eestats2_{0}.txt'.format(direction))
    eestats2 = ['vsearch', '--quiet', '--fastq_eestats2',
                '-', '--output', results]
    datafiles[direction]['eestats2'] = os.path.basename(results)

    return(datafiles, [(cmd, stats), (cmd, eestats), (cmd, eestats2)])


def _get_html(output_dir, datafiles):
    html = {}
    for direction in datafiles:
        html[direction] = {}
        for stats_type in datafiles[direction]:
            filename = datafiles[direction][stats_type]
            filename = os.path.join(output_dir, filename)
            data_df = pd.read_csv(filename, sep='\t')
            html[direction][stats_type] = q2templates.df_to_html(data_df,
                                                                 index=False)
    return(html)


def _fastq_stats(output_dir: str, sequences, threads, paired=False) -> None:
    # read manifest and add complete path
    manifest = pd.read_csv(os.path.join(str(sequences),
                                        sequences.manifest.pathspec),
                           header=0, comment='#')
    manifest.filename = manifest.filename.apply(
        lambda x: os.path.join(str(sequences), x))

    # get commands and filelist
    r_fwd = manifest[manifest.direction == 'forward'].filename.values.tolist()
    datafiles, cmds = _build_cmds(output_dir, r_fwd)
    if (paired):
        r_rev = manifest[
            manifest.direction == 'reverse'].filename.values.tolist()
        datafiles2, cmds2 = _build_cmds(output_dir, r_rev, 'reverse')
        cmds.extend(cmds2)
        datafiles.update(datafiles2)

    # multiprocessing
    cpus = cpu_count()
    if ((cpus < threads) or (threads == 0)):
        threads = cpus
    with Pool(processes=threads) as pool:
        pool.map(_get_stats, cmds)
        pool.close()

    html = _get_html(output_dir, datafiles)

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
                       sequences: SingleLanePerSamplePairedEndFastqDirFmt,
                       threads: int = 1
                       ) -> None:
    _fastq_stats(output_dir, sequences, threads, True)


def fastq_stats_single(output_dir: str,
                       sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                       threads: int = 1
                       ) -> None:
    _fastq_stats(output_dir, sequences, threads)
