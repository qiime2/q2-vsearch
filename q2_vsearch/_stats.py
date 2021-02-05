# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import fileinput
import pkg_resources
import subprocess
import pandas as pd
from multiprocessing import Pool, cpu_count

from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt
)
import q2templates

TEMPLATES = pkg_resources.resource_filename('q2_vsearch', 'assets')


def _get_stats_easy(cmds_packed) -> None:
    filelist, cmds = cmds_packed
    processes = []
    for cmd in cmds:
        processes.append(subprocess.Popen(cmd, stdin=subprocess.PIPE))

    with fileinput.input(files=filelist, mode='r',
                         openhook=fileinput.hook_compressed) as fh:
        for line in fh:
            for p in processes:
                p.stdin.write(line)

    for p in processes:
        p.stdin.close()
        p.wait()


def _build_cmds(output_dir: str, filelist, direction='forward'):
    datafiles = {direction: {}}

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

    return (datafiles, [(filelist, [stats, eestats, eestats2])])


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
    return html


def _fastq_stats(output_dir: str, sequences, threads) -> None:
    # read manifest
    manifest = sequences.manifest
    # check if paired reads available
    try:
        paired = manifest['reverse'][0] is not None
    except KeyError:
        paired = False

    # get commands and filelist
    datafiles, cmds = _build_cmds(output_dir, manifest['forward'].tolist())

    if (paired):
        datafiles_rev, cmds_rev = _build_cmds(output_dir,
                                              manifest['reverse'].tolist(),
                                              'reverse')
        datafiles.update(datafiles_rev)
        cmds.extend(cmds_rev)

    # multiprocessing
    cpus = cpu_count()
    try:
        if (cpus < threads):
            threads = cpus
    except TypeError:
        # QIIME itself checks for allowed input format and values
        # (see plugin_setup.py)
        threads = cpus

    if (threads < 4):  # read once and write once (3 or 6 times)
        # refactor into (filelist, [single_command]) to only spawn one
        # additional process in _get_stats_easy
        for cmd_packed in cmds:
            filelist, cmds = cmd_packed
            for cmd in cmds:
                _get_stats_easy((filelist, [cmd]))
    else:   # read once and write three times (1 or 2 times)
        # three additional processes spawn in _get_stats_easy
        # so only one (single) or two (paired) worker are needed
        jobs = 1
        if (paired and threads >= 8):   # parallel fwd/rev
            jobs = 2
        with Pool(processes=jobs) as pool:
            pool.map(_get_stats_easy, cmds)
            pool.close()

    html = _get_html(output_dir, datafiles)

    index = os.path.join(TEMPLATES, 'index.html')
    eestats_template = os.path.join(TEMPLATES, 'fastq_eestats.html')
    eestats2_template = os.path.join(TEMPLATES, 'fastq_eestats2.html')
    context = {
        'paired': paired,
        'datafiles': datafiles,
        'html': html,
        'tabs': [{'title': 'fastq_stats',
                  'url': 'index.html'},
                 {'title': 'fastq_eestats',
                  'url': 'fastq_eestats.html'},
                 {'title': 'fastq_eestats2',
                  'url': 'fastq_eestats2.html'}],
        }
    templates = [index, eestats_template, eestats2_template]
    q2templates.render(templates, output_dir, context=context)


def fastq_stats(output_dir: str,
                sequences: CasavaOneEightSingleLanePerSampleDirFmt,
                threads: int = 1
                ) -> None:
    _fastq_stats(output_dir, sequences, threads)
