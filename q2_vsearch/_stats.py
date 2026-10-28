# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import json
import fileinput
import importlib
import subprocess
import pandas as pd
from multiprocessing import Pool, cpu_count

from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt
)
import q2templates

TEMPLATES = importlib.resources.files('q2_vsearch') / 'assets'

# fastq_eestats2 reports reads kept at every length, for these max EE
EE_CUTOFFS = [0.5, 1, 2, 3]
LENGTH_ONLY = 'length only'
# Percent of reads kept at each length by each filter, as plotted
COLUMNS = ['direction', 'Length', 'series', 'PctReads']


def _get_stats_easy(cmds_packed) -> None:
    filelist, cmds = cmds_packed
    processes = []
    for cmd in cmds:
        processes.append(subprocess.Popen(cmd, stdin=subprocess.PIPE))

    with fileinput.input(files=filelist, mode='r',
                         openhook=fileinput.hook_compressed) as fh:
        for line in fh:
            encoded_line = line.encode('utf-8')
            for p in processes:
                p.stdin.write(encoded_line)

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
    eestats2 = ['vsearch', '--quiet', '--fastq_eestats2', '-',
                '--length_cutoffs', '1,*,1',
                '--ee_cutoffs', ','.join(str(ee) for ee in EE_CUTOFFS),
                '--output', results]
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


def _json_replace(json_obj, **values):
    """
    Replace each `{"{{REPLACE_PARAM}}": "some_key"}` in `json_obj` with
    `values["some_key"]`, as in q2-stats.
    """
    if type(json_obj) is dict and list(json_obj) == ["{{REPLACE_PARAM}}"]:
        return values[json_obj["{{REPLACE_PARAM}}"]]
    if type(json_obj) is list:
        return [_json_replace(x, **values) for x in json_obj]
    if type(json_obj) is dict:
        return {key: _json_replace(value, **values)
                for key, value in json_obj.items()}
    return json_obj


def _read_eestats(filename, direction):
    df = pd.read_csv(filename, sep='\t')
    df = df[['Pos', 'PctRecs', 'Low_Q', 'Med_Q', 'Mean_Q', 'Hi_Q',
             'Low_EE', 'Med_EE', 'Mean_EE', 'Hi_EE']]
    df.insert(0, 'direction', direction)
    return df


def _maxee_name(maxee):
    return 'max EE {0:g}'.format(maxee)


def _read_eestats2(filename, direction):
    """
    Parse the fixed-width table of fastq_eestats2 into the percent of reads
    kept at each length and max EE.
    """
    with open(filename) as fh:
        lines = fh.read().splitlines()
    header = next(line for line in lines if line.startswith('Length'))
    maxees = [float(ee) for ee in re.findall(r'MaxEE ([\d.]+)', header)]
    rows = []
    for line in lines:
        length = re.match(r'\s*(\d+)\s', line)
        kept = re.findall(r'\(\s*([\d.]+)%\)', line)
        for maxee, pct in zip(maxees, kept if length else []):
            rows.append({'direction': direction,
                         'Length': int(length.group(1)),
                         'series': _maxee_name(maxee),
                         'PctReads': float(pct)})
    return pd.DataFrame(rows, columns=COLUMNS)


def _section(lines, title):
    """The rows of a table in the fastq_stats log, after its header."""
    start = lines.index(title) + 1
    rows = []
    for line in lines[start:]:
        if not line.strip():
            break
        rows.append(line)
    return rows


def _read_fastq_stats(filename, direction):
    """
    Parse the read length distribution, and the percent of reads kept at
    each length when truncated at the first base at or under a Q score, from
    the fastq_stats log.
    """
    with open(filename) as fh:
        lines = fh.read().splitlines()

    lengths = []
    for line in _section(lines, 'Read length distribution'):
        row = re.match(r'\s*(?:>=)?\s*(\d+)\s+(\d+)\s+([\d.]+)%', line)
        if row:
            lengths.append({'direction': direction,
                            'Length': int(row.group(1)),
                            'Reads': int(row.group(2)),
                            'PctReads': float(row.group(3))})

    rows = _section(lines, 'Truncate at first Q')
    qs = [int(q) for q in re.findall(r'Q=(\d+)', rows[0])]
    truncq = []
    for line in rows[1:]:
        length = re.match(r'\s*(\d+)\s', line)
        if not length:
            continue
        for q, pct in zip(qs, re.findall(r'([\d.]+)%', line)):
            truncq.append({'direction': direction,
                           'Length': int(length.group(1)),
                           'series': 'first Q \u2264 {0}'.format(q),
                           'PctReads': float(pct)})
    lengths = pd.DataFrame(
        lengths, columns=['direction', 'Length', 'Reads', 'PctReads'])
    return lengths, pd.DataFrame(truncq, columns=COLUMNS)


def _load_spec(name, **values):
    with open(os.path.join(TEMPLATES, name)) as fh:
        spec = json.load(fh)
    return _json_replace(spec, **values)


def _get_specs(output_dir, datafiles):
    """The Vega specs of the Q scores and the expected errors tabs."""
    eestats, kept, lengths, truncq = [], [], [], []
    for direction in datafiles:
        files = datafiles[direction]
        ee = _read_eestats(os.path.join(output_dir, files['eestats']),
                           direction)
        eestats.append(ee)
        # Reads at least this long, so kept by truncating with no EE filter
        kept.append(pd.DataFrame({'direction': direction,
                                  'Length': ee['Pos'],
                                  'series': LENGTH_ONLY,
                                  'PctReads': ee['PctRecs']}))
        kept.append(_read_eestats2(
            os.path.join(output_dir, files['eestats2']), direction))
        length, trunc = _read_fastq_stats(
            os.path.join(output_dir, files['stats']), direction)
        lengths.append(length)
        truncq.append(trunc)
    eestats = pd.concat(eestats).to_dict('records')
    truncq = pd.concat(truncq)

    q_spec = _load_spec(
        'q_scores_spec.json',
        eestats=eestats,
        lengths=pd.concat(lengths).to_dict('records'),
        truncq=truncq.to_dict('records'),
        q_series=list(dict.fromkeys(truncq['series'])),
        directions=list(datafiles))
    ee_spec = _load_spec(
        'expected_errors_spec.json',
        eestats=eestats,
        kept=pd.concat(kept).to_dict('records'),
        maxees=[{'maxee': ee, 'series': _maxee_name(ee)}
                for ee in EE_CUTOFFS],
        series=[LENGTH_ONLY] + [_maxee_name(ee) for ee in EE_CUTOFFS],
        directions=list(datafiles))
    return q_spec, ee_spec


def _fastq_stats(output_dir: str, sequences, threads) -> None:
    # read manifest
    manifest = sequences.manifest
    # check if paired reads available
    try:
        paired = manifest['reverse'].iloc[0] is not None
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

    q_spec, ee_spec = _get_specs(output_dir, datafiles)
    index = os.path.join(TEMPLATES, 'index.html')
    ee_template = os.path.join(TEMPLATES, 'expected_errors.html')
    stats_template = os.path.join(TEMPLATES, 'fastq_stats.html')
    eestats_template = os.path.join(TEMPLATES, 'fastq_eestats.html')
    eestats2_template = os.path.join(TEMPLATES, 'fastq_eestats2.html')
    context = {
        'paired': paired,
        'datafiles': datafiles,
        'html': html,
        'q_spec': json.dumps(q_spec),
        'ee_spec': json.dumps(ee_spec),
        'tabs': [{'title': 'Q scores',
                  'url': 'index.html'},
                 {'title': 'Expected Error',
                  'url': 'expected_errors.html'},
                 {'title': 'fastq_stats',
                  'url': 'fastq_stats.html'},
                 {'title': 'fastq_eestats',
                  'url': 'fastq_eestats.html'},
                 {'title': 'fastq_eestats2',
                  'url': 'fastq_eestats2.html'}],
        }
    templates = [index, ee_template, stats_template, eestats_template,
                 eestats2_template]
    q2templates.util.copy_assets(os.path.join(TEMPLATES, 'vega'), output_dir)
    q2templates.render(templates, output_dir, context=context)


def fastq_stats(output_dir: str,
                sequences: CasavaOneEightSingleLanePerSampleDirFmt,
                threads: int = 1
                ) -> None:
    _fastq_stats(output_dir, sequences, threads)
