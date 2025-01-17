# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import yaml
from typing import List

import pandas as pd
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    FastqManifestFormat, YamlFormat)

from ._cluster_features import run_command


_mp_defaults = {
    'truncqual': None,
    'minlen': 1,
    'maxns': None,
    'allowmergestagger': False,
    'minovlen': 10,
    'maxdiffs': 10,
    'minmergelen': None,
    'maxmergelen': None,
    'maxee': None,
    'threads': 1
}


def merge_pairs(
    demultiplexed_seqs: SingleLanePerSamplePairedEndFastqDirFmt,
    truncqual: int = _mp_defaults['truncqual'],
    minlen: int = _mp_defaults['minlen'],
    maxns: int = _mp_defaults['maxns'],
    allowmergestagger: bool = _mp_defaults['allowmergestagger'],
    minovlen: int = _mp_defaults['minovlen'],
    maxdiffs: int = _mp_defaults['maxdiffs'],
    minmergelen: int = _mp_defaults['minmergelen'],
    maxmergelen: int = _mp_defaults['maxmergelen'],
    maxee: float = _mp_defaults['maxee'],
    threads: int = _mp_defaults['threads'],
) -> (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt
):
    _, merged, unmerged = _merge_pairs_w_command_output(
        demultiplexed_seqs, truncqual, minlen, maxns, allowmergestagger,
        minovlen, maxdiffs, minmergelen, maxmergelen, maxee, threads
    )

    return merged, unmerged


def _merge_pairs_w_command_output(
    demultiplexed_seqs: SingleLanePerSamplePairedEndFastqDirFmt,
    truncqual: int = _mp_defaults['truncqual'],
    minlen: int = _mp_defaults['minlen'],
    maxns: int = _mp_defaults['maxns'],
    allowmergestagger: bool = _mp_defaults['allowmergestagger'],
    minovlen: int = _mp_defaults['minovlen'],
    maxdiffs: int = _mp_defaults['maxdiffs'],
    minmergelen: int = _mp_defaults['minmergelen'],
    maxmergelen: int = _mp_defaults['maxmergelen'],
    maxee: float = _mp_defaults['maxee'],
    threads: int = _mp_defaults['threads'],
) -> (
    List[str],
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt
):
    # this function exists only to simplify unit testing

    # create formats
    merged = SingleLanePerSampleSingleEndFastqDirFmt()
    unmerged = SingleLanePerSamplePairedEndFastqDirFmt()

    # create manifests
    merged_manifest = FastqManifestFormat()
    merged_manifest_fh = merged_manifest.open()
    unmerged_manifest = FastqManifestFormat()
    unmerged_manifest_fh = unmerged_manifest.open()

    # write manifest headers
    _write_manifest_header(merged_manifest_fh, add_warning=True)
    _write_manifest_header(unmerged_manifest_fh)

    # generate input reads iterable
    manifest = pd.read_csv(
        os.path.join(
            str(demultiplexed_seqs), demultiplexed_seqs.manifest.pathspec
        ),
        header=0,
        comment='#'
    )

    manifest.filename = manifest.filename.apply(
        lambda x: os.path.join(str(demultiplexed_seqs), x)
    )

    phred_offset = yaml.load(
        open(os.path.join(
            str(demultiplexed_seqs), demultiplexed_seqs.metadata.pathspec
        )),
        Loader=yaml.SafeLoader
    )['phred-offset']

    id_to_fps = manifest.pivot(
        index='sample-id', columns='direction', values='filename'
    )

    for i, (sample_id, (fwd_fp, rev_fp)) in enumerate(id_to_fps.iterrows()):
        # The barcode id and lane number are not relevant for either format.
        # We might ultimately want to use a dir format other than these which
        # doesn't care about this information.
        # The read number (direction) is only relevant for the unmerged reads.

        gz_merged_path, fq_merged_path = _get_output_paths(
            merged, sample_id, i, 1
        )
        gz_unmerged_fwd_path, fq_unmerged_fwd_path = _get_output_paths(
            unmerged, sample_id, i, 1
        )
        gz_unmerged_rev_path, fq_unmerged_rev_path = _get_output_paths(
            unmerged, sample_id, i, 2
        )

        cmd = [
            'vsearch',
            '--fastq_mergepairs', fwd_fp,
            '--reverse', rev_fp,
            '--fastqout', fq_merged_path,
            '--fastqout_notmerged_fwd', fq_unmerged_fwd_path,
            '--fastqout_notmerged_rev', fq_unmerged_rev_path,
            '--fastq_ascii', str(phred_offset),
            '--fastq_minlen', str(minlen),
            '--fastq_minovlen', str(minovlen),
            '--fastq_maxdiffs', str(maxdiffs),
            '--fastq_qmin', '0',
            '--fastq_qminout', '0',
            '--fastq_qmax', '41',
            '--fastq_qmaxout', '41',
            '--fasta_width', '0'
        ]
        if truncqual is not None:
            cmd += ['--fastq_truncqual', str(truncqual)]
        if maxns is not None:
            cmd += ['--fastq_maxns', str(maxns)]
        if minmergelen is not None:
            cmd += ['--fastq_minmergelen', str(minmergelen)]
        if maxmergelen is not None:
            cmd += ['--fastq_maxmergelen', str(maxmergelen)]
        if maxee is not None:
            cmd += ['--fastq_maxee', str(maxee)]
        cmd += ['--threads', str(threads)]
        if allowmergestagger:
            cmd.append('--fastq_allowmergestagger')

        run_command(cmd)

        run_command([
            'gzip', fq_merged_path, fq_unmerged_fwd_path, fq_unmerged_rev_path
        ])

        merged_manifest_fh.write(
            '%s,%s,%s\n' % (sample_id, gz_merged_path.name, 'forward')
        )
        unmerged_manifest_fh.write(
            '%s,%s,%s\n' % (sample_id, gz_unmerged_fwd_path.name, 'forward')
        )
        unmerged_manifest_fh.write(
            '%s,%s,%s\n' % (sample_id, gz_unmerged_rev_path.name, 'reverse')
        )

    merged_manifest_fh.close()
    unmerged_manifest_fh.close()
    merged.manifest.write_data(merged_manifest, FastqManifestFormat)
    unmerged.manifest.write_data(unmerged_manifest, FastqManifestFormat)

    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': phred_offset}))
    merged.metadata.write_data(metadata, YamlFormat)
    unmerged.metadata.write_data(metadata, YamlFormat)

    return cmd, merged, unmerged


def _get_output_paths(format_, sample_id, barcode_id, direction):
    path = format_.sequences.path_maker(
        sample_id=sample_id,
        barcode_id=barcode_id,
        lane_number=1,
        read_number=direction
    )
    return path, str(path).strip('.gz')


def _write_manifest_header(manifest_fh, add_warning=False):
    manifest_fh.write('sample-id,filename,direction\n')
    if add_warning:
        manifest_fh.write('')
        manifest_fh.write('# direction is not meaningful for joined reads\n')
