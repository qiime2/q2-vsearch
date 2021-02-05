# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
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


_jp_defaults = {
    'truncqual': None,
    'minlen': 1,
    'maxns': None,
    'allowmergestagger': False,
    'minovlen': 10,
    'maxdiffs': 10,
    'minmergelen': None,
    'maxmergelen': None,
    'maxee': None,
    'qmin': 0,
    'qminout': 0,
    'qmax': 41,
    'qmaxout': 41,
    'threads': 1
}


def join_pairs(demultiplexed_seqs: SingleLanePerSamplePairedEndFastqDirFmt,
               truncqual: int = _jp_defaults['truncqual'],
               minlen: int = _jp_defaults['minlen'],
               maxns: int = _jp_defaults['maxns'],
               allowmergestagger: bool = _jp_defaults['allowmergestagger'],
               minovlen: int = _jp_defaults['minovlen'],
               maxdiffs: int = _jp_defaults['maxdiffs'],
               minmergelen: int = _jp_defaults['minmergelen'],
               maxmergelen: int = _jp_defaults['maxmergelen'],
               maxee: float = _jp_defaults['maxee'],
               qmin: int = _jp_defaults['qmin'],
               qminout: int = _jp_defaults['qminout'],
               qmax: int = _jp_defaults['qmax'],
               qmaxout: int = _jp_defaults['qmaxout'],
               threads: int = _jp_defaults['threads'],
               ) -> SingleLanePerSampleSingleEndFastqDirFmt:
    _, result = _join_pairs_w_command_output(
        demultiplexed_seqs, truncqual, minlen, maxns, allowmergestagger,
        minovlen, maxdiffs, minmergelen, maxmergelen, maxee, qmin, qminout,
        qmax, qmaxout, threads)
    return result


def _join_pairs_w_command_output(
        demultiplexed_seqs: SingleLanePerSamplePairedEndFastqDirFmt,
        truncqual: int = _jp_defaults['truncqual'],
        minlen: int = _jp_defaults['minlen'],
        maxns: int = _jp_defaults['maxns'],
        allowmergestagger: bool = _jp_defaults['allowmergestagger'],
        minovlen: int = _jp_defaults['minovlen'],
        maxdiffs: int = _jp_defaults['maxdiffs'],
        minmergelen: int = _jp_defaults['minmergelen'],
        maxmergelen: int = _jp_defaults['maxmergelen'],
        maxee: float = _jp_defaults['maxee'],
        qmin: int = _jp_defaults['qmin'],
        qminout: int = _jp_defaults['qminout'],
        qmax: int = _jp_defaults['qmax'],
        qmaxout: int = _jp_defaults['qmaxout'],
        threads: int = _jp_defaults['threads'],
        ) -> (List[str], SingleLanePerSampleSingleEndFastqDirFmt):
    # this function exists only to simplify unit testing

    result = SingleLanePerSampleSingleEndFastqDirFmt()

    manifest = pd.read_csv(
        os.path.join(str(demultiplexed_seqs),
                     demultiplexed_seqs.manifest.pathspec),
        header=0, comment='#')
    manifest.filename = manifest.filename.apply(
        lambda x: os.path.join(str(demultiplexed_seqs), x))

    phred_offset = yaml.load(open(
        os.path.join(str(demultiplexed_seqs),
                     demultiplexed_seqs.metadata.pathspec)))['phred-offset']

    id_to_fps = manifest.pivot(index='sample-id', columns='direction',
                               values='filename')

    output_manifest = FastqManifestFormat()
    output_manifest_fh = output_manifest.open()
    output_manifest_fh.write('sample-id,filename,direction\n')
    output_manifest_fh.write('# direction is not meaningful in this file '
                             'as these\n')
    output_manifest_fh.write('# data may be derived from forward, reverse, '
                             'or \n')
    output_manifest_fh.write('# joined reads\n')

    for i, (sample_id, (fwd_fp, rev_fp)) in enumerate(id_to_fps.iterrows()):
        # The barcode id, lane number and read number are not relevant
        # here. We might ultimately want to use a dir format other than
        # SingleLanePerSampleSingleEndFastqDirFmt which doesn't care
        # about this information. Similarly, the direction of the read
        # isn't relevant here anymore.
        path = result.sequences.path_maker(sample_id=sample_id,
                                           barcode_id=i,
                                           lane_number=1,
                                           read_number=1)
        uncompressed_path = str(path).strip('.gz')

        cmd = ['vsearch',
               '--fastq_mergepairs', fwd_fp,
               '--reverse', rev_fp,
               '--fastqout', uncompressed_path,
               '--fastq_ascii', str(phred_offset),
               '--fastq_minlen', str(minlen),
               '--fastq_minovlen', str(minovlen),
               '--fastq_maxdiffs', str(maxdiffs),
               '--fastq_qmin', str(qmin),
               '--fastq_qminout', str(qminout),
               '--fastq_qmax', str(qmax),
               '--fastq_qmaxout', str(qmaxout),
               '--minseqlength', '1',
               '--fasta_width', '0']
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
        run_command(['gzip', uncompressed_path])
        output_manifest_fh.write(
            '%s,%s,%s\n' % (sample_id, path.name, 'forward'))

    output_manifest_fh.close()
    result.manifest.write_data(output_manifest, FastqManifestFormat)

    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': phred_offset}))
    result.metadata.write_data(metadata, YamlFormat)

    return cmd, result
