# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile

from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences import QIIME1DemuxDirFmt
import biom
from biom.parse import parse_uc
import skbio.io

from ._cluster_features import run_command


def dereplicate_sequences(sequences: QIIME1DemuxDirFmt,
                          derep_prefix: bool = False,
                          min_seq_length: int = 1,
                          min_unique_size: int = 1
                          ) -> (biom.Table, DNAFASTAFormat):
    dereplicated_sequences = DNAFASTAFormat()
    with tempfile.NamedTemporaryFile(mode='w+') as out_uc:
        seqs_fp = '%s/seqs.fna' % str(sequences)
        cmd = ['vsearch',
               '--derep_prefix' if derep_prefix else '--derep_fulllength',
               seqs_fp,
               '--output', str(dereplicated_sequences),
               '--relabel_sha1', '--relabel_keep',
               '--uc', out_uc.name,
               '--xsize',
               '--minseqlength', str(min_seq_length),
               '--minuniquesize', str(min_unique_size),
               '--fasta_width', '0']
        run_command(cmd)
        out_uc.seek(0)
        table = parse_uc(out_uc)
    id_map = {e.metadata['description']: e.metadata['id']
              for e in skbio.io.read(str(dereplicated_sequences),
                                     constructor=skbio.DNA,
                                     format='fasta')}
    table.update_ids(id_map=id_map, axis='observation')
    return table, dereplicated_sequences
