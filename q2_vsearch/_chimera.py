# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile

import biom
import skbio.io
from q2_types.feature_data import DNAFASTAFormat

from ._cluster_features import _fasta_with_sizes, run_command
from ._format import UchimeStatsFmt


_uchime_defaults = {'dn': 1.4,
                    'mindiffs': 3,
                    'mindiv': 0.8,
                    'minh': 0.28,
                    'xn': 8.0}


def uchime_ref(sequences: DNAFASTAFormat,
               table: biom.Table,
               reference_sequences: DNAFASTAFormat,
               dn: float = _uchime_defaults['dn'],
               mindiffs: int = _uchime_defaults['mindiffs'],
               mindiv: float = _uchime_defaults['mindiv'],
               minh: float = _uchime_defaults['minh'],
               xn: float = _uchime_defaults['xn'],
               threads: int = 1) \
               -> (DNAFASTAFormat, DNAFASTAFormat, UchimeStatsFmt):
    cmd, chimeras, nonchimeras, uchime_stats = \
        _uchime_ref(sequences, table, reference_sequences, dn, mindiffs,
                    mindiv, minh, xn, threads)
    return chimeras, nonchimeras, uchime_stats


def _uchime_ref(sequences, table, reference_sequences, dn, mindiffs,
                mindiv, minh, xn, threads):
    # this function only exists to simplify testing
    chimeras = DNAFASTAFormat()
    nonchimeras = DNAFASTAFormat()
    uchime_stats = UchimeStatsFmt()
    with tempfile.NamedTemporaryFile() as fasta_with_sizes:
        with tempfile.NamedTemporaryFile() as temp_chimeras:
            _fasta_with_sizes(str(sequences), fasta_with_sizes.name, table)
            cmd = ['vsearch',
                   '--uchime_ref', fasta_with_sizes.name,
                   '--uchimeout', str(uchime_stats),
                   '--nonchimeras', str(nonchimeras),
                   '--chimeras', temp_chimeras.name,
                   '--dn', str(dn),
                   '--mindiffs', str(mindiffs),
                   '--mindiv', str(mindiv),
                   '--minh', str(minh),
                   '--xn', str(xn),
                   '--db', str(reference_sequences),
                   '--qmask', 'none',  # ensures no lowercase DNA chars
                   '--xsize',
                   '--threads', str(threads),
                   '--minseqlength', '1',
                   '--fasta_width', '0']
            run_command(cmd)
            # this processing step should be removed, pending fix of:
            # https://github.com/qiime2/q2-vsearch/issues/39
            _fix_chimera_ids(temp_chimeras, chimeras)

    return cmd, chimeras, nonchimeras, uchime_stats


def uchime_denovo(sequences: DNAFASTAFormat,
                  table: biom.Table,
                  dn: float = _uchime_defaults['dn'],
                  mindiffs: int = _uchime_defaults['mindiffs'],
                  mindiv: float = _uchime_defaults['mindiv'],
                  minh: float = _uchime_defaults['minh'],
                  xn: float = _uchime_defaults['xn']) \
                  -> (DNAFASTAFormat, DNAFASTAFormat, UchimeStatsFmt):
    cmd, chimeras, nonchimeras, uchime_stats = \
        _uchime_denovo(sequences, table, dn, mindiffs, mindiv, minh, xn)
    return chimeras, nonchimeras, uchime_stats


def _uchime_denovo(sequences, table, dn, mindiffs, mindiv, minh, xn):
    # this function only exists to simplify testing
    chimeras = DNAFASTAFormat()
    nonchimeras = DNAFASTAFormat()
    uchime_stats = UchimeStatsFmt()
    with tempfile.NamedTemporaryFile() as fasta_with_sizes:
        with tempfile.NamedTemporaryFile() as temp_chimeras:
            _fasta_with_sizes(str(sequences), fasta_with_sizes.name, table)
            cmd = ['vsearch',
                   '--uchime_denovo', fasta_with_sizes.name,
                   '--uchimeout', str(uchime_stats),
                   '--nonchimeras', str(nonchimeras),
                   '--chimeras', temp_chimeras.name,
                   '--dn', str(dn),
                   '--mindiffs', str(mindiffs),
                   '--mindiv', str(mindiv),
                   '--minh', str(minh),
                   '--xn', str(xn),
                   '--qmask', 'none',  # ensures no lowercase DNA chars
                   '--xsize',
                   '--minseqlength', '1',
                   '--fasta_width', '0']
            run_command(cmd)
            # this processing step should be removed, pending fix of:
            # https://github.com/qiime2/q2-vsearch/issues/39
            _fix_chimera_ids(temp_chimeras, chimeras)

    return cmd, chimeras, nonchimeras, uchime_stats


def _fix_chimera_ids(temp_chimeras, output_chimeras):
    # this processing function should be removed, pending fix of:
    # https://github.com/qiime2/q2-vsearch/issues/39
    temp_chimeras.seek(0)
    with open(str(output_chimeras), 'w') as out_fh:
        for seq in skbio.io.read(temp_chimeras, format='fasta',
                                 constructor=skbio.DNA):
            seq.metadata['id'] = seq.metadata['id'].rsplit(';', 1)[0]
            seq.write(out_fh)
