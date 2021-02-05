# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

import biom
import numpy as np

import qiime2
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin import ValidationError
from qiime2.util import redirected_stdio
from q2_types.feature_data import DNAFASTAFormat
from q2_vsearch._chimera import (uchime_denovo, _uchime_denovo,
                                 uchime_ref, _uchime_ref)
from q2_vsearch._format import UchimeStatsFmt
from .test_cluster_features import _read_seqs


class UchimeDenovoTests(TestPluginBase):

    package = 'q2_vsearch.tests'

    def setUp(self):
        super().setUp()
        input_sequences_fp = self.get_data_path('dna-sequences-3.fasta')
        self.input_sequences = DNAFASTAFormat(input_sequences_fp, mode='r')
        self.input_sequences_list = _read_seqs(self.input_sequences)

        self.input_table = biom.Table(np.array([[100, 101, 103],
                                                [99, 98, 99],
                                                [4, 5, 6],
                                                [2, 2, 2]]),
                                      ['feature1', 'feature2', 'feature3',
                                       'feature4'],
                                      ['sample1', 'sample2', 'sample3'])

    def test_uchime_denovo(self):
        with redirected_stdio(stderr=os.devnull):
            chime, nonchime, stats = uchime_denovo(
                sequences=self.input_sequences, table=self.input_table)

        obs_chime = _read_seqs(chime)
        exp_chime = [self.input_sequences_list[3]]
        self.assertEqual(obs_chime, exp_chime)

        # sequences are reverse-sorted by abundance in output
        obs_nonchime = _read_seqs(nonchime)
        exp_nonchime = [self.input_sequences_list[0],
                        self.input_sequences_list[1],
                        self.input_sequences_list[2]]
        self.assertEqual(obs_nonchime, exp_nonchime)

        with stats.open() as stats_fh:
            stats_text = stats_fh.read()
        self.assertTrue('feature1' in stats_text)
        self.assertTrue('feature2' in stats_text)
        self.assertTrue('feature3' in stats_text)
        self.assertTrue('feature4' in stats_text)
        stats_lines = [e for e in stats_text.split('\n')
                       if len(e) > 0]
        self.assertEqual(len(stats_lines), 4)

    def test_uchime_denovo_no_chimeras(self):
        input_table = biom.Table(np.array([[3, 4, 2],
                                           [1, 0, 0],
                                           [4, 5, 6],
                                           [2, 2, 2]]),
                                 ['feature1', 'feature2', 'feature3',
                                  'feature4'],
                                 ['sample1', 'sample2', 'sample3'])
        with redirected_stdio(stderr=os.devnull):
            chime, nonchime, stats = uchime_denovo(
                sequences=self.input_sequences, table=input_table)

        obs_chime = _read_seqs(chime)
        exp_chime = []
        self.assertEqual(obs_chime, exp_chime)

        # sequences are reverse-sorted by abundance in output
        obs_nonchime = _read_seqs(nonchime)
        exp_nonchime = [self.input_sequences_list[2],
                        self.input_sequences_list[0],
                        self.input_sequences_list[3],
                        self.input_sequences_list[1]]
        self.assertEqual(obs_nonchime, exp_nonchime)

        with stats.open() as stats_fh:
            stats_text = stats_fh.read()
        self.assertTrue('feature1' in stats_text)
        self.assertTrue('feature2' in stats_text)
        self.assertTrue('feature3' in stats_text)
        self.assertTrue('feature4' in stats_text)
        stats_lines = [e for e in stats_text.split('\n')
                       if len(e) > 0]
        self.assertEqual(len(stats_lines), 4)

    def test_uchime_denovo_no_chimeras_alt_params(self):

        with redirected_stdio(stderr=os.devnull):
            cmd, chime, nonchime, stats = _uchime_denovo(
                sequences=self.input_sequences, table=self.input_table,
                dn=42.42, mindiffs=4, mindiv=0.5, minh=0.42, xn=9.0)
        cmd = ' '.join(cmd)
        self.assertTrue('--dn 42.42' in cmd)
        self.assertTrue('--mindiffs 4' in cmd)
        self.assertTrue('--mindiv 0.5' in cmd)
        self.assertTrue('--minh 0.42' in cmd)
        self.assertTrue('--xn 9.0' in cmd)


class UchimeRefTests(TestPluginBase):

    package = 'q2_vsearch.tests'

    def setUp(self):
        super().setUp()
        input_sequences_fp = self.get_data_path('dna-sequences-3.fasta')
        self.input_sequences = DNAFASTAFormat(input_sequences_fp, mode='r')
        self.input_sequences_list = _read_seqs(self.input_sequences)

        self.input_table = biom.Table(np.array([[100, 101, 103],
                                                [99, 98, 99],
                                                [4, 5, 6],
                                                [2, 2, 2]]),
                                      ['feature1', 'feature2', 'feature3',
                                       'feature4'],
                                      ['sample1', 'sample2', 'sample3'])

    def test_uchime_ref(self):
        ref_sequences_fp = self.get_data_path('ref-sequences-3.fasta')
        ref_sequences = DNAFASTAFormat(ref_sequences_fp, mode='r')

        with redirected_stdio(stderr=os.devnull):
            chime, nonchime, stats = uchime_ref(
                sequences=self.input_sequences, table=self.input_table,
                reference_sequences=ref_sequences)

        obs_chime = _read_seqs(chime)
        exp_chime = [self.input_sequences_list[3]]
        self.assertEqual(obs_chime, exp_chime)

        # sequences are reverse-sorted by abundance in output
        obs_nonchime = _read_seqs(nonchime)
        exp_nonchime = [self.input_sequences_list[0],
                        self.input_sequences_list[1],
                        self.input_sequences_list[2]]
        self.assertEqual(obs_nonchime, exp_nonchime)

        with stats.open() as stats_fh:
            stats_text = stats_fh.read()
        self.assertTrue('feature1' in stats_text)
        self.assertTrue('feature2' in stats_text)
        self.assertTrue('feature3' in stats_text)
        self.assertTrue('feature4' in stats_text)
        stats_lines = [e for e in stats_text.split('\n')
                       if len(e) > 0]
        self.assertEqual(len(stats_lines), 4)

    def test_uchime_ref_no_chimeras(self):
        ref_sequences_fp = self.get_data_path('ref-sequences-4.fasta')
        ref_sequences = DNAFASTAFormat(ref_sequences_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            chime, nonchime, stats = uchime_ref(
                sequences=self.input_sequences, table=self.input_table,
                reference_sequences=ref_sequences)

        obs_chime = _read_seqs(chime)
        exp_chime = []
        self.assertEqual(obs_chime, exp_chime)

        # sequences are reverse-sorted by abundance in output
        obs_nonchime = _read_seqs(nonchime)
        exp_nonchime = [self.input_sequences_list[0],
                        self.input_sequences_list[1],
                        self.input_sequences_list[2],
                        self.input_sequences_list[3]]
        self.assertEqual(obs_nonchime, exp_nonchime)

        with stats.open() as stats_fh:
            stats_text = stats_fh.read()
        self.assertTrue('feature1' in stats_text)
        self.assertTrue('feature2' in stats_text)
        self.assertTrue('feature3' in stats_text)
        self.assertTrue('feature4' in stats_text)
        stats_lines = [e for e in stats_text.split('\n')
                       if len(e) > 0]
        self.assertEqual(len(stats_lines), 4)

    def test_uchime_ref_no_chimeras_alt_params(self):
        ref_sequences_fp = self.get_data_path('ref-sequences-4.fasta')
        ref_sequences = DNAFASTAFormat(ref_sequences_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            cmd, chime, nonchime, stats = _uchime_ref(
                sequences=self.input_sequences, table=self.input_table,
                reference_sequences=ref_sequences, dn=42.42, mindiffs=4,
                mindiv=0.5, minh=0.42, xn=9.0, threads=2)
        cmd = ' '.join(cmd)
        self.assertTrue('--dn 42.42' in cmd)
        self.assertTrue('--mindiffs 4' in cmd)
        self.assertTrue('--mindiv 0.5' in cmd)
        self.assertTrue('--minh 0.42' in cmd)
        self.assertTrue('--xn 9.0' in cmd)
        self.assertTrue('--threads 2' in cmd)


class UchimeStatsFmtTests(TestPluginBase):
    package = 'q2_vsearch.tests'

    def test_validate_positive(self):
        filepath = self.get_data_path('uchime-stats-1.txt')
        format = UchimeStatsFmt(filepath, mode='r')

        format.validate(level='min')
        format.validate(level='max')

    def test_validate_negative(self):
        filepath = self.get_data_path('uchime-stats-invalid-1.txt')
        format = UchimeStatsFmt(filepath, mode='r')

        with self.assertRaisesRegex(ValidationError, 'exactly 18'):
            format.validate(level='min')

        with self.assertRaisesRegex(ValidationError, 'exactly 18'):
            format.validate(level='max')

        filepath = self.get_data_path('uchime-stats-invalid-2.txt')
        format = UchimeStatsFmt(filepath, mode='r')

        # file appears valid with min
        format.validate(level='min')

        with self.assertRaisesRegex(ValidationError, 'exactly 18'):
            format.validate(level='max')

    def test_transform_to_metadata(self):
        filepath = self.get_data_path('uchime-stats-1.txt')
        format = UchimeStatsFmt(filepath, mode='r')
        transformer = self.get_transformer(UchimeStatsFmt, qiime2.Metadata)
        obs = transformer(format)
        self.assertEqual(obs.id_count, 6)
        self.assertEqual(obs.column_count, 17)
        self.assertEqual(obs.id_header, 'feature-id')
