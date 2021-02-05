# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

import skbio
import numpy as np
import pandas as pd
from qiime2.plugin.testing import TestPluginBase
from qiime2.util import redirected_stdio
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    FastqGzFormat)

from q2_vsearch._join_pairs import join_pairs, _join_pairs_w_command_output


class MergePairsTests(TestPluginBase):

    package = 'q2_vsearch.tests'

    def setUp(self):
        super().setUp()
        self.input_seqs = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path('demux-1'), 'r')

    def _parse_manifest(self, demultiplexed_seqs):
        return pd.read_csv(
                os.path.join(str(demultiplexed_seqs),
                             demultiplexed_seqs.manifest.pathspec),
                header=0, comment='#')

    def _test_manifest(self, demultiplexed_seqs):
        manifest = self._parse_manifest(demultiplexed_seqs)
        self.assertEqual(len(manifest), 3)
        self.assertEqual(list(manifest['sample-id']),
                         ['BAQ2687.1', 'BAQ3473.2', 'BAQ4697.2'])
        self.assertEqual(list(manifest['filename']),
                         ['BAQ2687.1_0_L001_R1_001.fastq.gz',
                          'BAQ3473.2_1_L001_R1_001.fastq.gz',
                          'BAQ4697.2_2_L001_R1_001.fastq.gz'])
        self.assertEqual(list(manifest['direction']),
                         ['forward', 'forward', 'forward'])

    def _test_seq_lengths(self, seq_lengths):
        self.assertTrue(seq_lengths.mean() > 200)
        for e in seq_lengths:
            # input reads are 151 bases, so all output must be longer
            self.assertTrue(e > 151)

    def test_join_pairs(self):

        with redirected_stdio(stderr=os.devnull):
            obs = join_pairs(self.input_seqs)

        # manifest is as expected
        self._test_manifest(obs)

        # expected number of fastq files are created
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # The following values were determined by running vsearch directly
        # with default parameters. It is possible that different versions of
        # vsearch will result in differences in these numbers, and that
        # the corresponding tests may therefore be too specific. We'll have
        # to adjust the tests if that's the case.
        default_exp_sequence_counts = {
            'BAQ2687.1_0_L001_R1_001.fastq.gz': 806,
            'BAQ3473.2_1_L001_R1_001.fastq.gz': 753,
            'BAQ4697.2_2_L001_R1_001.fastq.gz': 711,
        }
        for fastq_name, fastq_path in output_fastqs:
            seqs = skbio.io.read(str(fastq_path), format='fastq',
                                 compression='gzip', constructor=skbio.DNA)
            seqs = list(seqs)
            seq_lengths = np.asarray([len(s) for s in seqs])
            self._test_seq_lengths(seq_lengths)

            # expected number of sequences are joined
            self.assertEqual(
                len(seq_lengths),
                default_exp_sequence_counts[str(fastq_name)])

    def test_join_pairs_some_samples_w_no_joined_seqs(self):
        # minmergelen is set very high here, resulting in only one sequence
        # being joined across the three samples.
        with redirected_stdio(stderr=os.devnull):
            obs = join_pairs(self.input_seqs, minmergelen=279)

        # manifest is as expected
        self._test_manifest(obs)

        # expected number of fastq files are created
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # The following values were determined by running vsearch directly.
        exp_sequence_counts = {
            'BAQ2687.1_0_L001_R1_001.fastq.gz': 0,
            'BAQ3473.2_1_L001_R1_001.fastq.gz': 2,
            'BAQ4697.2_2_L001_R1_001.fastq.gz': 0,
        }

        for fastq_name, fastq_path in output_fastqs:
            with redirected_stdio(stderr=os.devnull):
                seqs = skbio.io.read(str(fastq_path), format='fastq',
                                     compression='gzip', constructor=skbio.DNA)
            seqs = list(seqs)
            seq_lengths = np.asarray([len(s) for s in seqs])

            # expected number of sequences are joined
            self.assertEqual(
                len(seq_lengths),
                exp_sequence_counts[str(fastq_name)])

    def test_join_pairs_all_samples_w_no_joined_seqs(self):
        # minmergelen is set very high here, resulting in no sequences
        # being joined across the three samples.
        with redirected_stdio(stderr=os.devnull):
            obs = join_pairs(self.input_seqs, minmergelen=500)

        # manifest is as expected
        self._test_manifest(obs)

        # expected number of fastq files are created
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        for fastq_name, fastq_path in output_fastqs:
            with redirected_stdio(stderr=os.devnull):
                seqs = skbio.io.read(str(fastq_path), format='fastq',
                                     compression='gzip', constructor=skbio.DNA)
            seqs = list(seqs)
            seq_lengths = np.asarray([len(s) for s in seqs])

            self.assertEqual(len(seq_lengths), 0)

    def test_join_pairs_alt_truncqual(self):
        with redirected_stdio(stderr=os.devnull):
            cmd, obs = _join_pairs_w_command_output(
                self.input_seqs, truncqual=5)

        # sanity check the output
        self._test_manifest(obs)
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # confirm altered parameter was passed to vsearch
        self.assertTrue('--fastq_truncqual 5' in ' '.join(cmd))

    def test_join_pairs_alt_minlen(self):
        with redirected_stdio(stderr=os.devnull):
            cmd, obs = _join_pairs_w_command_output(
                self.input_seqs, minlen=25)

        # sanity check the output
        self._test_manifest(obs)
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # confirm altered parameter was passed to vsearch
        self.assertTrue('--fastq_minlen 25' in ' '.join(cmd))

    def test_join_pairs_alt_maxns(self):
        with redirected_stdio(stderr=os.devnull):
            cmd, obs = _join_pairs_w_command_output(
                self.input_seqs, maxns=2)

        # sanity check the output
        self._test_manifest(obs)
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # confirm altered parameter was passed to vsearch
        self.assertTrue('--fastq_maxns 2' in ' '.join(cmd))

    def test_join_pairs_alt_allowmergestagger(self):
        with redirected_stdio(stderr=os.devnull):
            cmd, obs = _join_pairs_w_command_output(
                self.input_seqs, allowmergestagger=True)

        # sanity check the output
        self._test_manifest(obs)
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # confirm altered parameter was passed to vsearch
        self.assertTrue('--fastq_allowmergestagger' in cmd)

    def test_join_pairs_alt_minovlen(self):
        with redirected_stdio(stderr=os.devnull):
            cmd, obs = _join_pairs_w_command_output(
                self.input_seqs, minovlen=42)

        # sanity check the output
        self._test_manifest(obs)
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # confirm altered parameter was passed to vsearch
        self.assertTrue('--fastq_minovlen 42' in ' '.join(cmd))

    def test_join_pairs_alt_maxdiffs(self):
        with redirected_stdio(stderr=os.devnull):
            cmd, obs = _join_pairs_w_command_output(
                self.input_seqs, maxdiffs=2)

        # sanity check the output
        self._test_manifest(obs)
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # confirm altered parameter was passed to vsearch
        self.assertTrue('--fastq_maxdiffs 2' in ' '.join(cmd))

    def test_join_pairs_alt_minmergelen(self):
        with redirected_stdio(stderr=os.devnull):
            cmd, obs = _join_pairs_w_command_output(
                self.input_seqs, minmergelen=250)

        # sanity check the output
        self._test_manifest(obs)
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # confirm altered parameter was passed to vsearch
        self.assertTrue('--fastq_minmergelen 250' in ' '.join(cmd))

    def test_join_pairs_alt_maxmergelen(self):
        with redirected_stdio(stderr=os.devnull):
            cmd, obs = _join_pairs_w_command_output(
                self.input_seqs, maxmergelen=250)

        # sanity check the output
        self._test_manifest(obs)
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # confirm altered parameter was passed to vsearch
        self.assertTrue('--fastq_maxmergelen 250' in ' '.join(cmd))

    def test_join_pairs_alt_maxee(self):
        with redirected_stdio(stderr=os.devnull):
            cmd, obs = _join_pairs_w_command_output(
                self.input_seqs, maxee=25.0)

        # sanity check the output
        self._test_manifest(obs)
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # confirm altered parameter was passed to vsearch
        self.assertTrue('--fastq_maxee 25.0' in ' '.join(cmd))

    def test_join_pairs_alt_qmin(self):
        with redirected_stdio(stderr=os.devnull):
            cmd, obs = _join_pairs_w_command_output(
                self.input_seqs, qmin=-1)

        # sanity check the output
        self._test_manifest(obs)
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # confirm altered parameter was passed to vsearch
        self.assertTrue('--fastq_qmin -1' in ' '.join(cmd))

    def test_join_pairs_alt_qminout(self):
        with redirected_stdio(stderr=os.devnull):
            cmd, obs = _join_pairs_w_command_output(
                self.input_seqs, qminout=-1)

        # sanity check the output
        self._test_manifest(obs)
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # confirm altered parameter was passed to vsearch
        self.assertTrue('--fastq_qminout -1' in ' '.join(cmd))

    def test_join_pairs_alt_qmax(self):
        with redirected_stdio(stderr=os.devnull):
            cmd, obs = _join_pairs_w_command_output(
                self.input_seqs, qmax=40)

        # sanity check the output
        self._test_manifest(obs)
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # confirm altered parameter was passed to vsearch
        self.assertTrue('--fastq_qmax 40' in ' '.join(cmd))

    def test_join_pairs_alt_qmaxout(self):
        with redirected_stdio(stderr=os.devnull):
            cmd, obs = _join_pairs_w_command_output(
                self.input_seqs, qmaxout=40)

        # sanity check the output
        self._test_manifest(obs)
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # confirm altered parameter was passed to vsearch
        self.assertTrue('--fastq_qmaxout 40' in ' '.join(cmd))

    def test_join_pairs_alt_threads(self):
        with redirected_stdio(stderr=os.devnull):
            cmd, obs = _join_pairs_w_command_output(
                self.input_seqs, threads=2)

        # sanity check the output
        self._test_manifest(obs)
        output_fastqs = list(obs.sequences.iter_views(FastqGzFormat))
        self.assertEqual(len(output_fastqs), 3)

        # confirm altered parameter was passed to vsearch
        self.assertTrue('--threads 2' in ' '.join(cmd))
