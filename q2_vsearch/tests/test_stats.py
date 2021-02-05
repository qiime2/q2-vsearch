# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import glob

from qiime2.plugin.testing import TestPluginBase
from qiime2 import Artifact


class StatsTests(TestPluginBase):
    package = 'q2_vsearch.tests'

    def setUp(self):
        super().setUp()
        self.input_seqs_paired = Artifact.import_data(
            'SampleData[PairedEndSequencesWithQuality]',
            self.get_data_path('demux-1'))
        self.input_seqs_single = Artifact.import_data(
            'SampleData[SequencesWithQuality]',
            self.get_data_path('demux-1_se'))
        self.viz = self.plugin.visualizers['fastq_stats']

    def _test_fastq_stats(self, paired=False, threads=1):
        default_filelist = ['fastq_stats_forward.txt',
                            'fastq_eestats2_forward.txt',
                            'fastq_eestats_forward.txt']
        if (paired):
            default_filelist.extend(['fastq_stats_reverse.txt',
                                     'fastq_eestats2_reverse.txt',
                                     'fastq_eestats_reverse.txt'])
        default_filelist.sort()

        with tempfile.TemporaryDirectory() as output_dir:
            if (paired):
                self.input_seqs = self.input_seqs_paired
            else:
                self.input_seqs = self.input_seqs_single

            self.result = self.viz(self.input_seqs, threads)
            self.result.visualization.export_data(output_dir)

            pattern = output_dir + '/*.txt'
            filelist = [os.path.basename(x) for x in glob.glob(pattern)]
            filelist.sort()

            self.assertListEqual(default_filelist, filelist)

            for filename in filelist:
                with open(os.path.join(output_dir, filename),
                          'r') as inputfile:
                    default = inputfile.readlines()

                with open(os.path.join(output_dir, filename),
                          'r') as inputfile:
                    data = inputfile.readlines()

                if (filename.startswith('fastq_stats_')):
                    default = default[3:-4]
                    data = data[3:-4]

                self.assertListEqual(default, data)

    def test_fastq_stats_single(self):
        self._test_fastq_stats(paired=False, threads=1)

    def test_fastq_stats_paired(self):
        self._test_fastq_stats(paired=True, threads=1)
