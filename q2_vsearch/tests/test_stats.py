# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import json
import tempfile
import glob

from qiime2.plugin.testing import TestPluginBase
from qiime2 import Artifact

from q2_vsearch._stats import (
    _read_eestats, _read_eestats2, _read_fastq_stats, _json_replace,
    EE_CUTOFFS)


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
            self._test_plots(output_dir, paired)

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

    def _spec(self, output_dir, page):
        with open(os.path.join(output_dir, page)) as fh:
            html = fh.read()
        spec = json.loads(html.split(
            '<script id="spec" type="application/json">')[1].split(
            '</script>')[0])
        self.assertNotIn('{{REPLACE_PARAM}}', json.dumps(spec))
        return spec, {d['name']: d.get('values') for d in spec['data']}

    def _test_plots(self, output_dir, paired):
        for asset in ['js/vega.min.js', 'js/vega-embed.min.js']:
            self.assertTrue(os.path.exists(os.path.join(output_dir, asset)))
        directions = ['forward', 'reverse'] if paired else ['forward']

        spec, data = self._spec(output_dir, 'index.html')
        self.assertEqual(
            sorted({r['direction'] for r in data['eestats']}), directions)
        self.assertEqual(
            sorted({r['direction'] for r in data['lengths']}), directions)
        self.assertEqual(sorted({r['series'] for r in data['truncq']}),
                         ['first Q \u2264 10', 'first Q \u2264 15',
                          'first Q \u2264 20', 'first Q \u2264 5'])

        spec, data = self._spec(output_dir, 'expected_errors.html')
        # 151 bp reads: one row per length for each max EE and length alone
        self.assertEqual(len(data['kept']),
                         151 * (len(EE_CUTOFFS) + 1) * len(directions))
        self.assertEqual([r['maxee'] for r in data['maxees']], EE_CUTOFFS)

    def test_read_fastq_stats(self):
        lengths, truncq = _read_fastq_stats(
            self.get_data_path('stats/fastq_stats_forward.txt'), 'forward')
        self.assertEqual(lengths.to_dict('records'),
                         [{'direction': 'forward', 'Length': 151,
                           'Reads': 2357, 'PctReads': 100.0}])
        self.assertEqual(truncq['series'].unique().tolist(),
                         ['first Q \u2264 5', 'first Q \u2264 10',
                          'first Q \u2264 15', 'first Q \u2264 20'])
        self.assertEqual(truncq['Length'].max(), 151)

    def test_read_eestats2(self):
        obs = _read_eestats2(
            self.get_data_path('stats/fastq_eestats2_forward.txt'), 'forward')
        self.assertEqual(obs.shape, (9, 4))
        self.assertEqual(obs['Length'].unique().tolist(), [50, 100, 150])
        self.assertEqual(obs['series'].unique().tolist(),
                         ['max EE 0.5', 'max EE 1', 'max EE 2'])
        self.assertEqual(obs.iloc[-1].tolist(),
                         ['forward', 150, 'max EE 2', 98.3])

    def test_read_eestats(self):
        obs = _read_eestats(
            self.get_data_path('stats/fastq_eestats_forward.txt'), 'reverse')
        self.assertEqual(obs['direction'].unique().tolist(), ['reverse'])
        self.assertEqual(obs['Pos'].iloc[0], 1)
        self.assertIn('Med_EE', obs.columns)

    def test_json_replace(self):
        obs = _json_replace({'a': [{'{{REPLACE_PARAM}}': 'x'}, 1],
                             'b': {'{{REPLACE_PARAM}}': 'y'}},
                            x=[1, 2], y='z')
        self.assertEqual(obs, {'a': [[1, 2], 1], 'b': 'z'})
