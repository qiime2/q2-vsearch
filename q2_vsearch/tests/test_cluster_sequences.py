# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import skbio
import biom
import numpy as np

from qiime2.plugin.testing import TestPluginBase
from qiime2.util import redirected_stdio
from q2_types.per_sample_sequences import QIIME1DemuxDirFmt

from q2_vsearch._cluster_sequences import (dereplicate_sequences,
                                           _parse_uc)


class DereplicateSequences(TestPluginBase):

    package = 'q2_vsearch.tests'

    def test_dereplicate_sequences(self):
        input_sequences_fp = self.get_data_path('seqs-1')
        input_sequences = QIIME1DemuxDirFmt(input_sequences_fp, 'r')

        exp_table = biom.Table(np.array([[2, 1],
                                         [0, 1],
                                         [0, 2]]),
                               ['4574b947a0159c0da35a1f30f989681a1d9f64ef',
                                '1768cf7fca79f84d651b34d878de2492c6a7b971',
                                '16a1263bde4f2f99422630d1bb87935c4236d1ba'],
                               ['sample1', 's2'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, obs_sequences = dereplicate_sequences(
                sequences=input_sequences)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        # sequences are reverse-sorted by abundance in output
        obs_seqs = list(skbio.io.read(str(obs_sequences),
                        constructor=skbio.DNA, format='fasta'))
        exp_seqs = [skbio.DNA('AAACGTTACGGTTAACTATACATGCAGAAGACTAATCGG',
                              metadata={'id': ('4574b947a0159c0da35a1f30f'
                                               '989681a1d9f64ef'),
                                        'description': 'sample1_1'}),
                    skbio.DNA('ACGTACGTACGTACGTACGTACGTACGTACGTGCATGGTGCGACCG',
                              metadata={'id': ('16a1263bde4f2f99422630d1bb'
                                               '87935c4236d1ba'),
                                        'description': 's2_42'}),
                    skbio.DNA('AAACGTTACGGTTAACTATACATGCAGAAGACTA',
                              metadata={'id': ('1768cf7fca79f84d651b34d878d'
                                               'e2492c6a7b971'),
                                        'description': 's2_2'})]
        self.assertEqual(obs_seqs, exp_seqs)

    def test_dereplicate_sequences_underscores_in_ids(self):
        input_sequences_fp = self.get_data_path('seqs-2')
        input_sequences = QIIME1DemuxDirFmt(input_sequences_fp, 'r')

        exp_table = biom.Table(np.array([[2, 1],
                                         [0, 1],
                                         [0, 2]]),
                               ['4574b947a0159c0da35a1f30f989681a1d9f64ef',
                                '1768cf7fca79f84d651b34d878de2492c6a7b971',
                                '16a1263bde4f2f99422630d1bb87935c4236d1ba'],
                               ['sa_mple1', 's2'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, obs_sequences = dereplicate_sequences(
                sequences=input_sequences)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        # sequences are reverse-sorted by abundance in output
        obs_seqs = list(skbio.io.read(str(obs_sequences),
                        constructor=skbio.DNA, format='fasta'))
        exp_seqs = [skbio.DNA('AAACGTTACGGTTAACTATACATGCAGAAGACTAATCGG',
                              metadata={'id': ('4574b947a0159c0da35a1f30f'
                                               '989681a1d9f64ef'),
                                        'description': 'sa_mple1_1'}),
                    skbio.DNA('ACGTACGTACGTACGTACGTACGTACGTACGTGCATGGTGCGACCG',
                              metadata={'id': ('16a1263bde4f2f99422630d1bb'
                                               '87935c4236d1ba'),
                                        'description': 's2_42'}),
                    skbio.DNA('AAACGTTACGGTTAACTATACATGCAGAAGACTA',
                              metadata={'id': ('1768cf7fca79f84d651b34d878d'
                                               'e2492c6a7b971'),
                                        'description': 's2_2'})]
        self.assertEqual(obs_seqs, exp_seqs)

    def test_dereplicate_sequences_prefix(self):
        input_sequences_fp = self.get_data_path('seqs-1')
        input_sequences = QIIME1DemuxDirFmt(input_sequences_fp, 'r')

        exp_table = biom.Table(np.array([[2, 2],
                                        [2, 0]]),
                               ['4574b947a0159c0da35a1f30f989681a1d9f64ef',
                                '16a1263bde4f2f99422630d1bb87935c4236d1ba'],
                               ['s2', 'sample1'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, obs_sequences = dereplicate_sequences(
                sequences=input_sequences, derep_prefix=True)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        # sequences are reverse-sorted by abundance in output
        obs_seqs = list(skbio.io.read(str(obs_sequences),
                        constructor=skbio.DNA, format='fasta'))
        exp_seqs = [skbio.DNA('AAACGTTACGGTTAACTATACATGCAGAAGACTAATCGG',
                              metadata={'id': ('4574b947a0159c0da35a1f30f'
                                               '989681a1d9f64ef'),
                                        'description': 's2_1'}),
                    skbio.DNA('ACGTACGTACGTACGTACGTACGTACGTACGTGCATGGTGCGACCG',
                              metadata={'id': ('16a1263bde4f2f99422630d1bb'
                                               '87935c4236d1ba'),
                                        'description': 's2_42'})]
        self.assertEqual(obs_seqs, exp_seqs)


class ParseUc(TestPluginBase):
    # These tests and the test data below them is copied from the biom-format
    # project temporarily to fix a bug in handling of sample ids with
    # underscores in them (https://github.com/biocore/biom-format/issues/758).
    # This code will be contribued back upstream to the
    # biom-format project, and will be removed from this plugin when a
    # biom-format release is available that contains this fix.

    package = 'q2_vsearch.tests'

    def test_empty(self):
        """ empty uc file returns empty Table
        """
        actual = _parse_uc(uc_empty.split('\n'))
        expected = biom.Table(np.array([[]]),
                              observation_ids=[],
                              sample_ids=[])
        self.assertEqual(actual, expected)

    def test_minimal(self):
        """ single new seed observed
        """
        actual = _parse_uc(uc_minimal.split('\n'))
        expected = biom.Table(np.array([[1.0]]),
                              observation_ids=['f2_1539'],
                              sample_ids=['f2'])
        self.assertEqual(actual, expected)

    def test_lib_minimal(self):
        """ single library seed observed
        """
        actual = _parse_uc(uc_lib_minimal.split('\n'))
        expected = biom.Table(np.array([[1.0]]),
                              observation_ids=['295053'],
                              sample_ids=['f2'])
        self.assertEqual(actual, expected)

    def test_invalid(self):
        """ invalid query sequence identifier detected
        """
        self.assertRaises(ValueError, _parse_uc, uc_invalid_id.split('\n'))

    def test_seed_hits(self):
        """ multiple new seeds observed
        """
        actual = _parse_uc(uc_seed_hits.split('\n'))
        expected = biom.Table(np.array([[2.0, 1.0], [0.0, 1.0]]),
                              observation_ids=['f2_1539', 'f3_44'],
                              sample_ids=['f2', 'f3'])
        self.assertEqual(actual, expected)

    def test_mixed_hits(self):
        """ new and library seeds observed
        """
        actual = _parse_uc(uc_mixed_hits.split('\n'))
        expected = biom.Table(np.array([[2.0, 1.0], [0.0, 1.0], [1.0, 0.0]]),
                              observation_ids=['f2_1539', 'f3_44', '295053'],
                              sample_ids=['f2', 'f3'])
        self.assertEqual(actual, expected)

    # the following tests are new with respect to biom-format project 2.1.6

    def test_underscore_in_sample_id(self):
        """ single new seed observed for sample with underscores in id
        """
        actual = _parse_uc(uc_minimal_w_underscores.split('\n'))
        expected = biom.Table(np.array([[1.0]]),
                              observation_ids=['sample_id_w_underscores_42'],
                              sample_ids=['sample_id_w_underscores'])
        print(actual)
        print(expected)
        self.assertEqual(actual, expected)

    def test_uc_w_comments_and_blank_lines(self):
        """ uc contains comments and blank lines
        """
        actual = _parse_uc(uc_w_comments_and_blank_lines.split('\n'))
        expected = biom.Table(np.array([[1.0]]),
                              observation_ids=['f2_1539'],
                              sample_ids=['f2'])
        self.assertEqual(actual, expected)


# no hits or library seeds
uc_empty = """
"""

# label not in qiime post-split-libraries format
uc_invalid_id = """
S	0	133	*	*	*	*	*	1539	*
"""

# contains single new (de novo) seed hit
uc_minimal = """
S	0	133	*	*	*	*	*	f2_1539	*
"""

# contains single new (de novo) seed hit
uc_w_comments_and_blank_lines = """# sdfsdfsdf


# dfasdfsdf
# sdffsdfsd

S	0	133	*	*	*	*	*	f2_1539	*

# sdfsdfsdfsdfsdfsdsdfsf
# asdasddpeanutdfdffsdfsdfsdfsdsdfsd sdfsdfsdf sdfdsf

"""

# contains single seed hit for a sample with underscores in its id
uc_minimal_w_underscores = """
S	0	133	*	*	*	*	*	sample_id_w_underscores_42	*
"""

# contains single library (reference) seed hit
uc_lib_minimal = """
L	3	1389	*	*	*	*	*	295053	*
H	3	133	100.0	+	0	0	519I133M737I	f2_1539	295053
"""

# contains new seed (de novo) hits only
uc_seed_hits = """
S	0	133	*	*	*	*	*	f2_1539	*
H	0	141	100.0	+	0	0	133M8D	f3_42	f2_1539
H	0	141	100.0	+	0	0	133M8D	f2_43	f2_1539
S	0	133	*	*	*	*	*	f3_44	*
"""

# contains library (reference) and new seed (de novo) hits
uc_mixed_hits = """
S	0	133	*	*	*	*	*	f2_1539	*
H	0	141	100.0	+	0	0	133M8D	f3_42	f2_1539
H	0	141	100.0	+	0	0	133M8D	f2_43	f2_1539
S	0	133	*	*	*	*	*	f3_44	*
L	3	1389	*	*	*	*	*	295053	*
H	3	133	100.0	+	0	0	519I133M737I	f2_1539	295053
"""
