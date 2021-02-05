# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import sqlite3
import tempfile

import skbio
import biom
import numpy as np

from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase
from qiime2.util import redirected_stdio
from q2_types.feature_data import DNAFASTAFormat, DNAIterator

from q2_vsearch._cluster_features import (cluster_features_de_novo,
                                          cluster_features_closed_reference,
                                          _fasta_with_sizes,
                                          _fasta_from_sqlite,
                                          VSearchError)


def _read_seqs(seqs):
    if isinstance(seqs, DNAFASTAFormat):
        read_seqs = skbio.io.read(str(seqs), constructor=skbio.DNA,
                                  format='fasta')
    elif isinstance(seqs, str):
        read_seqs = skbio.io.read(seqs, constructor=skbio.DNA, format='fasta')
    else:
        read_seqs = seqs.view(DNAIterator)
    return list(read_seqs)


def _relabel_seqs(seqs, labels):
    for i in range(len(seqs)):
        seqs[i].metadata['id'] = labels[i]


class ClusterFeaturesDenovoTests(TestPluginBase):

    package = 'q2_vsearch.tests'

    def setUp(self):
        super().setUp()
        input_sequences_fp = self.get_data_path('dna-sequences-1.fasta')
        self.input_sequences = DNAFASTAFormat(input_sequences_fp, mode='r')
        self.input_table = biom.Table(np.array([[100, 101, 103],
                                                [1, 1, 2],
                                                [4, 5, 6],
                                                [7, 8, 9]]),
                                      ['feature1', 'feature2', 'feature3',
                                       'feature4'],
                                      ['sample1', 'sample2', 'sample3'])
        self.input_sequences_list = _read_seqs(self.input_sequences)

    def test_no_clustering(self):
        with redirected_stdio(stderr=os.devnull):
            obs_table, obs_sequences = cluster_features_de_novo(
                sequences=self.input_sequences, table=self.input_table,
                perc_identity=1.0)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(self.input_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, self.input_table)

        obs_seqs = _read_seqs(obs_sequences)

        # sequences are reverse-sorted by abundance in output
        exp_seqs = [self.input_sequences_list[0], self.input_sequences_list[3],
                    self.input_sequences_list[2], self.input_sequences_list[1]]
        self.assertEqual(obs_seqs, exp_seqs)

    def test_99_percent_clustering(self):
        exp_table = biom.Table(np.array([[104, 106, 109],
                                         [1, 1, 2],
                                         [7, 8, 9]]),
                               ['feature1', 'feature2',
                                'feature4'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, obs_sequences = cluster_features_de_novo(
                sequences=self.input_sequences, table=self.input_table,
                perc_identity=0.99)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        # sequences are reverse-sorted by abundance in output
        obs_seqs = _read_seqs(obs_sequences)
        exp_seqs = [self.input_sequences_list[0], self.input_sequences_list[3],
                    self.input_sequences_list[1]]
        self.assertEqual(obs_seqs, exp_seqs)

    def test_97_percent_clustering_feature1_most_abundant(self):
        exp_table = biom.Table(np.array([[111, 114, 118],
                                         [1, 1, 2]]),
                               ['feature1', 'feature2'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, obs_sequences = cluster_features_de_novo(
                sequences=self.input_sequences, table=self.input_table,
                perc_identity=0.97)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        # sequences are reverse-sorted by abundance in output
        obs_seqs = _read_seqs(obs_sequences)
        exp_seqs = [self.input_sequences_list[0], self.input_sequences_list[1]]
        self.assertEqual(obs_seqs, exp_seqs)

    def test_97_percent_clustering_feature3_most_abundant(self):
        input_table = biom.Table(np.array([[4, 5, 6],
                                           [1, 1, 2],
                                           [100, 101, 103],
                                           [7, 8, 9]]),
                                 ['feature1', 'feature2', 'feature3',
                                  'feature4'],
                                 ['sample1', 'sample2', 'sample3'])
        exp_table = biom.Table(np.array([[111, 114, 118],
                                         [1, 1, 2]]),
                               ['feature3', 'feature2'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, obs_sequences = cluster_features_de_novo(
                sequences=self.input_sequences, table=input_table,
                perc_identity=0.97)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        # sequences are reverse-sorted by abundance in output
        obs_seqs = _read_seqs(obs_sequences)
        exp_seqs = [self.input_sequences_list[2], self.input_sequences_list[1]]
        self.assertEqual(obs_seqs, exp_seqs)

    def test_97_percent_clustering_feature4_most_abundant(self):
        input_table = biom.Table(np.array([[4, 5, 6],
                                           [1, 1, 2],
                                           [7, 8, 9],
                                           [100, 101, 103]]),
                                 ['feature1', 'feature2', 'feature3',
                                  'feature4'],
                                 ['sample1', 'sample2', 'sample3'])
        exp_table = biom.Table(np.array([[111, 114, 118],
                                         [1, 1, 2]]),
                               ['feature4', 'feature2'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, obs_sequences = cluster_features_de_novo(
                sequences=self.input_sequences, table=input_table,
                perc_identity=0.97)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        # sequences are reverse-sorted by abundance in output
        obs_seqs = _read_seqs(obs_sequences)
        exp_seqs = [self.input_sequences_list[3], self.input_sequences_list[1]]
        self.assertEqual(obs_seqs, exp_seqs)

    def test_1_percent_clustering(self):
        exp_table = biom.Table(np.array([[112, 115, 120]]),
                               ['feature1'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, obs_sequences = cluster_features_de_novo(
                sequences=self.input_sequences, table=self.input_table,
                perc_identity=0.01)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        # sequences are reverse-sorted by abundance in output
        obs_seqs = _read_seqs(obs_sequences)
        exp_seqs = [self.input_sequences_list[0]]
        self.assertEqual(obs_seqs, exp_seqs)

    def test_short_sequences(self):
        input_sequences_fp = self.get_data_path('dna-sequences-short.fasta')
        input_sequences = DNAFASTAFormat(input_sequences_fp, mode='r')

        input_table = biom.Table(np.array([[0, 1, 3],
                                           [1, 1, 2]]),
                                 ['feature1', 'feature2'],
                                 ['sample1', 'sample2', 'sample3'])

        exp_table = biom.Table(np.array([[1, 2, 5]]),
                               ['feature1'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, obs_sequences = cluster_features_de_novo(
                sequences=input_sequences, table=input_table,
                perc_identity=0.01)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')

    def test_extra_features_in_sequences(self):
        input_table = biom.Table(np.array([[0, 1, 3], [1, 1, 2], [4, 5, 6]]),
                                 ['feature1', 'feature2', 'feature3'],
                                 ['sample1', 'sample2', 'sample3'])
        with self.assertRaisesRegex(ValueError,
                                    expected_regex='Feature feature4 is pre'):
            clustered_table, clustered_sequences = cluster_features_de_novo(
                sequences=self.input_sequences, table=input_table,
                perc_identity=1.0)

    def test_extra_features_in_table(self):
        input_table = biom.Table(np.array([[0, 1, 3],
                                           [1, 1, 2],
                                           [4, 5, 6],
                                           [7, 8, 9],
                                           [1, 1, 1]]),
                                 ['feature1', 'feature2', 'feature3',
                                  'feature4', 'feature5'],
                                 ['sample1', 'sample2', 'sample3'])
        with self.assertRaisesRegex(ValueError,
                                    expected_regex='Some feat.*feature5.*'):
            clustered_table, clustered_sequences = cluster_features_de_novo(
                sequences=self.input_sequences, table=input_table,
                perc_identity=1.0)

    def test_no_overlapping_feature_ids(self):
        input_table = biom.Table(np.array([[0, 1, 3],
                                           [1, 1, 2],
                                           [4, 5, 6],
                                           [7, 8, 9],
                                           [1, 1, 1]]),
                                 ['f1', 'f2', 'f3',
                                  'f4', 'f5'],
                                 ['sample1', 'sample2', 'sample3'])
        with self.assertRaisesRegex(ValueError,
                                    expected_regex='Feature feature1 is pre'):
            clustered_table, clustered_sequences = cluster_features_de_novo(
                sequences=self.input_sequences, table=input_table,
                perc_identity=1.0)


class ClusterFeaturesClosedReference(TestPluginBase):

    package = 'q2_vsearch.tests'

    def setUp(self):
        super().setUp()
        input_sequences_fp = self.get_data_path('dna-sequences-1.fasta')
        self.input_sequences = DNAFASTAFormat(input_sequences_fp, mode='r')
        ref_sequences_1_fp = self.get_data_path('ref-sequences-1.fasta')
        self.ref_sequences_1 = DNAFASTAFormat(ref_sequences_1_fp, mode='r')
        ref_sequences_2_fp = self.get_data_path('ref-sequences-2.fasta')
        self.ref_sequences_2 = DNAFASTAFormat(ref_sequences_2_fp, mode='r')
        self.input_table = biom.Table(np.array([[100, 101, 103],
                                                [1, 1, 2],
                                                [4, 5, 6],
                                                [7, 8, 9]]),
                                      ['feature1', 'feature2', 'feature3',
                                       'feature4'],
                                      ['sample1', 'sample2', 'sample3'])
        self.input_sequences_list = _read_seqs(self.input_sequences)

    def test_100_percent_clustering(self):
        # feature2 and feature3 don't cluster
        exp_table = biom.Table(np.array([[100, 101, 103],
                                         [7, 8, 9]]),
                               ['r1', 'r2'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, matched_seqs, unmatched_seqs = \
                    cluster_features_closed_reference(
                        sequences=self.input_sequences, table=self.input_table,
                        reference_sequences=self.ref_sequences_1,
                        perc_identity=1.0)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        obs_matched_seqs = _read_seqs(matched_seqs)
        # The rep seqs selected are feature1 and feature4, for r1 and r2,
        # respectively. Since no other features are in the cluster, there is
        # no count-based selection of the rep seq.
        exp_matched_seqs = [self.input_sequences_list[0],  # feature1
                            self.input_sequences_list[3]]  # feature4
        _relabel_seqs(exp_matched_seqs, ['r1', 'r2'])
        self.assertEqual(obs_matched_seqs, exp_matched_seqs)

        obs_unmatched_seqs = _read_seqs(unmatched_seqs)
        exp_unmatched_seqs = [self.input_sequences_list[2],  # feature3
                              self.input_sequences_list[1]]  # feature2
        self.assertEqual(obs_unmatched_seqs, exp_unmatched_seqs)

    def test_100_percent_clustering_strand(self):
        # feature2 and feature3 don't cluster
        exp_table = biom.Table(np.array([[100, 101, 103],
                                         [7, 8, 9]]),
                               ['r1', 'r2'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, matched_seqs, unmatched_seqs = \
                    cluster_features_closed_reference(
                        sequences=self.input_sequences, table=self.input_table,
                        reference_sequences=self.ref_sequences_2,
                        perc_identity=1.0, strand='both')
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        obs_matched_seqs = _read_seqs(matched_seqs)
        # The rep seqs selected are feature1 and feature4, for r1 and r2,
        # respectively. Since no other features are in the cluster, there is
        # no count-based selection of the rep seq.
        exp_matched_seqs = [self.input_sequences_list[0],  # feature1
                            self.input_sequences_list[3]]  # feature4
        _relabel_seqs(exp_matched_seqs, ['r1', 'r2'])
        self.assertEqual(obs_matched_seqs, exp_matched_seqs)

        obs_unmatched_seqs = _read_seqs(unmatched_seqs)
        exp_unmatched_seqs = [self.input_sequences_list[2],  # feature3
                              self.input_sequences_list[1]]  # feature2
        self.assertEqual(obs_unmatched_seqs, exp_unmatched_seqs)

    def test_no_matches(self):
        with self.assertRaisesRegex(VSearchError,
                                    expected_regex='No matches were iden'):
            with redirected_stdio(stderr=os.devnull):
                # self.ref_sequences_2 are rev comps of self.ref_sequences_1,
                # so if strand='both' is not passed, there should be no matches
                cluster_features_closed_reference(
                    sequences=self.input_sequences, table=self.input_table,
                    reference_sequences=self.ref_sequences_2,
                    perc_identity=1.0)

    def test_99_percent_clustering(self):
        # feature1 and feature3 cluster together; feature2 doesn't cluster at
        # all; feature4 clusters alone.
        exp_table = biom.Table(np.array([[104, 106, 109],
                                         [7, 8, 9]]),
                               ['r1', 'r2'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, matched_seqs, unmatched_seqs = \
                    cluster_features_closed_reference(
                        sequences=self.input_sequences, table=self.input_table,
                        reference_sequences=self.ref_sequences_1,
                        perc_identity=0.99)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        obs_matched_seqs = _read_seqs(matched_seqs)
        # The rep seqs selected are feature1 and feature4, for r1 and r2,
        # respectively. feature1 and feature3 are in the same cluster, but
        # feature1 is selected as the rep seq because it has a higher count.
        exp_matched_seqs = [self.input_sequences_list[0],  # feature1
                            self.input_sequences_list[3]]  # feature4
        _relabel_seqs(exp_matched_seqs, ['r1', 'r2'])
        self.assertEqual(obs_matched_seqs, exp_matched_seqs)

        obs_unmatched_seqs = _read_seqs(unmatched_seqs)
        exp_unmatched_seqs = [self.input_sequences_list[1]]  # feature2
        self.assertEqual(obs_unmatched_seqs, exp_unmatched_seqs)

    def test_97_percent_clustering(self):
        # feature1 and feature3 cluster together; feature2 doesn't cluster at
        # all; feature 4 clusters alone.
        exp_table = biom.Table(np.array([[104, 106, 109],
                                         [7, 8, 9]]),
                               ['r1', 'r2'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, matched_seqs, unmatched_seqs = \
                    cluster_features_closed_reference(
                        sequences=self.input_sequences, table=self.input_table,
                        reference_sequences=self.ref_sequences_1,
                        perc_identity=0.97)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        obs_matched_seqs = _read_seqs(matched_seqs)
        # The rep seqs selected are feature1 and feature4, for r1 and r2,
        # respectively. feature1 and feature3 are in the same cluster, but
        # feature1 is selected as the rep seq because it has a higher count.
        exp_matched_seqs = [self.input_sequences_list[0],  # feature1
                            self.input_sequences_list[3]]  # feature4
        _relabel_seqs(exp_matched_seqs, ['r1', 'r2'])
        self.assertEqual(obs_matched_seqs, exp_matched_seqs)

        obs_unmatched_seqs = _read_seqs(unmatched_seqs)
        exp_unmatched_seqs = [self.input_sequences_list[1]]  # feature2
        self.assertEqual(obs_unmatched_seqs, exp_unmatched_seqs)

    def test_1_percent_clustering(self):
        # feature1 and feature3 cluster together; feature2 and feature4
        # cluster together;
        exp_table = biom.Table(np.array([[104, 106, 109],
                                         [8, 9, 11]]),
                               ['r1', 'r2'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, matched_seqs, unmatched_seqs = \
                    cluster_features_closed_reference(
                        sequences=self.input_sequences, table=self.input_table,
                        reference_sequences=self.ref_sequences_1,
                        perc_identity=0.01)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        obs_matched_seqs = _read_seqs(matched_seqs)
        # The rep seqs selected are feature1 and feature4, for r1 and r2,
        # respectively. feature1 and feature3 are in the same cluster, but
        # feature1 is selected as the rep seq because it has a higher count.
        # Similarly, feature4 is selected as the cluster rep seq  because it
        # has a higher count.
        exp_matched_seqs = [self.input_sequences_list[0],  # feature1
                            self.input_sequences_list[3]]  # feature4
        _relabel_seqs(exp_matched_seqs, ['r1', 'r2'])
        self.assertEqual(obs_matched_seqs, exp_matched_seqs)

        # all sequences matched, so unmatched seqs is empty
        self.assertEqual(os.path.getsize(str(unmatched_seqs)), 0)

    def test_1_percent_clustering_alt_abundances(self):
        # feature1 and feature3 cluster together; feature2 and feature4
        # cluster together;
        input_table = biom.Table(np.array([[4, 5, 6],
                                           [7, 8, 9],
                                           [100, 101, 103],
                                           [1, 1, 2]]),
                                 ['feature1', 'feature2', 'feature3',
                                  'feature4'],
                                 ['sample1', 'sample2', 'sample3'])
        exp_table = biom.Table(np.array([[104, 106, 109],
                                         [8, 9, 11]]),
                               ['r1', 'r2'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, matched_seqs, unmatched_seqs = \
                    cluster_features_closed_reference(
                        sequences=self.input_sequences, table=input_table,
                        reference_sequences=self.ref_sequences_1,
                        perc_identity=0.01)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        obs_matched_seqs = _read_seqs(matched_seqs)
        # The rep seqs selected are feature3 and feature2, for r1 and r2,
        # respectively. feature1 and feature3 are in the same cluster, but
        # feature3 is selected as the rep seq because it has a higher count.
        # Similarly, feature2 is selected as the cluster rep seq  because it
        # has a higher count.
        exp_matched_seqs = [self.input_sequences_list[2],  # feature3
                            self.input_sequences_list[1]]  # feature2
        _relabel_seqs(exp_matched_seqs, ['r1', 'r2'])
        self.assertEqual(obs_matched_seqs, exp_matched_seqs)

        # all sequences matched, so unmatched seqs is empty
        self.assertEqual(os.path.getsize(str(unmatched_seqs)), 0)

    def test_short_sequences(self):
        input_sequences_fp = self.get_data_path('dna-sequences-short.fasta')
        input_sequences = DNAFASTAFormat(input_sequences_fp, mode='r')

        input_table = biom.Table(np.array([[0, 1, 3],
                                           [1, 1, 2]]),
                                 ['feature1', 'feature2'],
                                 ['sample1', 'sample2', 'sample3'])

        exp_table = biom.Table(np.array([[1, 2, 5]]),
                               ['r2'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, matched_seqs, unmatched_seqs = \
                cluster_features_closed_reference(
                    sequences=input_sequences, table=input_table,
                    reference_sequences=self.ref_sequences_1,
                    perc_identity=0.01)

        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')

    def test_extra_features_in_sequences(self):
        input_table = biom.Table(np.array([[0, 1, 3], [1, 1, 2], [4, 5, 6]]),
                                 ['feature1', 'feature2', 'feature3'],
                                 ['sample1', 'sample2', 'sample3'])
        with self.assertRaisesRegex(ValueError,
                                    expected_regex='Some feat.*feature4.*'):
            cluster_features_closed_reference(
                sequences=self.input_sequences, table=input_table,
                reference_sequences=self.ref_sequences_1, perc_identity=1.0)

    def test_extra_features_in_table(self):
        input_table = biom.Table(np.array([[0, 1, 3],
                                           [1, 1, 2],
                                           [4, 5, 6],
                                           [7, 8, 9],
                                           [1, 1, 1]]),
                                 ['feature1', 'feature2', 'feature3',
                                  'feature4', 'feature5'],
                                 ['sample1', 'sample2', 'sample3'])
        with self.assertRaisesRegex(ValueError,
                                    expected_regex='Some feat.*feature5.*'):
            cluster_features_closed_reference(
                sequences=self.input_sequences, table=input_table,
                reference_sequences=self.ref_sequences_1, perc_identity=1.0)

    def test_no_overlapping_feature_ids(self):
        input_table = biom.Table(np.array([[0, 1, 3],
                                           [1, 1, 2],
                                           [4, 5, 6],
                                           [7, 8, 9],
                                           [1, 1, 1]]),
                                 ['f1', 'f2', 'f3',
                                  'f4', 'f5'],
                                 ['sample1', 'sample2', 'sample3'])
        with self.assertRaisesRegex(ValueError,
                                    expected_regex='Some feat.*f1.*'):
            cluster_features_closed_reference(
                sequences=self.input_sequences, table=input_table,
                reference_sequences=self.ref_sequences_1, perc_identity=1.0)

    def test_features_with_same_counts(self):
        # feature1 and feature3 cluster into r1, feature2 and feature4 cluster
        # into r2. The features within a cluster have the same count, so this
        # test should ensure that the right rep seq is picked for each cluster.
        # The query in _fasta_from_sqlite should break ties by using the
        # first feature when sorting the tied features alphabetically by id.
        input_table = biom.Table(np.array([[4, 5, 6],
                                           [1, 2, 3],
                                           [4, 6, 5],
                                           [2, 1, 3]]),
                                 ['feature1', 'feature2', 'feature3',
                                  'feature4'],
                                 ['sample1', 'sample2', 'sample3'])
        exp_table = biom.Table(np.array([[8, 11, 11],
                                         [3, 3, 6]]),
                               ['r1', 'r2'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, matched_seqs, unmatched_seqs = \
                    cluster_features_closed_reference(
                        sequences=self.input_sequences, table=input_table,
                        reference_sequences=self.ref_sequences_1,
                        perc_identity=0.01)
        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        obs_matched_seqs = _read_seqs(matched_seqs)
        # The rep seqs selected are feature1 and feature2, for r1 and r2,
        # respectively. feature1 and feature3 are in the same cluster, but
        # feature1 is selected as the rep seq because it comes first
        # alphabetically, breaking the tie caused by the same counts.
        # Similarly, feature2 is selected as the cluster rep seq  because it
        # has a higher count.
        exp_matched_seqs = [self.input_sequences_list[0],  # feature1
                            self.input_sequences_list[1]]  # feature2
        _relabel_seqs(exp_matched_seqs, ['r1', 'r2'])
        self.assertEqual(obs_matched_seqs, exp_matched_seqs)

        # all sequences matched, so unmatched seqs is empty
        self.assertEqual(os.path.getsize(str(unmatched_seqs)), 0)


class ClusterFeaturesOpenReference(TestPluginBase):

    package = 'q2_vsearch.tests'

    def setUp(self):
        super().setUp()
        self.open_reference = \
            self.plugin.pipelines['cluster_features_open_reference']

        input_sequences_fp = self.get_data_path('dna-sequences-2.fasta')
        self.input_sequences = Artifact.import_data('FeatureData[Sequence]',
                                                    input_sequences_fp)

        ref_sequences_fp = self.get_data_path('ref-sequences-1.fasta')
        self.ref_sequences = Artifact.import_data('FeatureData[Sequence]',
                                                  ref_sequences_fp)

        input_table = biom.Table(np.array([[100, 101, 103],
                                           [1, 1, 2],
                                           [4, 5, 6],
                                           [7, 8, 9],
                                           [99, 98, 97]]),
                                 ['feature1', 'feature2', 'feature3',
                                  'feature4', 'feature5'],
                                 ['sample1', 'sample2', 'sample3'])
        self.input_table = Artifact.import_data('FeatureTable[Frequency]',
                                                input_table)
        self.input_sequences_list = _read_seqs(self.input_sequences)

    def test_100_percent_clustering(self):
        # feature1 clusters into r1 and feature4 clusters into r4 during
        # closed-ref clustering; feature2, feature3, and feature5 cluster into
        # their own clusters during de-novo clustering.
        exp_table = biom.Table(np.array([[100, 101, 103],
                                         [1, 1, 2],
                                         [4, 5, 6],
                                         [7, 8, 9],
                                         [99, 98, 97]]),
                               ['r1', 'feature2', 'feature3', 'r2',
                                'feature5'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, rep_seqs, new_ref_seqs = self.open_reference(
                sequences=self.input_sequences, table=self.input_table,
                reference_sequences=self.ref_sequences, perc_identity=1.0)

        obs_table = obs_table.view(biom.Table)
        obs_table_ids = set(obs_table.ids(axis='observation'))
        exp_table_ids = set(exp_table.ids(axis='observation'))
        self.assertEqual(obs_table_ids, exp_table_ids)

        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        obs_rep_seqs = _read_seqs(rep_seqs)
        exp_rep_seqs = [self.input_sequences_list[1],  # feature2
                        self.input_sequences_list[2],  # feature3
                        self.input_sequences_list[4],  # feature5
                        self.input_sequences_list[0],  # feature1
                        self.input_sequences_list[3]]  # feature4
        _relabel_seqs(exp_rep_seqs, ['feature2', 'feature3', 'feature5',
                                     'r1', 'r2'])
        self.assertEqual(obs_rep_seqs, exp_rep_seqs)

        obs_ref_seqs = _read_seqs(new_ref_seqs)
        ref_seqs = _read_seqs(self.ref_sequences)
        exp_ref_seqs = [self.input_sequences_list[1],  # feature2
                        self.input_sequences_list[2],
                        self.input_sequences_list[4]]  # feature3
        exp_ref_seqs = exp_ref_seqs + ref_seqs  # r1, r2
        self.assertEqual(obs_ref_seqs, exp_ref_seqs)

    def test_97_percent_clustering(self):
        # feature1 and feature3 clusters into r1 and feature4 clusters into r4
        # during closed-ref clustering; feature2 and feature5 clusters into a
        # cluster during de-novo clustering.
        exp_table = biom.Table(np.array([[104, 106, 109],
                                         [7, 8, 9],
                                         [100, 99, 99]]),
                               ['r1', 'r2', 'feature5'],
                               ['sample1', 'sample2', 'sample3'])

        with redirected_stdio(stderr=os.devnull):
            obs_table, rep_seqs, new_ref_seqs = self.open_reference(
                sequences=self.input_sequences, table=self.input_table,
                reference_sequences=self.ref_sequences, perc_identity=0.97)

        obs_table = obs_table.view(biom.Table)
        obs_table_ids = set(obs_table.ids(axis='observation'))
        exp_table_ids = set(exp_table.ids(axis='observation'))
        self.assertEqual(obs_table_ids, exp_table_ids)

        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        obs_rep_seqs = _read_seqs(rep_seqs)
        exp_rep_seqs = [self.input_sequences_list[4],  # feature5
                        self.input_sequences_list[0],  # feature1
                        self.input_sequences_list[3]]  # feature4
        _relabel_seqs(exp_rep_seqs, ['feature5', 'r1', 'r2'])
        self.assertEqual(obs_rep_seqs, exp_rep_seqs)

        obs_ref_seqs = _read_seqs(new_ref_seqs)
        ref_seqs = _read_seqs(self.ref_sequences)
        exp_ref_seqs = [self.input_sequences_list[4]]  # feature5
        exp_ref_seqs = exp_ref_seqs + ref_seqs  # r1, r2
        self.assertEqual(obs_ref_seqs, exp_ref_seqs)

    def test_skip_denovo(self):
        # feature1 and feature3 clusters into r1 and feature2 and feature4
        # clusters into r2 during closed-ref clustering; no unclustered
        # features so de-novo clustering is skipped.
        exp_table = biom.Table(np.array([[104, 106, 109],
                                         [107, 107, 108]]),
                               ['r1', 'r2'],
                               ['sample1', 'sample2', 'sample3'])
        with redirected_stdio(stderr=os.devnull):
            obs_table, rep_seqs, new_ref_seqs = self.open_reference(
                sequences=self.input_sequences, table=self.input_table,
                reference_sequences=self.ref_sequences, perc_identity=0.01)

        obs_table = obs_table.view(biom.Table)
        obs_table_ids = set(obs_table.ids(axis='observation'))
        exp_table_ids = set(exp_table.ids(axis='observation'))
        self.assertEqual(obs_table_ids, exp_table_ids)

        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        obs_rep_seqs = _read_seqs(rep_seqs)
        exp_rep_seqs = [self.input_sequences_list[0],  # feature1
                        self.input_sequences_list[4]]  # feature5
        _relabel_seqs(exp_rep_seqs, ['r1', 'r2'])
        self.assertEqual(obs_rep_seqs, exp_rep_seqs)

        obs_ref_seqs = _read_seqs(new_ref_seqs)
        # The returned "new" ref seqs should be the same as the original ref
        # seqs, because we skipped de-novo clustering.
        exp_ref_seqs = _read_seqs(self.ref_sequences)
        self.assertEqual(obs_ref_seqs, exp_ref_seqs)

    def test_skip_closed_reference(self):
        # feature1 and feature3 clusters into r1 and feature2 and feature4
        # clusters into r2 during closed-ref clustering; no unclustered
        # features so de-novo clustering is skipped.

        exp_table = biom.Table(np.array([[100, 101, 103],
                                         [1, 1, 2],
                                         [4, 5, 6],
                                         [7, 8, 9],
                                         [99, 98, 97]]),
                               ['feature1', 'feature2', 'feature3',
                                'feature4', 'feature5'],
                               ['sample1', 'sample2', 'sample3'])
        ref_sequences_fp = self.get_data_path('ref-sequences-2.fasta')
        ref_sequences = Artifact.import_data('FeatureData[Sequence]',
                                             ref_sequences_fp)

        with redirected_stdio(stderr=os.devnull):
            obs_table, rep_seqs, new_ref_seqs = self.open_reference(
                sequences=self.input_sequences, table=self.input_table,
                reference_sequences=ref_sequences, perc_identity=1.0)

        obs_table = obs_table.view(biom.Table)
        obs_table_ids = set(obs_table.ids(axis='observation'))
        exp_table_ids = set(exp_table.ids(axis='observation'))
        self.assertEqual(obs_table_ids, exp_table_ids)

        # order of identifiers is important for biom.Table equality
        obs_table = \
            obs_table.sort_order(exp_table.ids(axis='observation'),
                                 axis='observation')
        self.assertEqual(obs_table, exp_table)

        obs_rep_seqs = _read_seqs(rep_seqs)
        exp_rep_seqs = [self.input_sequences_list[0],  # feature1
                        self.input_sequences_list[4],  # feature5
                        self.input_sequences_list[3],  # feature4
                        self.input_sequences_list[2],  # feature3
                        self.input_sequences_list[1]]  # feature2
        self.assertEqual(obs_rep_seqs, exp_rep_seqs)

        obs_ref_seqs = _read_seqs(new_ref_seqs)
        # The returned "new" ref seqs should be the same as the original ref
        # seqs, plus the original rep seqs (since de-novo clustering happened
        # at 100% identity).
        exp_ref_seqs = self.input_sequences_list + _read_seqs(ref_sequences)
        self.assertEqual(obs_ref_seqs, exp_ref_seqs)


class PrivateFunctionTests(TestPluginBase):

    package = 'q2_vsearch.tests'

    def setUp(self):
        super().setUp()
        self.input_sequences_fp = self.get_data_path('dna-sequences-1.fasta')
        self.input_sequences = Artifact.import_data('FeatureData[Sequence]',
                                                    self.input_sequences_fp)
        self.input_table = biom.Table(np.array([[100, 101, 103],
                                                [1, 1, 2],
                                                [4, 5, 6],
                                                [7, 8, 9]]),
                                      ['feature1', 'feature2', 'feature3',
                                       'feature4'],
                                      ['sample1', 'sample2', 'sample3'])

    def test_fasta_with_sizes(self):
        with tempfile.NamedTemporaryFile() as output_sequences_f:
            _fasta_with_sizes(self.input_sequences_fp,
                              output_sequences_f.name,
                              self.input_table)

            obs_seqs = _read_seqs(output_sequences_f.name)
            input_seqs = _read_seqs(self.input_sequences_fp)

            self.assertEqual(len(obs_seqs), len(input_seqs))

            self.assertEqual(obs_seqs[0].metadata['id'], 'feature1;size=304')
            self.assertEqual(str(obs_seqs[0]), str(input_seqs[0]))
            self.assertEqual(obs_seqs[1].metadata['id'], 'feature2;size=4')
            self.assertEqual(str(obs_seqs[1]), str(input_seqs[1]))
            self.assertEqual(obs_seqs[2].metadata['id'], 'feature3;size=15')
            self.assertEqual(str(obs_seqs[2]), str(input_seqs[2]))
            self.assertEqual(obs_seqs[3].metadata['id'], 'feature4;size=24')
            self.assertEqual(str(obs_seqs[3]), str(input_seqs[3]))

    def test_fasta_from_sqlite(self):
        # artificially clustering feature1 and feature3 into r1, and
        # feature2 and feature4 into r2.
        conn = sqlite3.connect(':memory:')
        c = conn.cursor()
        c.execute('CREATE TABLE feature_cluster_map'
                  '(feature_id TEXT PRIMARY KEY,cluster_id TEXT NOT NULL, '
                  'count INTEGER);')
        s = [('feature1', 'r1', 204), ('feature2', 'r2', 4),
             ('feature3', 'r1', 15), ('feature4', 'r2', 24)]
        c.executemany('INSERT INTO feature_cluster_map VALUES (?, ?, ?);', s)
        conn.commit()

        with tempfile.NamedTemporaryFile() as output_sequences_f:
            _fasta_from_sqlite(conn, self.input_sequences_fp,
                               output_sequences_f.name)

            obs_seqs = _read_seqs(output_sequences_f.name)
        rep_seqs = _read_seqs(self.input_sequences)
        exp_seqs = [rep_seqs[0],  # feature1
                    rep_seqs[3]]  # feature4
        _relabel_seqs(exp_seqs, ['r1', 'r2'])
        self.assertEqual(obs_seqs, exp_seqs)

    def test_fasta_from_sqlite_same_clusters_different_rep_seq(self):
        # same as `test_fasta_from_sqlite`, but invert the selected rep seq
        # in each cluster by swapping the feature counts in the test data.
        conn = sqlite3.connect(':memory:')
        c = conn.cursor()
        c.execute('CREATE TABLE feature_cluster_map'
                  '(feature_id TEXT PRIMARY KEY,cluster_id TEXT NOT NULL, '
                  'count INTEGER);')
        s = [('feature1', 'r1', 15), ('feature2', 'r2', 24),
             ('feature3', 'r1', 204), ('feature4', 'r2', 4)]
        c.executemany('INSERT INTO feature_cluster_map VALUES (?, ?, ?);', s)
        conn.commit()

        with tempfile.NamedTemporaryFile() as output_sequences_f:
            _fasta_from_sqlite(conn, self.input_sequences_fp,
                               output_sequences_f.name)

            obs_seqs = _read_seqs(output_sequences_f.name)
        rep_seqs = _read_seqs(self.input_sequences)
        exp_seqs = [rep_seqs[2],  # feature3
                    rep_seqs[1]]  # feature2
        _relabel_seqs(exp_seqs, ['r1', 'r2'])
        self.assertEqual(obs_seqs, exp_seqs)

    def test_clusters_with_multiple_features_with_same_count(self):
        # feature1 and feature3 cluster into r1, feature2 and feature4 cluster
        # into r2. The features within a cluster have the same count, so this
        # test should ensure that the right rep seq is picked for each cluster.
        # The query in _fasta_from_sqlite should break ties by using the
        # first feature when sorting the tied features alphabetically by id.
        conn = sqlite3.connect(':memory:')
        c = conn.cursor()
        c.execute('CREATE TABLE feature_cluster_map'
                  '(feature_id TEXT PRIMARY KEY,cluster_id TEXT NOT NULL, '
                  'count INTEGER);')
        s = [('feature1', 'r1', 204), ('feature2', 'r2', 4),
             ('feature3', 'r1', 204), ('feature4', 'r2', 4)]
        c.executemany('INSERT INTO feature_cluster_map VALUES (?, ?, ?);', s)
        conn.commit()

        with tempfile.NamedTemporaryFile() as output_sequences_f:
            _fasta_from_sqlite(conn, self.input_sequences_fp,
                               output_sequences_f.name)

            obs_seqs = _read_seqs(output_sequences_f.name)
        rep_seqs = _read_seqs(self.input_sequences)
        exp_seqs = [rep_seqs[0],  # feature1
                    rep_seqs[1]]  # feature2
        _relabel_seqs(exp_seqs, ['r1', 'r2'])
        self.assertEqual(obs_seqs, exp_seqs)
