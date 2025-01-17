# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase


class TestUsageExample(TestPluginBase):
    package = 'q2_vsearch.tests'

    def test_usage_cluster_features_de_novo(self):
        self.execute_examples()
