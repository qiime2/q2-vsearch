# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin.model as model


class UchimeStatsFmt(model.TextFileFormat):
    def sniff(self):
        with open(str(self)) as fh:
            header = fh.readline().strip().split('\t')

        return len(header) == 18


UchimeStatsDirFmt = model.SingleFileDirectoryFormat(
    'UchimeStatsDirFmt', 'stats.tsv', UchimeStatsFmt)
