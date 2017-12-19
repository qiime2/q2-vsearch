# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import csv

import qiime2.plugin.model as model
from qiime2.plugin import ValidationError


class UchimeStatsFmt(model.TextFileFormat):
    def _validate_(self):
        with open(str(self)) as fh:
            csv_reader = csv.reader(fh)
            for fields in csv_reader:
                if len(fields) != 18:
                    raise ValidationError(
                        'Incorrect number of fields detected. Should be '
                        'exactly 18.')


UchimeStatsDirFmt = model.SingleFileDirectoryFormat(
    'UchimeStatsDirFmt', 'stats.tsv', UchimeStatsFmt)
