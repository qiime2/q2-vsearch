# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import qiime2

from .plugin_setup import plugin
from ._format import UchimeStatsFmt


_uchime_stats_header = ['score', 'Feature ID', 'A', 'B', 'T', 'idQM', 'idQA',
                        'idQB', 'idAB', 'idQT', 'LY', 'LN', 'LA', 'RY', 'RN',
                        'RA', 'div', 'YN']


def _stats_to_df(ff):
    df = pd.read_csv(str(ff), sep='\t')
    df.columns = _uchime_stats_header
    df = df.set_index('Feature ID')
    return df


@plugin.register_transformer
def _1(ff: UchimeStatsFmt) -> qiime2.Metadata:
    return qiime2.Metadata(_stats_to_df(ff))


@plugin.register_transformer
def _2(ff: UchimeStatsFmt) -> pd.DataFrame:
    return _stats_to_df(ff)
