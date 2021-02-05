# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections

import pandas as pd
import numpy as np
import qiime2

from .plugin_setup import plugin
from ._format import UchimeStatsFmt

# many of the numeric fields will contain * if a query is
# not chimeric, so being safe and making all fields other than
# score strings
_uchime_stats_header = collections.OrderedDict([
     ('score', np.number),
     ('feature-id', np.str),
     ('A', np.str),
     ('B', np.str),
     ('T', np.str),
     ('idQM', np.str),
     ('idQA', np.str),
     ('idQB', np.str),
     ('idAB', np.str),
     ('idQT', np.str),
     ('LY', np.str),
     ('LN', np.str),
     ('LA', np.str),
     ('RY', np.str),
     ('RN', np.str),
     ('RA', np.str),
     ('div', np.str),
     ('YN', np.str)])


def _stats_to_df(ff):
    df = pd.read_csv(str(ff), sep='\t', index_col='feature-id',
                     names=_uchime_stats_header.keys(),
                     dtype=_uchime_stats_header)
    return df


@plugin.register_transformer
def _1(ff: UchimeStatsFmt) -> qiime2.Metadata:
    return qiime2.Metadata(_stats_to_df(ff))


@plugin.register_transformer
def _2(ff: UchimeStatsFmt) -> pd.DataFrame:
    return _stats_to_df(ff)
