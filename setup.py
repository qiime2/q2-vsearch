# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

import versioneer

setup(
    name="q2-vsearch",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    description="QIIME 2 plugin for vsearch.",
    license='BSD-3-Clause',
    url="https://qiime2.org",
    entry_points={
        "qiime2.plugins":
        ["q2-vsearch=q2_vsearch.plugin_setup:plugin"]
    },
    package_data={
        'q2_vsearch': ['assets/*.html',
                       'citations.bib'],
        'q2_vsearch.tests': ['data/*',
                             'data/seqs-1/*',
                             'data/seqs-2/*',
                             'data/demux-1/*',
                             'data/demux-1_se/*',
                             'data/stats/*',
                             ]},
    zip_safe=False,
    )
