from setuptools import setup, find_packages
import re
import ast

# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('q2_vsearch/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    version = str(ast.literal_eval(hit))

setup(
    name="q2-vsearch",
    version=version,
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
    package_data={}
    )
