#!/usr/bin/env python

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

def get_version():
    with open(path.join(here, "trans_savnet/version.py"), encoding = 'utf-8') as hin:
        for line in hin:
            if line.startswith("__version__"):
                version = line.partition('=')[2]
                return version.strip().strip('\'"')
    raise ValueError('Could not find version.')

setup(
    name = 'savnet',
    version = get_version(),
    description='Python tools for evaluating trans-acting splicing changes caused by somatic mutations',
    url = 'https://github.com/friend1ws/TransSAVNet',
    author = 'Yuichi Shiraishi',
    author_email = 'friend1ws@gamil.com',
    license = 'GPLv3',

    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: Unix',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],

    packages = find_packages(exclude = ['tests']),

    install_requires = ["annot_utils>=0.2.1", "pysam>=0.9.0", "junc_utils>=0.4.1", "intron_retention_utils>=0.5.1", "savnet>=0.3.3b1"],
    entry_points = {'console_scripts': ['trans_savnet = trans_savnet:main']}

)


