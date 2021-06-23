#!/usr/bin/env python3

from distutils.core import setup

setup(
    name='MERFISH_probe_design',
    version='0.0.0',
    author='Xingjie Pan',
    author_email='xingjiepan@gmail.com',
    url='https://github.com/xingjiepan/MERFISH_probe_design',
    packages=[
        'MERFISH_probe_design'    
    ],
    install_requires=[
        'numpy',
        'matplotlib',
        'seaborn',
        'biopython',
    ],
    description='A python3 pipeline for MERFISH probe design.',
    long_description=open('README.md').read(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
    ],
)
