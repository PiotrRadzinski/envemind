#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup, find_packages
from distutils.core import Extension
import pathlib
import opentimspy

here = pathlib.Path(__file__).parent.resolve()

setup(
    name='envemind',
    version='0.0.1',
    description='Prediction of monoisotopic mass in mass spectra',
#    long_description=(here / 'README.md').read_text(encoding='utf-8'),
#    long_description_content_type='text/markdown',
    url='https://github.com/PiotrRadzinski/envemind',
    author='Piotr Radziński, Michał Piotr Startek',
    classifiers=[
            'Development Status :: 2 - Pre-Alpha',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3 :: Only',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
            ],
    keywords = 'Mass spectrometry monisotopic mass',
    packages=find_packages(),
    python_requires='>=3.6',
    install_requires='numpy scipy IsoSpecPy pyteomics'.split(),
#    entry_points={},
#    scripts=[''],
)
