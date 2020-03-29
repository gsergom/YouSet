#!/usr/bin/env python3
from distutils.core import setup

setup(
  name = 'YouSet',
  packages = ['YouSet','YouSet.modules'],
  version = '2.0',
  license='MIT',
  description = 'Superimposition macrocomplex builder',
  long_description = open('README.md').read(),
  long_description_content_type="text/markdown",
  author = 'Paula Lopez, Gerard Serrano, Laura Vila',
  url = 'https://github.com/gsergom/YouSet',
  download_url = 'https://github.com/gsergom/YouSet/releases/tag/2.0',
  keywords = ['macrocomplex', 'builder', 'bioinformatics'],
  install_requires=['biopython'])
