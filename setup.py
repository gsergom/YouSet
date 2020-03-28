#!/usr/bin/env python3
from distutils.core import setup

setup(
  name = 'YouSet',
  packages = ['YouSet'],
  version = '1.0',
  license='MIT',
  description = 'Superimposition macrocomplex builder',
  author = 'Paula Lopez, Gerard Serrano, Laura Vila',
  url = 'https://github.com/gsergom/YouSet',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/user/reponame/archive/v_01.tar.gz',    # I explain this later on
  keywords = ['macrocomplex', 'builder', 'bioinformatics'],
  install_requires=['biopython'])
