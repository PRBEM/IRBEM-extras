#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# setup.py to install invlib Python interface

__author__ = 'Steve Morley, Los Alamos National Lab (smorley@lanl.gov)'

from distutils.core import setup
import os, sys, shutil

#test for python version 2.x where x>=5
try:
    dum = sys.version_info
    assert dum[0]==2
    assert dum[1]>=5
except:
    raise Exception("""PyInvlib requires Python 2.X, where X>=5.""")

try:
    shutil.copy('invlib.so','pyinvlib/invlib.so')
    shutil.copy('specinv_test.in1','pyinvlib/specinv_test.in1')
except:
    raise RuntimeError('invlib.so not found, please compile INVLIB')

pkg_files = ['invlib.so', 'specinv_test.in1']
# run setup from distutil
setup(name='pyinvlib',
      version='0.1',
      description='PyInvlib: Python interface to INVLIB',
      author='Steve Morley',
      author_email='smorley@lanl.gov',
      requires=['numpy'],
      packages=['pyinvlib'],
      package_dir={'pyinvlib': 'pyinvlib'},
      package_data={'pyinvlib': pkg_files}
      )

os.remove('pyinvlib/invlib.so')
os.remove('pyinvlib/specinv_test.in1')