# Copyright 2012-2017 CNRS-UM LIRMM, CNRS-AIST JRL
#
# This file is part of SpaceVecAlg.
#
# SpaceVecAlg is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SpaceVecAlg is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SpaceVecAlg.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
try:
  from setuptools import setup
  from setuptools import Extension
except ImportError:
  from distutils.core import setup
  from distutils.extension import Extension

from Cython.Build import cythonize

import hashlib
import os
import subprocess

win32_build = os.name == 'nt'

this_path  = os.path.dirname(os.path.realpath(__file__))
with open(this_path + '/sva/__init__.py', 'w') as fd:
    fd.write('from .sva import *\n')

sha512 = hashlib.sha512()
src_files = ['sva/sva.pyx', 'sva/c_sva.pxd', 'sva/sva.pxd', 'include/sva_wrapper.hpp']
src_files = [ '{}/{}'.format(this_path, f) for f in src_files ]
for f in src_files:
  chunk = 2**16
  with open(f, 'r') as fd:
    while True:
      data = fd.read(chunk)
      if data:
        sha512.update(data.encode('ascii'))
      else:
        break
version_hash = sha512.hexdigest()[:7]

class pkg_config(object):
  def __init__(self):
    self.compile_args = '@SVA_COMPILE_FLAGS@'.split(';')
    self.include_dirs = '@SVA_INCLUDE_DIRECTORIES@'.split(';')
    self.library_dirs = '@SVA_LINK_FLAGS@'.split(';')
    self.libraries = '@SVA_LINK_LIBRARIES@'.split(';')
    self.found = True

python_libs = []
python_lib_dirs = []
python_others = []
if not win32_build:
  tokens = subprocess.check_output(['python-config', '--ldflags']).split()
  tokens = [ token.decode('ascii') for token in tokens ]
  for token in tokens:
    flag = token[:2]
    value = token[2:]
    if flag == '-l':
      python_libs.append(value)
    elif flag == '-L':
      python_lib_dirs.append(value)
    elif token[:1] == '-':
      python_others.append(token)

config = pkg_config()

config.compile_args.append('-std=c++11')
for o in python_others:
  config.compile_args.append(o)
config.include_dirs.append(os.getcwd() + "/include")
if not win32_build:
  config.library_dirs.extend(python_lib_dirs)
  config.libraries.extend(python_libs)
else:
  config.compile_args.append("-DWIN32")

def GenExtension(name, pkg, ):
  pyx_src = name.replace('.', '/')
  cpp_src = pyx_src + '.cpp'
  pyx_src = pyx_src + '.pyx'
  ext_src = pyx_src
  if pkg.found:
    return Extension(name, [ext_src], extra_compile_args = pkg.compile_args, include_dirs = pkg.include_dirs, library_dirs = pkg.library_dirs, libraries = pkg.libraries)
  else:
    print("Failed to find {}".format(pkg.name))
    return None

extensions = [
  GenExtension('sva.sva', config)
]

extensions = [ x for x in extensions if x is not None ]
packages = ['sva']
data = ['__init__.py', 'c_sva.pxd', 'sva.pxd']

extensions = cythonize(extensions)

setup(
    name = 'sva',
    version='1.0.0-{}'.format(version_hash),
    ext_modules = extensions,
    packages = packages,
    package_data = { 'sva': data }
)