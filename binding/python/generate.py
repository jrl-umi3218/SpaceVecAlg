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

from pybindgen import *
import sys

def import_eigen3_types(mod):
  mod.add_class('Vector3d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('Vector6d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

  mod.add_class('Matrix3d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('Matrix6d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')


def add_motion_vec(mod):
  mv = mod.add_class('MotionVec')
  mv.add_constructor([])
  mv.add_constructor([param('Eigen::Vector6d', 'vec')])
  mv.add_constructor([param('Eigen::Vector3d', 'angular'),
                      param('Eigen::Vector3d', 'linear')])

  mv.add_method('angular', retval('Eigen::Vector3d'), [], is_const=True)
  mv.add_method('linear', retval('Eigen::Vector3d'), [], is_const=True)
  mv.add_method('vector', retval('Eigen::Vector6d'), [], is_const=True)

  mv.add_binary_numeric_operator('+')
  mv.add_binary_numeric_operator('-')
  mv.add_binary_numeric_operator('*', left_cppclass=Parameter.new('double', 'scalar'))

def add_force_vec(mod):
  fv = mod.add_class('ForceVec')
  fv.add_constructor([])
  fv.add_constructor([param('Eigen::Vector6d', 'vec')])
  fv.add_constructor([param('Eigen::Vector3d', 'couple'),
                      param('Eigen::Vector3d', 'force')])

  fv.add_method('couple', retval('Eigen::Vector3d'), [], is_const=True)
  fv.add_method('force', retval('Eigen::Vector3d'), [], is_const=True)
  fv.add_method('vector', retval('Eigen::Vector6d'), [], is_const=True)

  fv.add_binary_numeric_operator('+')
  fv.add_binary_numeric_operator('-')
  fv.add_binary_numeric_operator('*', left_cppclass=Parameter.new('double', 'scalar'))

def add_rb_inertia(mod):
  rb = mod.add_class('RBInertia')
  rb.add_constructor([])
  rb.add_constructor([param('double', 'mass'),
                      param('Eigen::Vector3d', 'momentum'),
                      param('Eigen::Matrix3d', 'inertia_matrix')])

  rb.add_method('mass', retval('double'), [])
  rb.add_method('momentum', retval('Eigen::Vector3d'), [])
  rb.add_method('inertia', retval('Eigen::Matrix3d'), [])
  rb.add_method('matrix', retval('Eigen::Matrix6d'), [])

  rb.add_binary_numeric_operator('+')
  rb.add_binary_numeric_operator('*', left_cppclass=Parameter.new('double', 'scalar'))

def add_ab_inertia(mod):
  ab = mod.add_class('ABInertia')
  ab.add_constructor([])
  ab.add_constructor([param('Eigen::Matrix3d', 'mass_matrix'),
                      param('Eigen::Matrix3d', 'generalized_inertia_matrix'),
                      param('Eigen::Matrix3d', 'inertia_matrix')])

  ab.add_method('massMatrix', retval('Eigen::Matrix3d'), [])
  ab.add_method('gInertia', retval('Eigen::Matrix3d'), [])
  ab.add_method('inertia', retval('Eigen::Matrix3d'), [])
  ab.add_method('matrix', retval('Eigen::Matrix6d'), [])

  ab.add_binary_numeric_operator('+')
  ab.add_binary_numeric_operator('*', left_cppclass=Parameter.new('double', 'scalar'))

if __name__ == '__main__':
  if len(sys.argv) < 2:
    sys.exit(1)

  sva = Module('_spacevecalg', cpp_namespace='::sva')
  sva.add_include('<SpaceVecAlg>')

  # import Eigen3 types
  import_eigen3_types(sva)

  add_motion_vec(sva)
  add_force_vec(sva)
  add_rb_inertia(sva)
  add_ab_inertia(sva)

  with open(sys.argv[1], 'w') as f:
    sva.generate(f)

