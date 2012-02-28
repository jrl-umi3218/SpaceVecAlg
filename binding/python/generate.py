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

  mod.add_class('Quaterniond', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

def build_motion_vec(mv):
  mv.add_copy_constructor()
  mv.add_constructor([])
  mv.add_constructor([param('Eigen::Vector6d', 'vec')])
  mv.add_constructor([param('Eigen::Vector3d', 'angular'),
                      param('Eigen::Vector3d', 'linear')])

  mv.add_method('angular', retval('Eigen::Vector3d'), [], is_const=True)
  mv.add_method('linear', retval('Eigen::Vector3d'), [], is_const=True)
  mv.add_method('vector', retval('Eigen::Vector6d'), [], is_const=True)

  mv.add_function_as_method('cross', retval('sva::MotionVec'),
                            [param('sva::MotionVec', 'm1'), param('sva::MotionVec', 'm2')])
  mv.add_function_as_method('crossDual', retval('sva::ForceVec'),
                            [param('sva::MotionVec', 'm1'), param('sva::ForceVec', 'm2')])
  mv.add_function_as_method('dot', retval('double'),
                            [param('sva::MotionVec', 'm1'), param('sva::ForceVec', 'm2')])

  mv.add_binary_numeric_operator('+')
  mv.add_binary_numeric_operator('-')
  mv.add_binary_numeric_operator('*', left_cppclass=Parameter.new('double', 'scalar'))

def build_force_vec(fv):
  fv.add_copy_constructor()
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

def build_rb_inertia(rb):
  rb.add_copy_constructor()
  rb.add_constructor([])
  rb.add_constructor([param('double', 'mass'),
                      param('Eigen::Vector3d', 'momentum'),
                      param('Eigen::Matrix3d', 'inertia_matrix')])

  rb.add_method('mass', retval('double'), [], is_const=True)
  rb.add_method('momentum', retval('Eigen::Vector3d'), [], is_const=True)
  rb.add_method('inertia', retval('Eigen::Matrix3d'), [], is_const=True)
  rb.add_method('matrix', retval('Eigen::Matrix6d'), [], is_const=True)

  rb.add_binary_numeric_operator('+')
  rb.add_binary_numeric_operator('*', left_cppclass=Parameter.new('double', 'scalar'))
  rb.add_binary_numeric_operator('*', ReturnValue.new('sva::ForceVec'), right=param('sva::MotionVec', 'mv'))

def build_ab_inertia(ab):
  ab.add_constructor([])
  ab.add_copy_constructor()
  ab.add_constructor([param('Eigen::Matrix3d', 'mass_matrix'),
                      param('Eigen::Matrix3d', 'generalized_inertia_matrix'),
                      param('Eigen::Matrix3d', 'inertia_matrix')])

  ab.add_method('massMatrix', retval('Eigen::Matrix3d'), [], is_const=True)
  ab.add_method('gInertia', retval('Eigen::Matrix3d'), [], is_const=True)
  ab.add_method('inertia', retval('Eigen::Matrix3d'), [], is_const=True)
  ab.add_method('matrix', retval('Eigen::Matrix6d'), [], is_const=True)

  ab.add_binary_numeric_operator('+')
  ab.add_binary_numeric_operator('*', left_cppclass=Parameter.new('double', 'scalar'))
  rb.add_binary_numeric_operator('+', right=param('sva::RBInertia', 'rb'))
  rb.add_binary_numeric_operator('*', ReturnValue.new('sva::ForceVec'), right=param('sva::MotionVec', 'mv'))

def build_p_transform(mod):
  pt.add_constructor([])
  pt.add_copy_constructor()
  pt.add_constructor([param('Eigen::Matrix3d', 'rot'),
                      param('Eigen::Vector3d', 'trans')])
  pt.add_constructor([param('Eigen::Quaterniond', 'rot'),
                      param('Eigen::Vector3d', 'trans')])
  pt.add_constructor([param('Eigen::Quaterniond', 'rot')])
  pt.add_constructor([param('Eigen::Matrix3d', 'rot')])
  pt.add_constructor([param('Eigen::Vector3d', 'trans')])

  pt.add_method('rotation', retval('Eigen::Matrix3d'), [], is_const=True)
  pt.add_method('translation', retval('Eigen::Vector3d'), [], is_const=True)
  pt.add_method('matrix', retval('Eigen::Matrix6d'), [], is_const=True)
  pt.add_method('dualMatrix', retval('Eigen::Matrix6d'), [], is_const=True)

  pt.add_method('inv', retval('sva::PTransform'), [])

  pt.add_binary_numeric_operator('*')

  pt.add_binary_numeric_operator('*', ReturnValue.new('sva::MotionVec'), right=param('sva::MotionVec', 'mv'))
  pt.add_function_as_method('invMul', retval('sva::MotionVec'),
                            [param('sva::PTransform', 'p1'), param('sva::MotionVec', 'p2')])

  pt.add_function_as_method('dualMul', retval('sva::ForceVec'),
                            [param('sva::PTransform', 'p1'), param('sva::ForceVec', 'p2')])
  pt.add_function_as_method('transMul', retval('sva::ForceVec'),
                            [param('sva::PTransform', 'p1'), param('sva::ForceVec', 'p2')])

  pt.add_function_as_method('dualMul', retval('sva::RBInertia'),
                            [param('sva::PTransform', 'p1'), param('sva::RBInertia', 'p2')])
  pt.add_function_as_method('transMul', retval('sva::RBInertia'),
                            [param('sva::PTransform', 'p1'), param('sva::RBInertia', 'p2')])

  pt.add_function_as_method('dualMul', retval('sva::ABInertia'),
                            [param('sva::PTransform', 'p1'), param('sva::ABInertia', 'p2')])
  pt.add_function_as_method('transMul', retval('sva::ABInertia'),
                            [param('sva::PTransform', 'p1'), param('sva::ABInertia', 'p2')])


if __name__ == '__main__':
  if len(sys.argv) < 2:
    sys.exit(1)

  sva = Module('_spacevecalg', cpp_namespace='::sva')
  sva.add_include('<SpaceVecAlg>')

  # import Eigen3 types
  import_eigen3_types(sva)

  mv = sva.add_class('MotionVec')
  fv = sva.add_class('ForceVec')
  rb = sva.add_class('RBInertia')
  ab = sva.add_class('ABInertia')
  pt = sva.add_class('PTransform')

  build_motion_vec(mv)
  build_force_vec(fv)
  build_rb_inertia(rb)
  build_ab_inertia(ab)
  build_p_transform(pt)

  with open(sys.argv[1], 'w') as f:
    sva.generate(f)

