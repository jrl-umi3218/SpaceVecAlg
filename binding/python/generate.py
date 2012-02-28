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

if __name__ == '__main__':
  if len(sys.argv) < 2:
    sys.exit(1)

  sva = Module('_spacevecalg', cpp_namespace='::sva')
  sva.add_include('<SpaceVecAlg>')

  # import Eigen3 types
  import_eigen3_types(sva)

  add_motion_vec(sva)

  with open(sys.argv[1], 'w') as f:
    sva.generate(f)

