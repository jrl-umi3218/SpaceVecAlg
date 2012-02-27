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

from pybindgen.typehandlers.base import ReturnValue, Parameter, ReverseWrapperBase, ForwardWrapperBase

# Eigen type support based on python object type name.
# It's ugly but is the only way i have found to easly support
# Eigen3ToPython type without embedded them in this binding.
def createEigenType(eType):
  class Vector6dParam(Parameter):
    DIRECTIONS = [Parameter.DIRECTION_IN]
    CTYPES = ['Eigen::%s&' % eType, 'Eigen::%s' % eType]

    # TODO implement c to python
    def convert_c_to_python(self, wrapper):
      assert isinstance(wrapper, ReverseWrapperBase)
      wrapper.build_params.add_parameter('O', [self.value])

    def convert_python_to_c(self, wrapper):
      assert isinstance(wrapper, ForwardWrapperBase)
      name = wrapper.declarations.declare_variable('PyBindGenLambdaWithType', self.name,
                                                   '{0, "_eigen3.%s"}' % eType)
      wrapper.parse_params.add_parameter('O&', ['pyBindGenLambdaconverter', '&'+name], self.name)
      wrapper.call_params.append('*(Eigen::%s*)%s.lambda->obj' % (eType, name))

def define_eigen_types(mod):
  v3 = mod.add_class('Vector3d', foreign_cpp_namespace='Eigen')
  v6= mod.add_class('Vector6d', foreign_cpp_namespace='Eigen')

  m3 = mod.add_class('Matrix3d', foreign_cpp_namespace='Eigen')
  m6 = mod.add_class('Matrix6d', foreign_cpp_namespace='Eigen')

  types = [v3, v6, m3, m6]

  def make_tp_name(type):
    type.slots['tp_name'] = '_eigen3.%s' % type.name

  for t in types:
    make_tp_name(t)

def add_motion_vec(mod):
  mv = mod.add_class('MotionVec')
  mv.add_constructor([])
  mv.add_constructor([param('Eigen::Vector6d', 'vec')])

if __name__ == '__main__':
  if len(sys.argv) < 2:
    sys.exit(1)

  sva = Module('_spacevecalg', cpp_namespace='::sva')
  sva.add_include('<SpaceVecAlg>')
  sva.add_include('"PyBindGenLambda.h"')

  # create _eigen3 bridge
  createEigenType('Vector3d')
  createEigenType('Vector6d')

  createEigenType('Matrix3d')
  createEigenType('Matrix6d')

  add_motion_vec(sva)


  with open(sys.argv[1], 'w') as f:
    sva.generate(f)

