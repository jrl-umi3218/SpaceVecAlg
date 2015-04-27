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

from _spacevecalg import *
import eigen3 as e3

# for compatibility purpose default SpaceVecAlg type
# are referencing double implementation

MotionVec = MotionVecd
ForceVec = ForceVecd
RBInertia = RBInertiad
ABInertia = ABInertiad
PTransform = PTransformd


# pickle, cPickle and copy support for all sva type

def matrix3dFromList(matList):
  """Convert a 9 double list into a eigen3.Matrix3d"""
  m3d = e3.Matrix3d()
  for i in xrange(len(m3d)):
    m3d[i] = matList[i]
  return m3d


def motionvecConstructor(angular, linear):
  return MotionVecd(e3.Vector3d(*angular), e3.Vector3d(*linear))

def motionvecPickle(mv):
  return motionvecConstructor, (list(mv.angular()), list(mv.linear()))


def forcevecConstructor(couple, force):
  return ForceVecd(e3.Vector3d(*couple), e3.Vector3d(*force))

def forcevecPickle(fv):
  return forcevecConstructor, (list(fv.couple()), list(fv.force()))


def ptransformConstructor(rotation, translation):
  return PTransformd(matrix3dFromList(rotation), e3.Vector3d(*translation))

def ptransformPickle(pt):
  return ptransformConstructor, (list(pt.rotation()), list(pt.translation()))


def rbinertiaConstructor(mass, momentum, inertia):
  return RBInertiad(mass, e3.Vector3d(*momentum), matrix3dFromList(inertia))

def rbinertiaPickle(rb):
  return rbinertiaConstructor, (rb.mass(), list(rb.momentum()),
                                list(rb.inertia()))


def abinertiaConstructor(M, H, I):
  return ABInertiad(matrix3dFromList(M), matrix3dFromList(H),
                    matrix3dFromList(I))

def abinertiaPickle(ab):
  return abinertiaConstructor, (list(ab.massMatrix()), list(ab.gInertia()),
                                list(ab.inertia()))

def copy_reg_pickle():
  # python 2 and 3 support
  # first try to import copyreg (python 3)
  # if the import fail we import copy_reg (python 2)
  try:
    import copyreg
  except ImportError:
    import copy_reg as copyreg

  copyreg.pickle(MotionVecd, motionvecPickle)
  copyreg.pickle(ForceVecd, forcevecPickle)
  copyreg.pickle(PTransformd, ptransformPickle)
  copyreg.pickle(RBInertiad, rbinertiaPickle)
  copyreg.pickle(ABInertiad, abinertiaPickle)
