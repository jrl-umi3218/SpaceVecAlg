// This file is part of SpaceVecAlg.
// 
// SpaceVecAlg is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// SpaceVecAlg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with SpaceVecAlg.  If not, see <http://www.gnu.org/licenses/>.

#pragma once
#include <string>

/// Generic structure to use pyBindGen Type (v0.15).
typedef struct {
  PyObject_HEAD
  void* obj;
  PyBindGenWrapperFlags flags:8;
} PyBindGenLambda;

/// Structure with type name.
typedef struct {
  PyBindGenLambda* lambda;
  const char* type;
} PyBindGenLambdaWithType;

/// Check if object have the same type as param->type.
int pyBindGenLambdaconverter(PyObject* object, PyBindGenLambdaWithType* param)
{
  if(strcmp(object->ob_type->tp_name, param->type) != 0)
  {
    std::string err = std::string("argument must be ") + param->type + std::string(", not ")
      + object->ob_type->tp_name;
    PyErr_SetString(PyExc_TypeError, err.c_str());
    return 0;
  }
  else
  {
    param->lambda = (PyBindGenLambda*)object;
    return 1;
  }
}

