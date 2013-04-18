%module ODEFieldSolver
%{
#include <arrayobject.h>
#include <sstream>
#include "ODEFieldSolver.h"
%}


%init %{
import_array();
%}

%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (int n, double* array) {
$1 = PyArray_Check($input) ? 1 : 0;
}

%typemap(in) (int n, double* array){
  if (!PyArray_Check($input)) {
    PyErr_SetString(PyExc_TypeError, "Not a NumPy array");
    return NULL; ;
  }
  PyArrayObject* pyarray;
  pyarray = (PyArrayObject*)$input;
  if (!(PyArray_TYPE(pyarray) == NPY_DOUBLE)) {
    PyErr_SetString(PyExc_TypeError, "Not a NumPy array of doubles");
    return NULL; ;
  }
  $1 = int(pyarray->dimensions[0]);
  $2 = (double*)pyarray->data;
}

/* Wrap ODEFieldSolver::get_s in a Python function */
%rename (_get_s) ODEFieldSolver::get_s;

%extend ODEFieldSolver{
 %pythoncode%{
    def get_s(self):
      import numpy as np
      a = np.zeros(self.size())
      self._get_s(a)
      return a
 %}
}

%include std_string.i
%include "ODEFieldSolver.h"






