cimport numpy as cnp
import numpy as np

cdef extern from "ODEFieldSolver.h":

    cppclass ODEFieldSolver_ "ODEFieldSolver":
        ODEFieldSolver()
        void redim(int n)
        int size()
        double g(double s, double u)
        void set_dt(double dt)
        void set_IC(int n, double* array)
        void set_const_IC(double s0)
        void set_u (int n, double* array)
        void set_IC(double val)
        void set_u (double const_value)
        void advance_one_timestep()
        void get_s  (int n, double* array)

cdef class ODEFieldSolver:
    cdef ODEFieldSolver_  *wrapped

    def __init__(self):
        self.wrapped = new ODEFieldSolver_()

    def __dealloc__(self):
        if self.wrapped != NULL:
            del self.wrapped

    def redim(self, n):
        self.wrapped.redim(n)

    def size(self):
        return self.wrapped.size()

    def g(self, s, u):
        return self.wrapped.g(s, u)

    def set_dt(self, dt):
        self.wrapped.set_dt(dt)

    def set_IC(self, cnp.ndarray[double, ndim=1, mode='c'] array):
        if array.shape[0] != self.wrapped.size():
            raise ValueError('incorrect dimension on array')
        self.wrapped.set_IC(array.shape[0], &array[0])

    # FIXME: what about function overloading?
    def set_u(self, cnp.ndarray[double, ndim=1, mode='c'] array):
        if array.shape[0] != self.wrapped.size():
            raise ValueError('incorrect dimension on array')
        self.wrapped.set_IC(array.shape[0], &array[0])

    def advance_one_timestep(self):
        self.wrapped.advance_one_timestep()


    def set_const_IC(self, s0):
        self.wrapped.set_const_IC(s0)

    def get_s(self,
              cnp.ndarray[double, ndim=1, mode='c'] out=None):
        if out is None:
            out = np.empty(self.wrapped.size(), dtype=np.double)
        elif out.shape[0] != self.wrapped.size():
            raise ValueError('incorrect dimension on out')
        self.wrapped.get_s(out.shape[0], &out[0])
        return out

