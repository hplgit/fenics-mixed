from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("ODEFieldSolver", 
                             ["ODEFieldSolver.pyx", "../ODEFieldSolver.cpp"],
                             include_dirs=['..'],
                             language='c++')]
)
