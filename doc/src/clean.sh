#!/bin/sh -x
doconce clean
rm -f *.py
rm -rf doc/._part*.rst
rm -f automake*
src=src-fenics-mixed
cd $src/ODEFieldSolver
rm -f app
cd swig
rm -f ODEFieldSolver.py ODEFieldSolver_wrap.cxx
cd ../f2py-c
sh clean.sh
cd ../f2py-cpp
sh clean.sh
