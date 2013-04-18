#!/bin/sh -x
sh ./clean.sh
name=ODEFieldSolvercpp

f2py -m $name -h $name.pyf --overwrite-signature signatures_capi2cpp.f

f2py -c --fcompiler=gfortran -I.. --build-dir tmp1 \
     -DF2PY_REPORT_ON_ARRAY_COPY=1 $name.pyf ../ODEFieldSolver.cpp capi2cpp.cpp

python -c "import $name; print dir($name); print $name.__doc__"
