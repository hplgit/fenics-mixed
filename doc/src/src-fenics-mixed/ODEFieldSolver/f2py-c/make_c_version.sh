#!/bin/sh
sh ./clean.sh
name=ODEFieldSolverc

f2py -m $name -h $name.pyf --overwrite-signature signatures_c_version.f

f2py -c --fcompiler=gfortran -I.. --build-dir tmp1 \
     -DF2PY_REPORT_ON_ARRAY_COPY=1 $name.pyf c_version.c

python -c "import $name; print dir($name); print $name.__doc__"
