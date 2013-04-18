#!/bin/sh
f2py -c -m bcmodelf77 ../PMODEL.f
# Test if module is correctly built
python -c 'import bcmodelf77
print bcmodelf77.__doc__
print bcmodelf77.pmodel.__doc__'
cp bcmodelf77.so ../../



