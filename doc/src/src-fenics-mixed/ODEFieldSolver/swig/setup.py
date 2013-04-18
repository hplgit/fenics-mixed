import os, numpy
from distutils.core import setup, Extension
name = 'ODEFieldSolver'
swig_cmd ='swig -python   -c++  -I.. -I.  %s.i' % name
os.system(swig_cmd)
sources = ['../%s.cpp' % name, '%s_wrap.cxx' % name]
setup(name=name,
      ext_modules= [
          Extension('_' + name,
                    sources,
                    include_dirs=['..',
                                  numpy.get_include() + "/numpy"])])



