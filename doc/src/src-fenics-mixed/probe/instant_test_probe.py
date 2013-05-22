
from dolfin import *
import time
import numpy
import os
import instant

code = """
namespace dolfin
{
  class Probe
  {
  public:
    Probe(const Array<double>& point, const FunctionSpace& V);
    void eval(const Function& u);
    std::vector<double> get_values(std::size_t component);

    std::size_t num_components() {return value_size_loc;};
    std::size_t number_of_eval_calls() {return _probes[0].size();};
    std::vector<double> get_point();
    void erase(std::size_t i);
    void clear();

  private:
    std::vector<std::vector<double> > basis_matrix;
    std::vector<double> coefficients;
    double _x[3];
    boost::shared_ptr<const FiniteElement> _element;
    Cell* dolfin_cell;
    UFCCell* ufc_cell;
    std::size_t value_size_loc;
    std::vector<std::vector<double> > _probes;
  };
}
"""

additional_decl = """
%init%{
import_array();
%}

// Include global SWIG interface files:
// Typemaps, shared_ptr declarations, exceptions, version
%include <boost_shared_ptr.i>

// Global typemaps and forward declarations
%include "dolfin/swig/typemaps/includes.i"
%include "dolfin/swig/forwarddeclarations.i"

// Global exceptions
%include <exception.i>

// Local shared_ptr declarations
%shared_ptr(dolfin::Function)
%shared_ptr(dolfin::FunctionSpace)

// %import types from submodule function of SWIG module function
%import(module="dolfin.cpp.function") "dolfin/function/Function.h"
%import(module="dolfin.cpp.function") "dolfin/function/FunctionSpace.h"

%feature("autodoc", "1");
"""

system_headers = ['numpy/arrayobject.h',
                  'dolfin/function/Function.h',
                  'dolfin/function/FunctionSpace.h']
swigargs = ['-c++', '-fcompact', '-O', '-I.', '-small']
cmake_packages = ['DOLFIN']
sources = ["Probe.cpp"]
source_dir = "Probe"
include_dirs = [".", os.path.abspath("Probe")]

compiled_module = instant.build_module(
    code=code,
    source_directory=source_dir,
    additional_declarations=additional_decl,
    system_headers=system_headers,
    include_dirs=include_dirs,
    swigargs=swigargs,
    sources=sources,
    cmake_packages=cmake_packages)

mesh = UnitCubeMesh(10, 10, 10)
V = FunctionSpace(mesh, 'CG', 1)

x = numpy.array((0.5, 0.5, 0.5))
probe = compiled_module.Probe(x, V)

# Just create some random data to be used for probing
u0 = interpolate(Expression('x[0]'), V)
probe.eval(u0)
print "The number of probes is ", probe.number_of_eval_calls()
print "The value at ", x, " is ", probe.get_values(0)




