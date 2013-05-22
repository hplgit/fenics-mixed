
"This demo program demonstrates how to include additional C++ code in DOLFIN."

# Copyright (C) 2013 Kent-Andre Mardal, Mikael Mortensen, Johan Hake
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2013-04-02


from dolfin import *
import numpy
import os

header_file = open("Probe/Probe.h", "r")
code = header_file.read()
header_file.close()
probe_module = compile_extension_module(
    code=code, source_directory="Probe", sources=["Probe.cpp"],
    include_dirs=[".", os.path.abspath("Probe")])

mesh = UnitCubeMesh(10, 10, 10)
V = FunctionSpace(mesh, 'CG', 1)

x = numpy.array((0.5, 0.5, 0.5))
probe = probe_module.Probe(x, V)

u0 = interpolate(Expression('x[0]'), V)
# Fast evaluation of U0 at x:
probe.eval(u0)
print "The number of probes is ", probe.number_of_eval_calls()
print "The value at ", x, " is ", probe.get_values(0)




