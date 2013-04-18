"""Demo of the Python interface to the pure C version of ODEFieldSolver."""
import ODEFieldSolvercpp as solver 
import numpy

s0 = numpy.array([0, 1, 4], float)  # ensure float elements!
u = numpy.array([0, 1, 1], float)
n = s0.size
s = numpy.zeros(n)

solver.set_ic_and_dt(s0, dt=0.1)
for n in range(1, 8):
    solver.set_u(u)
    solver.advance_one_timestep()
    s = solver.get_s(s)
    print n, s

