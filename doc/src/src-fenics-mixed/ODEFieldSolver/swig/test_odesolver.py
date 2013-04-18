import numpy
import ODEFieldSolver
from numpy import random

N = 12 
dt = 0.01

S = numpy.zeros(N)
odesolver = ODEFieldSolver.ODEFieldSolver()
odesolver.redim(N)
odesolver.set_dt(dt)
s_ic = random.random(N)
odesolver.set_IC(s_ic)

t = 0
T = 4.5  
dt = 0.1 
while t<=T: 
  odesolver.set_u(numpy.ones(N))
  odesolver.advance_one_timestep();
  t += dt 
  

s = odesolver.get_s()

import pylab
pylab.plot(s_ic)
pylab.plot(s)
pylab.show()



