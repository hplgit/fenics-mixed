import NSSolver
from BCModel import BCModel

solver = NSSolver.NSSolver()
solver.setIC()
t = 0
dt = 0.01   # fluid solver time step
T = 5.0     # end time of simulation

C = 0.127
R = 5.43
N = 1000    # use N steps in the ODE solver in [t,t+dt]

num_outlets = 2  # no out outflow boundaries

# Create an ODE model for the pressure at each outlet boundary
pmodels = [BCModel(C, R, N, dt) for i in range(0,num_outlets)]

P_ = [16000, 100]   # start values for outlet pressures

while t < T:
    t += dt

    # Compute u_ and p_ using known outlet pressures P_
    solver.advance_one_time_step(P_, t)
    # Compute the flux at outlet boundaries
    Q = solver.flux()

    # Advance outlet pressure boundary condition to the
    # next time step (for each outlet boundary)
    # (pmodels returns a vector of size N containg the 
    # the solution between [t, t+dt]. 
    # We take the last one with [-1])
    for i in range(0, num_outlets):
        P_[i] = pmodels[i](P_[i], Q[i])[-1]

