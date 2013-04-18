import dolfin
import numpy

def F(S, U):
  if isinstance(S, dolfin.Function) and isinstance(U, dolfin.Function):
    from dolfin import *
    f = sin(S)*exp(U)
    return f
  if isinstance(S, numpy.ndarray) and isinstance(U, numpy.ndarray):
    from numpy import *
    f = sin(S)*exp(U)
    return f

class ParabolicSolver:
  def __init__(self, N, dt):
    """Set up PDE problem for NxN mesh and time step dt."""
    from dolfin import UnitSquareMesh, FunctionSpace, TrialFunction, \
         TestFunction, Function, dx, dot, grad

    mesh = UnitSquareMesh(N,N)
    self.V = V = FunctionSpace(mesh, "Lagrange", 1)

    u = TrialFunction(V)
    v = TestFunction(V)

    a = u*v*dx + dt*dot(grad(u), grad(v))*dx

    self.a = a
    self.dt = dt
    self.mesh = mesh
    self.U = Function(V)

  def advance_one_timestep(self, f, u_1):
    """
    Solve the PDE for one time step.
    f: the source term in the PDE.
    u_1: solution at the previous time step.
    """
    from dolfin import TestFunction, dx, solve

    V, a, dt = self.V, self.a, self.dt  # strip off self prefix
    v = TestFunction(V)
    L = (u_1 + dt*f)*v*dx

    solve(self.a == L, self.U)
    return self.U

import dolfin
import numpy

N = 12     # mesh partition
dt = 0.01  # time step
parabolicsolver = ParabolicSolver(N, dt)
U1 = dolfin.Function(parabolicsolver.V)
U0 = dolfin.Function(parabolicsolver.V)
U0.vector()[:] = numpy.random.random(parabolicsolver.V.dim())

Q = dolfin.FunctionSpace(parabolicsolver.mesh, "DG", 0)
S0_ex = dolfin.Expression("x[0]")
S0 = dolfin.interpolate(S0_ex, Q)
S1 = dolfin.Function(Q)

import ODEFieldSolver  # import module wrapping the ODE solver
odesolver = ODEFieldSolver.ODEFieldSolver()
odesolver.redim(S0.vector().size(0))
odesolver.set_IC(S0.vector().array())
plot = True

for i in range(0, 23):  # time loop
  f = F(S0, U0)
  U1 = parabolicsolver.advance_one_timestep(f, U0)

  U1c = dolfin.project(U1, Q)

  odesolver.set_u(U1c.vector().array())
  odesolver.advance_one_timestep()
  S1.vector()[:] = odesolver.get_s()

  U0 = U1
  S0 = S1

  if plot:
    dolfin.plot(U1, title="U")
    dolfin.plot(S1, title="S")
    dolfin.interactive()


