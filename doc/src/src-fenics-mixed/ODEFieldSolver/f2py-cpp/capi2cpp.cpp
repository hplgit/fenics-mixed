#include "ODEFieldSolver.h"

ODEFieldSolver solver = ODEFieldSolver();

extern "C" {

void set_ic_and_dt(int n, double* s0, double dt)
{
  solver.redim(n);
  solver.set_dt(dt);
  solver.set_IC(n, s0);
}

void set_const_ic_and_dt(int n, double s0, double dt)
{
  solver.redim(n);
  solver.set_dt(dt);
  solver.set_const_IC(s0);
}

void set_u(int n, double* u)
{
  solver.set_u(n, u);
}

void advance_one_timestep()
{
  solver.advance_one_timestep();
}

void get_s(int n, double* s)
{
  solver.get_s(n, s);
}

}
