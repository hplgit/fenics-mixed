#include <ODEFieldSolver.h>
#include <iostream>
#include <assert.h>

ODEFieldSolver::ODEFieldSolver()
{ n = 0; s = 0; s_1 =0 ; u = 0; }

void ODEFieldSolver::redim(int n_)
{
  n = n_;
  s = new double[n];
  s_1 = new double[n];
  u = new double[n];
}

int ODEFieldSolver::size() { return n; }

ODEFieldSolver::~ODEFieldSolver()
{
  delete [] s;
  delete [] s_1;
  delete [] u;
}

double ODEFieldSolver:: g(double s, double u)
{
  return s*u*(1 - s);
}

void ODEFieldSolver::advance_one_timestep()
{
  // Use the Forward Euler time integration for simplicity
  for (int i=0; i<n; i++)
  {
    s[i] = s_1[i] + dt*g(s_1[i], u[i]);
  }
  // Update for nODEFieldSolver time step
  for (int i=0; i<n; i++) { s_1[i] = s[i]; }
}

void ODEFieldSolver::set_IC(int n_, double* s0)
{
  assert(n == n_);
  for (int i=0; i<n; i++) { s_1[i] = s0[i]; }
}

void ODEFieldSolver::set_const_IC(double s0)
{
  for (int i=0; i<n; i++) { s_1[i] = s0; }
}

void ODEFieldSolver::set_IC(double s0)
{
  set_const_IC(s0); 
}



void ODEFieldSolver::set_u(int n_, double* u_)
{
  assert(n == n_);
  for (int i=0; i<n; i++) { u[i] = u_[i]; }
}

void ODEFieldSolver::set_u(double u_)
{
  for (int i=0; i<n; i++) { u[i] = u_; }
}

void ODEFieldSolver::set_dt(double dt_)
{
  dt = dt_;
}

void ODEFieldSolver::get_s(int n_, double* s_)
{
  assert(n == n_);
  n_ = n;
  for (int i=0; i<n; i++) { s_[i] = s[i]; }
}
