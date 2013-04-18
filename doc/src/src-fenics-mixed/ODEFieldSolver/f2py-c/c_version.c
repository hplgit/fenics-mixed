#include <stdlib.h>
#include <stdio.h>

/* global variables */
double* s;
double* s_1;
double* u;
double dt;
int n;

void redim(int n_)
{
  n = n_;
  s =   malloc(sizeof(double)*n);
  s_1 = malloc(sizeof(double)*n);
  u =   malloc(sizeof(double)*n);
}

void deallocate()
{
  free(s); free(s_1); free(u);
}

/* Note: do not mix upper and lower case letters as in set_IC_...
   This leads to undefined symbols when f2py compiles the code.
*/

void set_ic_and_dt(int n_, double* s0, double dt_)
{
  int i;
  redim(n_);
  dt = dt_;
  for (i=0; i<n; i++) {
    s_1[i] = s0[i];
  }
}

void set_const_ic_and_dt(int n_, double s0, double dt_)
{
  int i;
  redim(n_);
  dt = dt;
  for (i=0; i<n; i++) {
    s_1[i] = s0;
  }
}

void set_u(int n_, double* u_)
{
  int i;
  for (i=0; i<n; i++) {
    u[i] = u_[i];
  }
}

double g(double s_, double u_) {
  /* return s_*u_*(1 - s_); */
  return s_;
}

void advance_one_timestep()
{
  /* Use the Forward Euler time integration for simplicity */
  int i;
  for (i=0; i<n; i++) {
    s[i] = s_1[i] + dt*g(s_1[i], u[i]);
    /* For debugging: */
    /* printf("i=%d, s_1=%g, dt=%g, g=%g, s=%g\n",
       i, s_1[i], dt, g(s_1[i], u[i]), s[i]); */
  }
  /* Update for next time step */
  for (i=0; i<n; i++) { s_1[i] = s[i]; }

}

void get_s(int n_, double* s_)
{
  int i;
  for (i=0; i<n; i++) {
    s_[i] = s[i];
  }
}
