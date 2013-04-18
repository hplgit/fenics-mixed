#include <ext.h>
#include <iostream>

int main() {
  ODEFieldSolver solver = ODEFieldSolver();
  int n = 3;
  double dt = 0.1;
  double T = 2;
  double* s = new double[n];
  solver.redim(n);               // allocate data
  solver.set_dt(dt);             // set time step
  solver.set_IC(0.1);            // set initial conditions
  double t = 0;
    while (t <= T) {
      solver.set_u(0.2);         // set current environment
      solver.advance_one_timestep();
      solver.get_s(n, s);
      for (int i=0; i<n; i++) {
	std::cout << "s(" << t << ", " << i << ")=" << s[i] << std::endl;
      }
      t += dt; 
    }
  return 0;
}
