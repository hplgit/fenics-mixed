

class ODEFieldSolver {
  int n;                // no of points (or regions)
  double* s;            // discrete values of unknown s)
  double* s_1;          // s at the previous time level
  double* u;            // discrete values of external field u
  double dt;            // time step size

  public:
  ODEFieldSolver();
 ~ODEFieldSolver();
  void redim(int n);    // allocate data structures
  int size();           // return the no of points/regions
  virtual double g(double s, double u);
  void set_dt(double dt);
  void set_IC(int n, double* array);
  void set_const_IC(double s0);
  void set_u (int n, double* array);
  void set_IC(double const_value);
  void set_u (double const_value);
  void advance_one_timestep();

  void get_s  (int n, double* array);
};



