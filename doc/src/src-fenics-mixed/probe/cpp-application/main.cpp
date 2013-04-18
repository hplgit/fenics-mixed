#include <dolfin.h>
#include "Probe.h"
#include "Lagrange1.h"

using namespace dolfin; 

class U0 : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = x[0];  
  }
}; 

int main(){ 
  UnitCubeMesh mesh = UnitCubeMesh(10, 10, 10);
  Lagrange1::FunctionSpace V(mesh);
  Array<double> x(3); 
  x[0] = 0.5; x[1] = 0.5; x[2] = 0.5;  
  Probe probe(x, V); 
  U0 u0; 
  Function u(V); 
  u.interpolate(u0); 
  probe.eval(u); 
  std::cout <<" number of probes "<< probe.value_size() <<std::endl; 
  std::cout <<" evaluation "<< probe.get_probe(0)[0] <<std::endl; 
}


