#ifndef __PROBE_H
#define __PROBE_H

#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>

namespace dolfin
{
  class Probe
  {
  public:
    Probe(const Array<double>& point, const FunctionSpace& V);
    void eval(const Function& u);
    std::vector<double> get_values(std::size_t component);

    std::size_t num_components() {return value_size_loc;};
    std::size_t number_of_eval_calls() {return _probes[0].size();};
    std::vector<double> get_point();
    void erase(std::size_t i);
    void clear();

  private:
    std::vector<std::vector<double> > basis_matrix;
    std::vector<double> coefficients;
    double _x[3];
    boost::shared_ptr<const FiniteElement> _element;
    Cell* dolfin_cell;
    UFCCell* ufc_cell;
    std::size_t value_size_loc;
    std::vector<std::vector<double> > _probes;
  };
}
#endif
