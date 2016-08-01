#ifndef PARAMETER_READER_H
#define PARAMETER_READER_H

#include "includes.h"
#include <deal.II/base/parameter_handler.h>

// Parameter class
class ParameterReader : public Subscriptor
{
  public:
    ParameterReader(ParameterHandler &);
    void read_parameters(const std::string);

  private:
    void declare_parameters();
    ParameterHandler &prm;
};

#endif
