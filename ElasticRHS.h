#ifndef ELASTIC_RHS_H
#define ELASTIC_RHS_H

#include "includes.h"

template <int dim>
class ElasticRHS : public Function<dim>
{
   public:
   ElasticRHS();
   virtual void vector_value(const Point<dim> &p, Vector<double> &values) const;
   virtual void vector_value_list(const std::vector<Point<dim> > &points, std::vector<Vector<double> > &value_list) const;
   
};

// template class implementation
#include "ElasticRHS.tcc"

#endif
