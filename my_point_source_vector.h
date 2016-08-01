#ifndef MY_POINT_SOURCE_VECTOR_H
#define MY_POINT_SOURCE_VECTOR_H

#include "includes.h"

template <int dim, int spacedim> 
void my_point_source_vector (const hp::MappingCollection<dim,spacedim> &mapping,
                             const hp::DoFHandler<dim,spacedim>        &dof_handler,
                             const Point<spacedim>                     &p,
                             const Point<dim>                          &orientation,
                             Vector<double>                            &rhs_vector);

template <int dim, int spacedim> 
void my_point_source_vector (const hp::DoFHandler<dim,spacedim>   &dof_handler,
                             const Point<spacedim>                &p,
                             const Point<dim>                     &orientation,
                             Vector<double>                       &rhs_vector);

// template functions implementation
#include "my_point_source_vector.tcc"

#endif
