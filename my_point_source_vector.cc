#include "my_point_source_vector.h"

template <int dim, int spacedim>
void my_point_source_vector (const hp::MappingCollection<dim,spacedim> &mapping,
                             const hp::DoFHandler<dim,spacedim>        &dof_handler,
                             const Point<spacedim>                     &p,
                             const Point<dim>                          &orientation,
                             Vector<double>                            &rhs_vector)
{
    Assert (rhs_vector.size() == dof_handler.n_dofs(),
            ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
    //Assert (dof_handler.get_fe().n_components() == dim,
    //        ExcMessage ("This function only works for vector-valued finite elements."));

    rhs_vector = 0;

    std::pair<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
    cell_point = GridTools::find_active_cell_around_point (mapping, dof_handler, p);

    Quadrature<dim> q(GeometryInfo<dim>::project_to_unit_cell(cell_point.second));

    const FEValuesExtractors::Vector vec (0);
    FEValues<dim> fe_values(mapping[cell_point.first->active_fe_index()],
                            cell_point.first->get_fe(), q, UpdateFlags(update_values));
    fe_values.reinit(cell_point.first);

    const unsigned int dofs_per_cell = cell_point.first->get_fe().dofs_per_cell;

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    cell_point.first->get_dof_indices (local_dof_indices);

    for (unsigned int i=0; i<dofs_per_cell; i++)
        rhs_vector(local_dof_indices[i]) =  orientation * fe_values[vec].value(i,0);
}



template <int dim, int spacedim>
void my_point_source_vector (const hp::DoFHandler<dim,spacedim>   &dof_handler,
                             const Point<spacedim>                &p,
                             const Point<dim>                     &orientation,
                             Vector<double>                       &rhs_vector)
{
    my_point_source_vector(hp::StaticMappingQ1<dim>::mapping_collection,
                           dof_handler,
                           p, orientation, rhs_vector);
}

