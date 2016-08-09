template <int dim>
void AcoustoElastic<dim>::my_sparsity_pattern()
{
   deallog << "Building sparsity pattern... ";
   Timer timer;
   timer.start();

   CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());

   QGauss<dim> elastic_quadrature(degree+1);
   QGauss<dim> acoustic_quadrature(degree+1);

   hp::QCollection<dim> q_collection;
   q_collection.push_back(elastic_quadrature);
   q_collection.push_back(acoustic_quadrature);

   //hp::FEValues<dim> hp_fe_values(fe_collection, q_collection, update_values | update_gradients | update_quadrature_points | update_JxW_values);
   hp::FEValues<dim> hp_fe_values(fe_collection, q_collection, update_default);

   QGauss<dim-1> common_face_quadrature(degree+1);

   //FEFaceValues<dim> elastic_fe_face_values(elastic_fe, common_face_quadrature, update_values);
   FEFaceValues<dim> elastic_fe_face_values(elastic_fe, common_face_quadrature, update_default);
   //FEFaceValues<dim> acoustic_fe_face_values(acoustic_fe, common_face_quadrature, update_values | update_normal_vectors | update_JxW_values);
   FEFaceValues<dim> acoustic_fe_face_values(acoustic_fe, common_face_quadrature, update_default);

   //const unsigned int elastic_dofs_per_cell = elastic_fe.dofs_per_cell;
   const unsigned int acoustic_dofs_per_cell = acoustic_fe.dofs_per_cell;

   typename hp::DoFHandler<dim>::active_cell_iterator
      cell=dof_handler.begin_active(),
      endc=dof_handler.end();

   std::vector<types::global_dof_index> local_dof_indices;
   std::vector<types::global_dof_index> neighbor_dof_indices(acoustic_dofs_per_cell);

   for(; cell!=endc; ++cell)
   {
      hp_fe_values.reinit(cell);
      //const FEValues<dim>& fe_values = hp_fe_values.get_present_fe_values();
      //const unsigned int n_q_points=fe_values.n_quadrature_points;

      local_dof_indices.resize(cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);

      if(cell_is_in_solid_domain(cell))
      {
         const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
         Assert(dofs_per_cell == elastic_dofs_per_cell, ExcInternalError());

         for(unsigned int i=0; i<dofs_per_cell; ++i)
            for(unsigned int j=0; j<dofs_per_cell; ++j)
               c_sparsity.add(local_dof_indices[i],local_dof_indices[j]);


         for(unsigned int f=0; f < GeometryInfo<dim>::faces_per_cell; ++f)
            if(cell->at_boundary(f)==false && cell_is_in_fluid_domain(cell->neighbor(f)))
            {
               elastic_fe_face_values.reinit(cell,f);
               acoustic_fe_face_values.reinit(cell->neighbor(f),
                           cell->neighbor_of_neighbor(f));

               Assert(elastic_fe_face_values.n_quadrature_points==
                  acoustic_fe_face_values.n_quadrature_points,
                  ExcInternalError());

               cell->neighbor(f)->get_dof_indices(neighbor_dof_indices);

               for(unsigned int i=0; i<elastic_fe_face_values.dofs_per_cell; ++i)
                  for(unsigned int j=0; j<acoustic_fe_face_values.dofs_per_cell; ++j)
                  {
                     c_sparsity.add(local_dof_indices[i],neighbor_dof_indices[j]);
                     c_sparsity.add(neighbor_dof_indices[j],local_dof_indices[i]);
                  }
            }
      }
      else
      {
         const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
         Assert(dofs_per_cell == acoustic_dofs_per_cell, ExcInternalError());

         for(unsigned int i=0; i<dofs_per_cell; ++i)
            for(unsigned int j=0; j<dofs_per_cell; ++j)
               c_sparsity.add(local_dof_indices[i],local_dof_indices[j]);

         for(unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
         {
            if((cell->face(face)->at_boundary()) //&& (cell->face(face)->boundary_indicator()==0))
            || (cell->at_boundary(face)==false && cell_is_in_solid_domain(cell->neighbor(face))))
            {
               acoustic_fe_face_values.reinit(cell,face);

               for(unsigned int i=0; i<acoustic_fe_face_values.dofs_per_cell; ++i)
                  for(unsigned int j=0; j<acoustic_fe_face_values.dofs_per_cell; ++j)
                     c_sparsity.add(local_dof_indices[i],local_dof_indices[j]);
            }

         }
      }
   }

  // DoFTools::make_flux_sparsity_pattern(dof_handler,c_sparsity);

   constraints.condense(c_sparsity);
   sparsity_pattern.copy_from(c_sparsity);

   //std::ofstream sp_out("sparsity_pattern.gpl");
  // sparsity_pattern.print_gnuplot(sp_out);

   timer.stop ();
   deallog << "done (" << timer.wall_time() << "s)" << std::endl;
}

