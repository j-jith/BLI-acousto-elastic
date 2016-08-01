// #include "AcoustoElastic.h"

template <int dim>
AcoustoElastic<dim>::AcoustoElastic(ParameterHandler &param, const unsigned int degree, bool a_and_d) : 
               degree(degree),
               assemble_and_dump(a_and_d),
               prm(param),
               triangulation(Triangulation<dim>::maximum_smoothing),
               elastic_fe(FE_Q<dim>(degree), dim*2, FE_Nothing<dim>(), 2),
               acoustic_fe(FE_Nothing<dim>(), dim*2, FE_Q<dim>(degree), 2),
               dof_handler(triangulation)
{
   fe_collection.push_back(elastic_fe);
   fe_collection.push_back(acoustic_fe);

   point_load = false; // point load not applied by default
}


template <int dim>
bool AcoustoElastic<dim>::cell_is_in_solid_domain(const typename hp::DoFHandler<dim>::cell_iterator &cell)
{
   return (cell->material_id() == solid_domain_id);
}

template <int dim>
bool AcoustoElastic<dim>::cell_is_in_fluid_domain(const typename hp::DoFHandler<dim>::cell_iterator &cell)
{
   return (cell->material_id() == fluid_domain_id);
}

template <int dim>
void AcoustoElastic<dim>::make_grid()
{

   //READ MESH FILE, SET cell->set_material_id($$_domain_id)
   GridIn<dim> grid_in;
   grid_in.attach_triangulation(triangulation);
   deallog<<"Opening mesh file."<<std::endl; 
   std::ifstream input_file("grid.msh");
   if(input_file.is_open())
   {
      grid_in.read_msh(input_file);
      deallog<<"Mesh file read successfully."<<std::endl;   
   }
   else
   {
      deallog<<"Unable to open mesh file."<<std::endl;
   }

   //triangulation.refine_global(1);

   GridOut grid_out;
   std::ofstream outmesh("grid_out.msh");
   grid_out.write_msh(triangulation,outmesh);

/*
   for(typename Triangulation<dim>::active_cell_iterator cell=dof_handler.begin_active(); cell!=dof_handler.end(); ++cell)
   {
      if((std::fabs(cell->center()[dim-1]) > 0.006) 
         && (std::fabs(cell->center()[dim-1]) < 0.014)
         && (std::sqrt(pow(cell->center()[0],2)+pow(cell->center()[1],2)) < 0.175)
         )
         cell->set_material_id(solid_domain_id);
      else
         cell->set_material_id(fluid_domain_id);
   }
*/
}

template <int dim>
void AcoustoElastic<dim>::set_active_fe_indices()
{
   for(typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell!=dof_handler.end(); ++cell)
   {
      if(cell_is_in_solid_domain(cell))
         cell->set_active_fe_index(0);
      else if(cell_is_in_fluid_domain(cell))
         cell->set_active_fe_index(1);
      else
         Assert(false,ExcNotImplemented());
   }
}


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


template <int dim>
void AcoustoElastic<dim>::setup_system()
{
   set_active_fe_indices();
   
   dof_handler.distribute_dofs(fe_collection);

   //set constraints here
   constraints.clear();

   // Setting fixed BC at boundary 
   const FEValuesExtractors::Vector displacements(0);
   VectorTools::interpolate_boundary_values(dof_handler, 
                  SOLID_FIXED_BOUNDARY, 
                  ZeroFunction<dim>(dim*2+2), 
                  constraints, 
                  fe_collection.component_mask(displacements));

   constraints.close();
   
   deallog << " No. of active cells: "
      << triangulation.n_active_cells()
      << std::endl
      << " Total no. of cells: "
      << triangulation.n_cells()
      << std::endl;

   my_sparsity_pattern();
      
   system_matrix.reinit(sparsity_pattern);
   system_matrix_aux.reinit(sparsity_pattern);
   mass_matrix.reinit(sparsity_pattern);
   damp_matrix.reinit(sparsity_pattern);
   stif_matrix.reinit(sparsity_pattern);
   system_rhs.reinit(dof_handler.n_dofs());
   solution.reinit(dof_handler.n_dofs());
   global_ktr.reinit(dof_handler.n_dofs());
}

template <int dim>
void AcoustoElastic<dim>::set_point_load(Point<dim> location, Point<dim> direction, double magnitude)
{
    point_load = true;

    point_load_loc = location;
    point_load_dir = direction;
    point_load_mag = magnitude;

}

template <int dim>
void AcoustoElastic<dim>::assemble_system(bool step1)
{
   deallog << "Assembling system... ";
   Timer timer;
   timer.start ();

   mass_matrix=0;
   damp_matrix=0;
   stif_matrix=0;

   QGauss<dim> elastic_quadrature(degree+1);
   QGauss<dim> acoustic_quadrature(degree+1);

   hp::QCollection<dim> q_collection;
   q_collection.push_back(elastic_quadrature);
   q_collection.push_back(acoustic_quadrature);

   hp::FEValues<dim> hp_fe_values(fe_collection, q_collection, update_values | update_gradients | update_quadrature_points | update_JxW_values);

   QGauss<dim-1> common_face_quadrature(degree+1);

   FEFaceValues<dim> elastic_fe_face_values(elastic_fe, common_face_quadrature, update_values);
   FEFaceValues<dim> acoustic_fe_face_values(acoustic_fe, common_face_quadrature, update_values | update_normal_vectors | update_JxW_values);

   const unsigned int elastic_dofs_per_cell = elastic_fe.dofs_per_cell;
   const unsigned int acoustic_dofs_per_cell = acoustic_fe.dofs_per_cell;

   FullMatrix<double> cell_mass_matrix;
   FullMatrix<double> cell_damp_matrix;
   FullMatrix<double> cell_stif_matrix;
   FullMatrix<double> cell_interface_matrix_EA(elastic_dofs_per_cell,acoustic_dofs_per_cell);
   FullMatrix<double> cell_interface_matrix_AE(acoustic_dofs_per_cell,elastic_dofs_per_cell);
   Vector<double> cell_rhs;

   std::vector<types::global_dof_index> local_dof_indices;
   std::vector<types::global_dof_index> neighbor_dof_indices(acoustic_dofs_per_cell);

   // Global Point load vector
   //ElasticRHS<dim> elastic_rhs;
   // Point<dim> f_p(0.120,0.0,0.011);
   // Point<dim> f_n(0,0,1);
   // my_point_source_vector(dof_handler, f_p, f_n, system_rhs);
   // system_rhs*=1e5;
   if(point_load)
   {
        my_point_source_vector(dof_handler, point_load_loc, point_load_dir, system_rhs);
        system_rhs*=point_load_mag;
   }
   
   const FEValuesExtractors::Vector displacementsR(0);
   const FEValuesExtractors::Vector displacementsI(dim);
   Vector<double> phiI_i_v(dim);
   Vector<double> phiR_i_v(dim);

   const FEValuesExtractors::Scalar pressureR(2*dim);
   const FEValuesExtractors::Scalar pressureI(2*dim+1);

   typename hp::DoFHandler<dim>::active_cell_iterator
      cell=dof_handler.begin_active(),
      endc=dof_handler.end();
   
   for(; cell!=endc; ++cell)
   {
      hp_fe_values.reinit(cell);

      const FEValues<dim>& fe_values = hp_fe_values.get_present_fe_values();
      const unsigned int n_q_points=fe_values.n_quadrature_points;
      
      cell_mass_matrix.reinit(cell->get_fe().dofs_per_cell, cell->get_fe().dofs_per_cell);
      cell_damp_matrix.reinit(cell->get_fe().dofs_per_cell, cell->get_fe().dofs_per_cell);
      cell_stif_matrix.reinit(cell->get_fe().dofs_per_cell, cell->get_fe().dofs_per_cell);
      cell_rhs.reinit(cell->get_fe().dofs_per_cell);
     
//------------------------------ELASTIC ASSEMBLY----------------------------------------------------- 
      if(cell_is_in_solid_domain(cell))
      {
         const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
         Assert(dofs_per_cell == elastic_dofs_per_cell, ExcInternalError());

         //std::vector<Vector<double> > elastic_rhs_values(n_q_points, Vector<double>(dim));
          //elastic_rhs.vector_value_list(fe_values.get_quadrature_points(),elastic_rhs_values);
         
         for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
         {
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
               const SymmetricTensor<2,dim> phiR_i_symmgrad = fe_values[displacementsR].symmetric_gradient(i,q_point);
               const SymmetricTensor<2,dim> phiI_i_symmgrad = fe_values[displacementsI].symmetric_gradient(i,q_point);
               const double phiR_i_div = fe_values[displacementsR].divergence(i,q_point);
               const double phiI_i_div = fe_values[displacementsI].divergence(i,q_point);
               
               const Tensor<1,dim> phiR_i = fe_values[displacementsR].value(i,q_point);
               const Tensor<1,dim> phiI_i = fe_values[displacementsI].value(i,q_point);
               
               for (unsigned int j=0; j<dofs_per_cell; ++j)
               {
                  const SymmetricTensor<2,dim> phiR_j_symmgrad = fe_values[displacementsR].symmetric_gradient(j,q_point);
                  const SymmetricTensor<2,dim> phiI_j_symmgrad = fe_values[displacementsI].symmetric_gradient(j,q_point);
                  const double phiR_j_div = fe_values[displacementsR].divergence(j,q_point);
                  const double phiI_j_div = fe_values[displacementsI].divergence(j,q_point);

                  const Tensor<1,dim> phiR_j = fe_values[displacementsR].value(j,q_point);
                  const Tensor<1,dim> phiI_j = fe_values[displacementsI].value(j,q_point);
                  
                  cell_mass_matrix(i,j)+=(-rho_s*phiR_i*phiR_j
                                     -rho_s*phiI_i*phiI_j)*fe_values.JxW(q_point);

                  cell_stif_matrix(i,j)+= ( phiR_i_div*phiR_j_div*lambda_s
                                          + 2*(phiR_i_symmgrad*phiR_j_symmgrad)*mu_s)*fe_values.JxW(q_point)
                                          + ( phiI_i_div*phiI_j_div*lambda_s
                                          + 2*(phiI_i_symmgrad*phiI_j_symmgrad)*mu_s)*fe_values.JxW(q_point);
                  
               }
                  
            }
         }
         // Asssemble RHS
        // for(unsigned int i=0; i<dofs_per_cell; ++i)
        //  {
        //     for(unsigned int q_point=0; q_point<n_q_points; ++q_point)
        //        {
        //           const Tensor<1,dim> phiR_i = fe_values[displacementsR].value(i,q_point);
        //           phiR_i.unroll(phiR_i_v);

        //           //cell_rhs(i) += (phiR_i_v * elastic_rhs_values[q_point]) * rho_s *fe_values.JxW(q_point);
        //           cell_rhs(i) += (phiR_i_v * elastic_rhs_values[q_point]) * rho_s *fe_values.JxW(q_point);

        //           //const Tensor<1,dim> phiI_i = fe_values[displacementsI].value(i,q_point);
        //           //phiI_i.unroll(phiI_i_v);
        //        
        //           //cell_rhs(i) += (phiI_i_v * elastic_rhs_values[q_point]) * rho_values[q_point] *fe_values.JxW(q_point);
        //        }
        //  }
//--------------------------------------------------------------------------------------------------------------
      }
      else
      {
//-----------------------------------------ACOUSTIC ASSEMBLY-----------------------------------------------------
         const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
         Assert(dofs_per_cell == acoustic_dofs_per_cell, ExcInternalError());
         
         for(unsigned int i=0; i<dofs_per_cell; ++i)
         {
            for(unsigned int j=0; j<dofs_per_cell; ++j)
            {
               for(unsigned int q_point=0; q_point<n_q_points; ++q_point)
               {     
                  cell_mass_matrix(i,j) += (1/c_0*1/c_0*
                             fe_values[pressureR].value(i,q_point) *
                             fe_values[pressureR].value(j,q_point) *
                             fe_values.JxW(q_point)) 
                             + (1/c_0*1/c_0*
                             fe_values[pressureI].value(i,q_point) *
                             fe_values[pressureI].value(j,q_point) *
                             fe_values.JxW(q_point));
                  
                  cell_stif_matrix(i,j)+= -(fe_values[pressureR].gradient(i,q_point) *
                             fe_values[pressureR].gradient(j,q_point) *
                             fe_values.JxW(q_point))
                            - (fe_values[pressureI].gradient(i,q_point) *
                             fe_values[pressureI].gradient(j,q_point) *
                             fe_values.JxW(q_point));
               }
            }
         }
        
         // if step 1 of BLI sol. procedure, assemble admittance with k_t^2/k_0^2 = 0.5
         if(step1)
         {
            //Boundary integrals
            double bnd_A0=0;
            for(unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
            {
              if((cell->face(face)->at_boundary()) //&& (cell->face(face)->boundary_indicator()==0))
                || (cell->at_boundary(face)==false && cell_is_in_solid_domain(cell->neighbor(face))))
              {
                  bnd_A0=(1/sqrt(2))*(0.5*sqrt(mu/rho_0/c_0) + (gamma-1)*sqrt(kappa/rho_0/c_0/cp));

                  acoustic_fe_face_values.reinit(cell,face);

                  for(unsigned int i=0; i<acoustic_fe_face_values.dofs_per_cell; ++i)
                      for(unsigned int j=0; j<acoustic_fe_face_values.dofs_per_cell; ++j)
                          if(acoustic_fe.has_support_on_face(i,face) && acoustic_fe.has_support_on_face(j,face))
                              for(unsigned int q_point=0; q_point<acoustic_fe_face_values.n_quadrature_points; ++q_point)
                              {
                                  cell_mass_matrix(i,j) +=  acoustic_fe_face_values[pressureR].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureR].value(j,q_point) *
                                                            (fit_c2/c_0/sqrt(c_0))*bnd_A0*
                                                            acoustic_fe_face_values.JxW(q_point)
                                                            +
                                                            acoustic_fe_face_values[pressureR].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureI].value(j,q_point) *
                                                            (fit_c2/c_0/sqrt(c_0))*bnd_A0*
                                                            acoustic_fe_face_values.JxW(q_point)
                                                            +
                                                            acoustic_fe_face_values[pressureI].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureI].value(j,q_point) *
                                                            (fit_c2/c_0/sqrt(c_0))*bnd_A0*
                                                            acoustic_fe_face_values.JxW(q_point)
                                                            -
                                                            acoustic_fe_face_values[pressureI].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureR].value(j,q_point) *
                                                            (fit_c2/c_0/sqrt(c_0))*bnd_A0*
                                                            acoustic_fe_face_values.JxW(q_point);

                                  cell_damp_matrix(i,j) +=  acoustic_fe_face_values[pressureR].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureR].value(j,q_point) *
                                                            (fit_c1/c_0/sqrt(c_0))*bnd_A0*
                                                            acoustic_fe_face_values.JxW(q_point)
                                                            +
                                                            acoustic_fe_face_values[pressureR].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureI].value(j,q_point) *
                                                            (fit_c1/c_0/sqrt(c_0))*bnd_A0*
                                                            acoustic_fe_face_values.JxW(q_point)
                                                            +
                                                            acoustic_fe_face_values[pressureI].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureI].value(j,q_point) *
                                                            (fit_c1/c_0/sqrt(c_0))*bnd_A0*
                                                            acoustic_fe_face_values.JxW(q_point)
                                                            -
                                                            acoustic_fe_face_values[pressureI].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureR].value(j,q_point) *
                                                            (fit_c1/c_0/sqrt(c_0))*bnd_A0*
                                                            acoustic_fe_face_values.JxW(q_point);

                                  cell_stif_matrix(i,j) +=  acoustic_fe_face_values[pressureR].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureR].value(j,q_point) *
                                                            (fit_c0/c_0/sqrt(c_0))*bnd_A0*
                                                            acoustic_fe_face_values.JxW(q_point)
                                                            +
                                                            acoustic_fe_face_values[pressureR].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureI].value(j,q_point) *
                                                            (fit_c0/c_0/sqrt(c_0))*bnd_A0*
                                                            acoustic_fe_face_values.JxW(q_point)
                                                            +
                                                            acoustic_fe_face_values[pressureI].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureI].value(j,q_point) *
                                                            (fit_c0/c_0/sqrt(c_0))*bnd_A0*
                                                            acoustic_fe_face_values.JxW(q_point)
                                                            -
                                                            acoustic_fe_face_values[pressureI].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureR].value(j,q_point) *
                                                            (fit_c0/c_0/sqrt(c_0))*bnd_A0*
                                                            acoustic_fe_face_values.JxW(q_point);
                                              
                              }
              }
            }
        }



//--------------------------------------------------------------------------------------------------------------

      }

      local_dof_indices.resize(cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
      
      //cell_mass_matrix*=omega*omega;
      constraints.distribute_local_to_global(cell_mass_matrix,cell_rhs,
                     local_dof_indices,
                     mass_matrix,system_rhs);

      //cell_damp_matrix*=omega;
      constraints.distribute_local_to_global(cell_damp_matrix,
                     local_dof_indices,
                     local_dof_indices,
                     damp_matrix);

      constraints.distribute_local_to_global(cell_stif_matrix,
                     local_dof_indices,
                     local_dof_indices,
                     stif_matrix);


//-----------------------------------------INTERFACE ASSEMBLY---------------------------------------------------

      if(cell_is_in_solid_domain(cell))
         for(unsigned int f=0; f < GeometryInfo<dim>::faces_per_cell; ++f)
            if(cell->at_boundary(f)==false && cell_is_in_fluid_domain(cell->neighbor(f)))
            {
               elastic_fe_face_values.reinit(cell,f);
               acoustic_fe_face_values.reinit(cell->neighbor(f),
                           cell->neighbor_of_neighbor(f));

               Assert(elastic_fe_face_values.n_quadrature_points==
                  acoustic_fe_face_values.n_quadrature_points,
                  ExcInternalError());

               const unsigned int n_face_quadrature_points=
                        elastic_fe_face_values.n_quadrature_points;

               cell_interface_matrix_EA=0;
               cell_interface_matrix_AE=0;
               for(unsigned int q_point=0; q_point<n_face_quadrature_points; ++q_point)
               {
                  const Tensor<1,dim> normal_vector = acoustic_fe_face_values.normal_vector(q_point);

                  for(unsigned int i=0; i<elastic_fe_face_values.dofs_per_cell; ++i)
                  {
                     const Tensor<1,dim> d_phiR_i = elastic_fe_face_values[displacementsR].value(i,q_point);
                     const Tensor<1,dim> d_phiI_i = elastic_fe_face_values[displacementsI].value(i,q_point);

                     for(unsigned int j=0; j<acoustic_fe_face_values.dofs_per_cell; ++j)
                     {
                        const double p_phiR_j = acoustic_fe_face_values[pressureR].value(j,q_point);
                        const double p_phiI_j = acoustic_fe_face_values[pressureI].value(j,q_point);

                        cell_interface_matrix_EA(i,j) += -d_phiR_i*(p_phiR_j*normal_vector)*acoustic_fe_face_values.JxW(q_point);
                        cell_interface_matrix_EA(i,j) += -d_phiI_i*(p_phiI_j*normal_vector)*acoustic_fe_face_values.JxW(q_point);

                        cell_interface_matrix_AE(j,i) += rho_0*(d_phiR_i*normal_vector)*p_phiR_j*acoustic_fe_face_values.JxW(q_point);
                        cell_interface_matrix_AE(j,i) += rho_0*(d_phiI_i*normal_vector)*p_phiI_j*acoustic_fe_face_values.JxW(q_point);

                     }
                  }
               }
                //cell_interface_matrix.print(deallog); deallog<<std::endl;
               cell->neighbor(f)->get_dof_indices(neighbor_dof_indices);
               constraints.distribute_local_to_global(cell_interface_matrix_EA,
                                                      local_dof_indices,
                                                      neighbor_dof_indices,
                                                      stif_matrix);
               //cell_interface_matrix_AE*=omega*omega;
               constraints.distribute_local_to_global(cell_interface_matrix_AE,
                                                      neighbor_dof_indices,
                                                      local_dof_indices,
                                                      mass_matrix);
//---------------------------------------------------------------------------------------------------------------
            }


   }

   timer.stop ();
   deallog << "done (" << timer.wall_time() << "s)" << std::endl;
}


template <int dim>
void AcoustoElastic<dim>::assemble_system_aux()
{
   deallog << "Assembling system(aux)... ";
   Timer timer;
   timer.start ();

   system_matrix_aux=0;
   global_ktr=0;

   QGauss<dim> elastic_quadrature(degree+1);
   QGauss<dim> acoustic_quadrature(degree+1);

   hp::QCollection<dim> q_collection;
   q_collection.push_back(elastic_quadrature);
   q_collection.push_back(acoustic_quadrature);

   hp::FEValues<dim> hp_fe_values(fe_collection, q_collection, update_values | update_gradients | 
                                    update_hessians | update_quadrature_points | update_JxW_values);

   QGauss<dim-1> common_face_quadrature(degree+1);

   FEFaceValues<dim> elastic_fe_face_values(elastic_fe, common_face_quadrature, update_values);

   FEFaceValues<dim> acoustic_fe_face_values(acoustic_fe, common_face_quadrature, update_values | 
                                                update_normal_vectors | update_hessians | 
                                                update_quadrature_points | update_JxW_values);

   //const unsigned int elastic_dofs_per_cell = elastic_fe.dofs_per_cell;
   const unsigned int acoustic_dofs_per_cell = acoustic_fe.dofs_per_cell;

   FullMatrix<double> cell_matrix;
   Vector<double> cell_ktr;

   std::vector<types::global_dof_index> local_dof_indices;
   std::vector<types::global_dof_index> neighbor_dof_indices(acoustic_dofs_per_cell);

   //ElasticRHS<dim> elastic_rhs;
   
   const FEValuesExtractors::Vector displacementsR(0);
   const FEValuesExtractors::Vector displacementsI(dim);

   const FEValuesExtractors::Scalar pressureR(2*dim);
   const FEValuesExtractors::Scalar pressureI(2*dim+1);

   typename hp::DoFHandler<dim>::active_cell_iterator
      cell=dof_handler.begin_active(),
      endc=dof_handler.end();
   
   for(; cell!=endc; ++cell)
   {

      if(cell_is_in_fluid_domain(cell))
      {
         hp_fe_values.reinit(cell);

         const FEValues<dim>& fe_values = hp_fe_values.get_present_fe_values();
         const unsigned int n_q_points=fe_values.n_quadrature_points;
      
         cell_matrix.reinit(cell->get_fe().dofs_per_cell, cell->get_fe().dofs_per_cell);
         cell_ktr.reinit(cell->get_fe().dofs_per_cell);
     
         //Boundary integrals

         std::vector<double> local_sol_lap_R;   
         std::vector<double> local_sol_lap_I;   
         std::vector<double> local_sol_R;   
         std::vector<double> local_sol_I;   
         std::vector<double> k_t_ratio;   
         double k_tmp_R, k_tmp_I;
         std::vector<double> bnd_A0;

        for(unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; ++face)
        {
            if((cell->face(face)->at_boundary()) //&& (cell->face(face)->boundary_indicator()==0))
            || (cell->at_boundary(face)==false && cell_is_in_solid_domain(cell->neighbor(face))))
            {
                //bnd_A0=(1/sqrt(2))*(sqrt(mu/rho_0/c_0) + (gamma-1)*sqrt(kappa/rho_0/c_0/cp));

                acoustic_fe_face_values.reinit(cell,face);

               local_sol_lap_R.resize(n_q_points);
               local_sol_lap_I.resize(n_q_points);
               local_sol_R.resize(n_q_points);
               local_sol_I.resize(n_q_points);
               k_t_ratio.resize(n_q_points);
               bnd_A0.resize(n_q_points);

               fe_values[pressureR].get_function_laplacians(solution,local_sol_lap_R);
               fe_values[pressureI].get_function_laplacians(solution,local_sol_lap_I);
               fe_values[pressureR].get_function_values(solution,local_sol_R);
               fe_values[pressureI].get_function_values(solution,local_sol_I);
               
               for(unsigned int q_point=0; q_point<acoustic_fe_face_values.n_quadrature_points; ++q_point)
               {
                  k_tmp_R = -(local_sol_lap_R[q_point]*local_sol_R[q_point] 
                            + local_sol_lap_I[q_point]*local_sol_I[q_point])
                            /(pow(local_sol_R[q_point],2)+pow(local_sol_I[q_point],2))/pow(omega/c_0,2); 

                  k_tmp_I = -(local_sol_lap_I[q_point]*local_sol_R[q_point] 
                            - local_sol_lap_R[q_point]*local_sol_I[q_point])
                            /(pow(local_sol_R[q_point],2)+pow(local_sol_I[q_point],2))/pow(omega/c_0,2);

                  k_t_ratio[q_point]=std::min(std::sqrt(pow(k_tmp_R,2)+pow(k_tmp_I,2)),1.0);
                  //k_t_ratio[q_point]=0.5;
                  bnd_A0[q_point]=(1/sqrt(2))*(k_t_ratio[q_point]*sqrt(mu/rho_0/c_0) + (gamma-1)*sqrt(kappa/rho_0/c_0/cp));
               }

                for(unsigned int i=0; i<acoustic_fe_face_values.dofs_per_cell; ++i)
               {
                  for(unsigned int q=0; q<acoustic_fe_face_values.n_quadrature_points; ++q)
                           cell_ktr(i)+=k_t_ratio[q];      
         
                    for(unsigned int j=0; j<acoustic_fe_face_values.dofs_per_cell; ++j)
                        if(acoustic_fe.has_support_on_face(i,face) && acoustic_fe.has_support_on_face(j,face))
                            for(unsigned int q_point=0; q_point<acoustic_fe_face_values.n_quadrature_points; ++q_point)
                            {

                                cell_matrix(i,j) += omega*omega*(acoustic_fe_face_values[pressureR].value(i,q_point) *
                                                          acoustic_fe_face_values[pressureR].value(j,q_point) *
                                                          (fit_c2/c_0/sqrt(c_0))*bnd_A0[q_point]* 
                                                            acoustic_fe_face_values.JxW(q_point)
                                                        +
                                                          acoustic_fe_face_values[pressureR].value(i,q_point) *
                                                          acoustic_fe_face_values[pressureI].value(j,q_point) *
                                                          (fit_c2/c_0/sqrt(c_0))*bnd_A0[q_point]* 
                                                            acoustic_fe_face_values.JxW(q_point)
                                                        +
                                                          acoustic_fe_face_values[pressureI].value(i,q_point) *
                                                          acoustic_fe_face_values[pressureI].value(j,q_point) *
                                                          (fit_c2/c_0/sqrt(c_0))*bnd_A0[q_point]* 
                                                          acoustic_fe_face_values.JxW(q_point)
                                                        -
                                                          acoustic_fe_face_values[pressureI].value(i,q_point) *
                                                          acoustic_fe_face_values[pressureR].value(j,q_point) *
                                                          (fit_c2/c_0/sqrt(c_0))*bnd_A0[q_point]* 
                                                          acoustic_fe_face_values.JxW(q_point));

                                cell_matrix(i,j) +=  omega*(acoustic_fe_face_values[pressureR].value(i,q_point) *
                                                          acoustic_fe_face_values[pressureR].value(j,q_point) *
                                                          (fit_c1/c_0/sqrt(c_0))*bnd_A0[q_point]* 
                                                          acoustic_fe_face_values.JxW(q_point)
                                                        +
                                                          acoustic_fe_face_values[pressureR].value(i,q_point) *
                                                          acoustic_fe_face_values[pressureI].value(j,q_point) *
                                                          (fit_c1/c_0/sqrt(c_0))*bnd_A0[q_point]* 
                                                          acoustic_fe_face_values.JxW(q_point)
                                                        +
                                                          acoustic_fe_face_values[pressureI].value(i,q_point) *
                                                          acoustic_fe_face_values[pressureI].value(j,q_point) *
                                                          (fit_c1/c_0/sqrt(c_0))*bnd_A0[q_point]* 
                                                          acoustic_fe_face_values.JxW(q_point)
                                                        -
                                                          acoustic_fe_face_values[pressureI].value(i,q_point) *
                                                          acoustic_fe_face_values[pressureR].value(j,q_point) *
                                                          (fit_c1/c_0/sqrt(c_0))*bnd_A0[q_point]* 
                                                          acoustic_fe_face_values.JxW(q_point));

                                cell_matrix(i,j) +=  acoustic_fe_face_values[pressureR].value(i,q_point) *
                                                          acoustic_fe_face_values[pressureR].value(j,q_point) *
                                                          (fit_c0/c_0/sqrt(c_0))*bnd_A0[q_point]* 
                                                          acoustic_fe_face_values.JxW(q_point)
                                                        +
                                                          acoustic_fe_face_values[pressureR].value(i,q_point) *
                                                          acoustic_fe_face_values[pressureI].value(j,q_point) *
                                                          (fit_c0/c_0/sqrt(c_0))*bnd_A0[q_point]* 
                                                          acoustic_fe_face_values.JxW(q_point)
                                                        +
                                                          acoustic_fe_face_values[pressureI].value(i,q_point) *
                                                          acoustic_fe_face_values[pressureI].value(j,q_point) *
                                                          (fit_c0/c_0/sqrt(c_0))*bnd_A0[q_point]* 
                                                          acoustic_fe_face_values.JxW(q_point)
                                                        -
                                                          acoustic_fe_face_values[pressureI].value(i,q_point) *
                                                          acoustic_fe_face_values[pressureR].value(j,q_point) *
                                                          (fit_c0/c_0/sqrt(c_0))*bnd_A0[q_point]* 
                                                          acoustic_fe_face_values.JxW(q_point);

                                          
                            }
               }
            }
        }



//--------------------------------------------------------------------------------------------------------------
      local_dof_indices.resize(cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
      
      //cell_mass_matrix*=omega*omega;
      constraints.distribute_local_to_global(cell_matrix,
                     local_dof_indices,
                     local_dof_indices,
                     system_matrix_aux);

      constraints.distribute_local_to_global(cell_ktr,local_dof_indices,global_ktr);

      }

   }

   timer.stop ();
   deallog << "done (" << timer.wall_time() << "s)" << std::endl;
}


template <int dim>
void AcoustoElastic<dim>::dump_matrices()
{
    deallog << "Dumping matrices... ";
    Timer timer;
    timer.start ();

    std::ofstream mass_out("mass_out.txt");
    std::ofstream damp_out("damp_out.txt");
    std::ofstream stif_out("stif_out.txt");
    std::ofstream rhs_out("rhs_out.txt");
    mass_matrix.print(mass_out);
    damp_matrix.print(damp_out);
    stif_matrix.print(stif_out);
    system_rhs.print(rhs_out);

   timer.stop ();
   deallog << "done (" << timer.wall_time() << "s)" << std::endl;
}


template <int dim>
void AcoustoElastic<dim>::solve()
{
   system_matrix=0;
   system_matrix.add(omega*omega,mass_matrix);
   system_matrix.add(omega,damp_matrix);
   system_matrix.add(1.0,stif_matrix);
   system_matrix.add(1.0,system_matrix_aux);
    
    deallog << "Solving linear system (ndof="<< dof_handler.n_dofs()<<")...";
    Timer timer;
    timer.start ();
    
   SparseDirectUMFPACK A_direct;
   A_direct.initialize(system_matrix);
   A_direct.vmult (solution, system_rhs);

    //SolverControl solver_control(1e4,1e-10);

    //SolverCG<> solver(solver_control);
    //SolverGMRES<> solver(solver_control);

    //solver.solve(system_matrix,solution,system_rhs,PreconditionIdentity());
    
    timer.stop ();
    deallog << "done (" << timer.wall_time() << "s)" << std::endl;
}

template <int dim>
void AcoustoElastic<dim>::read_input(unsigned int n_sol)
{
   //const unsigned int N_dofs=dof_handler.n_dofs();
   const std::string input_sol_name = "input/solution_";
    const std::string input_ext = ".dat";
    std::string input_num;
    std::string input_file;

   input_num = patch::to_string(n_sol);
   input_file = input_sol_name+input_num+input_ext;
   std::ifstream infile(input_file.c_str());

   deallog <<"Reading " <<input_file.c_str()<<std::endl;

   solution.block_read(infile);

   infile.close();
}

template <int dim>
void AcoustoElastic<dim>::read_pod()
{
   const unsigned int N_dofs=dof_handler.n_dofs();
   const std::string input_mod_name = "pod_input/solution_";
    const std::string input_ext = ".dat";
    std::string input_num;
    std::string input_file;

   std::vector<Vector<double> > input_sols_tmp(N_pod, Vector<double> (N_dofs));

   for(unsigned int i=0; i<N_pod; i++)
   {
      input_num = patch::to_string(i);
      input_file = input_mod_name+input_num+input_ext;
      std::ifstream infile(input_file.c_str());

      deallog <<"Reading " <<input_file.c_str()<<std::endl;

      input_sols_tmp[i].block_read(infile);

      infile.close();

      //Subtract mean
      //input_sols_tmp[i].add(-input_sols_tmp[i].mean_value());

      input_sols_tmp[i]*=std::sqrt((f_1-f_0)/(n_f-1));
   }
  
   deallog<<"Populating correlation matrix...";
   Timer timer;
   timer.restart();
   double lambda_max=0;
   //double lambda_max_tmp=0;
   FullMatrix<double> cor_matrix_full(N_pod,N_pod);
   for(unsigned int i=0; i<N_pod; i++)
   {
      //lambda_max_tmp=0;
      for(unsigned int j=i; j<N_pod; j++)
      {
         cor_matrix_full(i,j)=(input_sols_tmp[i]*input_sols_tmp[j])/N_pod;
         cor_matrix_full(j,i)=cor_matrix_full(i,j);
         //lambda_max_tmp+=abs(cor_matrix(i,j));
      }
      //if(lambda_max_tmp>lambda_max)
        // lambda_max=lambda_max_tmp;

      //lambda_max+=cor_matrix_full(i,i);
   }
   //cor_matrix_full.symmetrize();
   std::ofstream cormat_outfile("cormat.dat");
   cor_matrix_full.print(cormat_outfile,12,5);
   cormat_outfile.close();

   lambda_max=cor_matrix_full.trace();
   LAPACKFullMatrix<double> cor_matrix(N_pod,N_pod);
   cor_matrix=cor_matrix_full;

   timer.stop();
   deallog<<"done("<<timer.wall_time()<<" s)"<<std::endl;
   deallog<<"lambda_max = "<<lambda_max<<std::endl;

   deallog<<"Computing eigevalues"<<std::endl;
   timer.restart();
   Vector<double> cor_evals(N_pod);
   FullMatrix<double> cor_evecs(N_pod, N_pod);
   cor_matrix.compute_eigenvalues_symmetric(0, lambda_max, 0, cor_evals, cor_evecs);
   deallog<<cor_evals.size()<<" eigenvalues found"<<std::endl;

   //cor_evals.print();
   std::ofstream eigen_file("eigen.dat");
   cor_evals.print(eigen_file, 6, true, false);
   eigen_file.close();

   unsigned int n_pod=0;
   unsigned int n_pod_0=0;
   double full_sum=cor_evals.mean_value()*cor_evals.size();
   double partial_sum=0;
   //double pod_thresh=2.5e-15;

   while(n_pod_0<cor_evals.size() && partial_sum<pod_thresh)
   {
      partial_sum+=cor_evals(n_pod_0)/full_sum;
      n_pod_0+=1;
   }
   deallog<<"Eigenvalues to be discarded (eps="<<pod_thresh<<") = "<<n_pod_0<<std::endl;

   //while(partial_sum<pod_thresh && n_pod<cor_evals.size())
   while(n_pod<(cor_evals.size()-n_pod_0))
   //while((cor_evals(cor_evals.size()-1-n_pod)/cor_evals(cor_evals.size()-1)>1e-16) && n_pod<cor_evals.size())
   {
      //partial_sum += cor_evals(cor_evals.size()-1-n_pod)/full_sum;
      deallog<<cor_evals(cor_evals.size()-1-n_pod)<<std::endl;
      for(unsigned int i=0; i<N_pod; i++)
         cor_evecs(i,cor_evals.size()-1-n_pod)/=std::sqrt(cor_evals(cor_evals.size()-1-n_pod));
      n_pod+=1;
   }
   timer.stop();
   deallog<<"done("<<timer.wall_time()<<" s)"<<std::endl;
   deallog<<"No. of POD modes = "<<n_pod<<std::endl;
  
   deallog<<"Pushing and writing POD modes...";
   timer.restart();
   pod_modes.resize(n_pod);
   std::string pod_name;
   for(unsigned int j=0; j<n_pod; j++)
   {
      pod_modes[j].reinit(dof_handler.n_dofs());
      for(unsigned int i=0; i<N_pod; i++)
      {
         pod_modes[j].add(cor_evecs(i,cor_evals.size()-1-j),input_sols_tmp[i]);
      }
      pod_name="binary_om/pod_mode_"+patch::to_string(j)+".dat";
      std::ofstream pod_outfile(pod_name.c_str());
      pod_modes[j].block_write(pod_outfile);
      pod_outfile.close();
   }
   timer.stop();
   deallog<<"done("<<timer.wall_time()<<" s)"<<std::endl;


   N_pod=n_pod;
   red_mass_matrix.reinit(N_pod,N_pod);
   red_stif_matrix.reinit(N_pod,N_pod);
   red_damp_matrix.reinit(N_pod,N_pod);
   red_rhs.reinit(N_pod);
   Vector<double> tmp_vec_m(N_dofs);
   Vector<double> tmp_vec_c(N_dofs);
   Vector<double> tmp_vec_k(N_dofs);

   deallog<<"Populating reduced matrix...";
   Timer timer1;
   timer1.start();
   red_mass_matrix=0;
   red_stif_matrix=0;
   red_damp_matrix=0;
   red_rhs=0;
   for(unsigned int j=0; j<N_pod; j++)
   {
      mass_matrix.vmult(tmp_vec_m,pod_modes[j]);
      damp_matrix.vmult(tmp_vec_c,pod_modes[j]);
      stif_matrix.vmult(tmp_vec_k,pod_modes[j]);
      for(unsigned int i=0; i<N_pod; i++)
      {
         //red_matrix(i,j) = system_matrix.matrix_scalar_product(arnoldi_vecs[i],arnoldi_vecs[j]); 
         red_mass_matrix(i,j)=tmp_vec_m*pod_modes[i];
         red_damp_matrix(i,j)=tmp_vec_c*pod_modes[i];
         red_stif_matrix(i,j)=tmp_vec_k*pod_modes[i];
      }
      red_rhs(j) = system_rhs*pod_modes[j];
      //std::cout<<"\n"<<i<<std::endl;

   }
   timer1.stop ();
   deallog << "done ("  << timer1.wall_time() << "s)" << std::endl;
}

template <int dim>
void AcoustoElastic<dim>::pod_aux()
{
   const unsigned int N_dofs=dof_handler.n_dofs();
   red_matrix_aux.reinit(N_pod,N_pod);
   Vector<double> tmp_vec(N_dofs);

   deallog<<"Populating reduced matrix(aux)...";
   Timer timer1;
   timer1.start();
   red_matrix_aux=0;
   for(unsigned int j=0; j<N_pod; j++)
   {
      system_matrix_aux.vmult(tmp_vec,pod_modes[j]);
      for(unsigned int i=0; i<N_pod; i++)
      {
         red_matrix_aux(i,j)=tmp_vec*pod_modes[i];
      }
   }
   timer1.stop ();
   deallog << "done ("  << timer1.wall_time() << "s)" << std::endl;

}

template <int dim>
void AcoustoElastic<dim>::solve_reduced()
{


   const unsigned int N_dofs=dof_handler.n_dofs();

   unsigned int N_size=0;

   N_size=N_pod;

   FullMatrix<double> red_matrix(N_size,N_size);
   Vector<double> red_solution(N_size);
   //Vector<double> tmp_vec(N_dofs);

   red_matrix=0;
   red_matrix.add(omega*omega,red_mass_matrix);
   red_matrix.add(omega,red_damp_matrix);
   red_matrix.add(1.0,red_stif_matrix);
   red_matrix.add(1.0,red_matrix_aux);


   deallog<<"Solving reduced system...";
   red_solution=0;
   red_matrix.gauss_jordan();
   red_matrix.vmult(red_solution, red_rhs);
   deallog << "done" << std::endl;

   solution=0;
   for(unsigned int i=0; i<N_size; i++)
   {
      solution.add(red_solution(i),pod_modes[i]);
   } 
}


template <int dim>
void AcoustoElastic<dim>::run(unsigned int step)
{
   deallog << "Initialising and reading grid... ";

   make_grid();


   f_0=prm.get_double("f_0");
   f_1=prm.get_double("f_1");

   if(f_1>f_0)
   {
      n_f=prm.get_double("n_f");
      sweep=true;
   }
   else
   {
      omega=f_0*2*M_PI;
      sweep=false;
   }
   
   rho_s=prm.get_double("rho_s");
   
   nu_s=prm.get_double("nu_s");
   
   E_s=prm.get_double("E_s");
   
   mu_s=E_s/2/(1+nu_s);
   lambda_s=E_s*nu_s/(1+nu_s)/(1-2*nu_s);

   rho_0=prm.get_double("rho_0");
   c_0=prm.get_double("c_0");
   mu=prm.get_double("mu");
   kappa=prm.get_double("kappa");
   cp=prm.get_double("cp");
   gamma=prm.get_double("gamma");

   fit_c2=prm.get_double("fit_c2");
   fit_c1=prm.get_double("fit_c1");
   fit_c0=prm.get_double("fit_c0");

   N_pod=prm.get_integer("N_pod");
   pod_thresh=prm.get_double("pod_thresh");

   //Output parameters for troubleshooting
   std::ofstream prm_out("prm_out.txt");
   prm.print_parameters(prm_out, ParameterHandler::Text);

   setup_system();

   if(step == BLI_STEP_1) // BLI acousto-elastic step-1
        assemble_system(true);
   else // BLI acousto-elastic step-2 OR plain acousto-elastic (no admittance)
       assemble_system(false);

   //read_pod();

   if(assemble_and_dump)
   {
      dump_matrices();
   }
   else
   {
      if(sweep)
      {
         for(unsigned int n_sweep=0; n_sweep<n_f; n_sweep++)
         {
            omega=(f_0+n_sweep*(f_1-f_0)/(n_f-1))*2*M_PI;
            deallog <<"Frequency: " <<omega/2/M_PI<<" Hz"<<std::endl;

            if(step == BLI_STEP_2) // BLI acousto-elastic step-2
            {
                read_input(n_sweep);
                assemble_system_aux();
                //output_ktr(n_sweep);
            }

            solve();
            //output_results(n_sweep);
            //pod_aux();
            //solve_reduced();
            //output_binary(n_sweep);
            output_results(n_sweep);
         }
      }
      else
      {
         deallog <<"Frequency: " <<omega/2/M_PI<<" Hz"<<std::endl;
         solve();
         output_results(0); 
      }
   }
}

template <int dim>
void AcoustoElastic<dim>::output_binary(unsigned int num)
{
   const std::string output_name = "binary_op/solution_";
    const std::string output_num = patch::to_string(num);
    //const std::string output_ext = ".gpl";
    const std::string output_ext = ".dat";
    const std::string output_file = output_name+output_num+output_ext;

    std::ofstream outfile(output_file.c_str());

   solution.block_write(outfile);
}

template <int dim>
void AcoustoElastic<dim>::output_ktr(unsigned int num)
{
   const std::string output_name = "ktr_op/ktr_";
    const std::string output_num = patch::to_string(num);
    //const std::string output_ext = ".gpl";
    const std::string output_ext = ".dat";
    const std::string output_file = output_name+output_num+output_ext;

    std::ofstream outfile(output_file.c_str());

   global_ktr.block_write(outfile);
}

template <int dim>
void AcoustoElastic<dim>::output_results (unsigned int num)
{
    deallog << "Generating output... ";
    Timer timer;
    timer.start ();
    
    DataOut<dim,hp::DoFHandler<dim> > data_out;
    data_out.attach_dof_handler (dof_handler);
    std::vector<std::string> solution_names;
    switch (dim)
      {
      case 1:
        solution_names.push_back ("displacement_real");
        solution_names.push_back ("displacement_imag");
        solution_names.push_back ("pressure_real");
        solution_names.push_back ("pressure_imag");
        break;
      case 2:
        solution_names.push_back ("x_displacement_real");
        solution_names.push_back ("y_displacement_real");
        solution_names.push_back ("x_displacement_imag");
        solution_names.push_back ("y_displacement_imag");
        solution_names.push_back ("pressure_real");
        solution_names.push_back ("pressure_imag");
        break;
      case 3:
        solution_names.push_back ("x_displacement_real");
        solution_names.push_back ("y_displacement_real");
        solution_names.push_back ("z_displacement_real");
        solution_names.push_back ("x_displacement_imag");
        solution_names.push_back ("y_displacement_imag");
        solution_names.push_back ("z_displacement_imag");
        solution_names.push_back ("pressure_real");
        solution_names.push_back ("pressure_imag");
        break;
      default:
        Assert (false, ExcNotImplemented());
      }
    data_out.add_data_vector (solution, solution_names);
    data_out.build_patches ();

   const std::string output_name = "sweep/solution_";
    const std::string output_num = patch::to_string(num);
    //const std::string output_ext = ".gpl";
    const std::string output_ext = ".vtu";
    const std::string output_file = output_name+output_num+output_ext;

    std::ofstream outfile(output_file.c_str());
    //data_out.write_gnuplot(outfile);
    data_out.write_vtu(outfile);
    deallog<<"Output written to "<<output_file<<std::endl;

    timer.stop ();
    deallog << "done ("
            << timer()
            << "s)"
            << std::endl;
}
