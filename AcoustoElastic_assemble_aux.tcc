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

