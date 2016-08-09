template <int dim>
void AcoustoElastic<dim>::assemble_system()
{
    deallog << "Assembling Acousto-Elastic system... ";
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

