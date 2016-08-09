template <int dim>
void AcoustoElastic<dim>::assemble_impedance()
{
    switch(solution_step)
    {
        case BLI_STEP_1:
            assemble_impedance_step1_unsplit();
            break;
        case BLI_STEP_1_APPROX:
            assemble_impedance_step1_split();
            break;
        case BLI_STEP_2:
            assemble_impedance_step2_unsplit();
            break;
        case BLI_STEP_2_APPROX:
            assemble_impedance_step2_split();
            break;
        default:
            deallog << "Undefined solution step. Impedance assembly not done." << std::endl;
    }
}


template <int dim>
void AcoustoElastic<dim>::assemble_impedance_step1_unsplit()
{
    deallog << "Assembling impedance(Step 1 - Unsplit)... ";
    Timer timer;
    timer.start ();

    system_matrix_aux=0;

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


    FullMatrix<double> cell_matrix;

    std::vector<types::global_dof_index> local_dof_indices;
    std::vector<types::global_dof_index> neighbor_dof_indices(acoustic_fe.dofs_per_cell);

    //ElasticRHS<dim> elastic_rhs;

    const FEValuesExtractors::Vector displacementsR(0);
    const FEValuesExtractors::Vector displacementsI(dim);

    const FEValuesExtractors::Scalar pressureR(2*dim);
    const FEValuesExtractors::Scalar pressureI(2*dim+1);

    typename hp::DoFHandler<dim>::active_cell_iterator
        cell=dof_handler.begin_active(),
        endc=dof_handler.end();

    double bnd_A0=(1/sqrt(2))*(0.5*sqrt(mu/rho_0/c_0) + (gamma-1)*sqrt(kappa/rho_0/c_0/cp));

    for(; cell!=endc; ++cell)
    {

        if(cell_is_in_fluid_domain(cell))
        {
            hp_fe_values.reinit(cell);

            const FEValues<dim>& fe_values = hp_fe_values.get_present_fe_values();
            const unsigned int n_q_points=fe_values.n_quadrature_points;

            cell_matrix.reinit(cell->get_fe().dofs_per_cell, cell->get_fe().dofs_per_cell);

            //Boundary integrals
            for(unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; ++face)
            {
                if((cell->face(face)->at_boundary()) //&& (cell->face(face)->boundary_indicator()==0))
                    || (cell->at_boundary(face)==false && cell_is_in_solid_domain(cell->neighbor(face))))
                    {

                        acoustic_fe_face_values.reinit(cell,face);

                        for(unsigned int i=0; i<acoustic_fe_face_values.dofs_per_cell; ++i)
                        {
                            for(unsigned int j=0; j<acoustic_fe_face_values.dofs_per_cell; ++j)
                                if(acoustic_fe.has_support_on_face(i,face) && acoustic_fe.has_support_on_face(j,face))
                                    for(unsigned int q_point=0; q_point<acoustic_fe_face_values.n_quadrature_points; ++q_point)
                                    {

                                        cell_matrix(i,j) += acoustic_fe_face_values[pressureR].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureR].value(j,q_point) *
                                                            (omega/c_0*sqrt(omega/c_0))*bnd_A0* 
                                                            acoustic_fe_face_values.JxW(q_point)
                                                            +
                                                            acoustic_fe_face_values[pressureR].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureI].value(j,q_point) *
                                                            (omega/c_0/sqrt(omega/c_0))*bnd_A0* 
                                                            acoustic_fe_face_values.JxW(q_point)
                                                            +
                                                            acoustic_fe_face_values[pressureI].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureI].value(j,q_point) *
                                                            (omega/c_0/sqrt(omega/c_0))*bnd_A0* 
                                                            acoustic_fe_face_values.JxW(q_point)
                                                            -
                                                            acoustic_fe_face_values[pressureI].value(i,q_point) *
                                                            acoustic_fe_face_values[pressureR].value(j,q_point) *
                                                            (omega/c_0/sqrt(omega/c_0))*bnd_A0* 
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

        }

    }

    timer.stop ();
    deallog << "done (" << timer.wall_time() << "s)" << std::endl;
}

template <int dim>
void AcoustoElastic<dim>::assemble_impedance_step1_split()
{
    deallog << "Assembling impedance(Step 1 - Split)... ";
    Timer timer;
    timer.start ();

    system_matrix_aux=0;

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


    FullMatrix<double> cell_mass_matrix;
    FullMatrix<double> cell_damp_matrix;
    FullMatrix<double> cell_stif_matrix;
    Vector<double> cell_rhs;

    std::vector<types::global_dof_index> local_dof_indices;
    std::vector<types::global_dof_index> neighbor_dof_indices(acoustic_fe.dofs_per_cell);

    //ElasticRHS<dim> elastic_rhs;

    const FEValuesExtractors::Vector displacementsR(0);
    const FEValuesExtractors::Vector displacementsI(dim);

    const FEValuesExtractors::Scalar pressureR(2*dim);
    const FEValuesExtractors::Scalar pressureI(2*dim+1);

    typename hp::DoFHandler<dim>::active_cell_iterator
        cell=dof_handler.begin_active(),
        endc=dof_handler.end();

    double bnd_A0=(1/sqrt(2))*(0.5*sqrt(mu/rho_0/c_0) + (gamma-1)*sqrt(kappa/rho_0/c_0/cp));

    for(; cell!=endc; ++cell)
    {

        if(cell_is_in_fluid_domain(cell))
        {
            hp_fe_values.reinit(cell);

            const FEValues<dim>& fe_values = hp_fe_values.get_present_fe_values();
            const unsigned int n_q_points=fe_values.n_quadrature_points;

            cell_mass_matrix.reinit(cell->get_fe().dofs_per_cell, cell->get_fe().dofs_per_cell);
            cell_damp_matrix.reinit(cell->get_fe().dofs_per_cell, cell->get_fe().dofs_per_cell);
            cell_stif_matrix.reinit(cell->get_fe().dofs_per_cell, cell->get_fe().dofs_per_cell);
            cell_rhs.reinit(cell->get_fe().dofs_per_cell);

            //Boundary integrals
            for(unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; ++face)
            {
                if((cell->face(face)->at_boundary()) //&& (cell->face(face)->boundary_indicator()==0))
                    || (cell->at_boundary(face)==false && cell_is_in_solid_domain(cell->neighbor(face))))
                    {

                        acoustic_fe_face_values.reinit(cell,face);

                        for(unsigned int i=0; i<acoustic_fe_face_values.dofs_per_cell; ++i)
                        {
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

        }

    }

    timer.stop ();
    deallog << "done (" << timer.wall_time() << "s)" << std::endl;
}


template <int dim>
void AcoustoElastic<dim>::assemble_impedance_step2_unsplit()
{
}

template <int dim>
void AcoustoElastic<dim>::assemble_impedance_step2_split()
{
}
