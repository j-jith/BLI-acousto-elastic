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
void AcoustoElastic<dim>::read_params()
{
    f_0=prm.get_double("f_0");
    f_1=prm.get_double("f_1");
    n_f=prm.get_integer("n_f");

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


    solution_step = prm.get_integer("solution_step");

    if(solution_step == BLI_STEP_1_APPROX || solution_step == BLI_STEP_2_APPROX)
    {
        fit_c2=prm.get_double("fit_c2");
        fit_c1=prm.get_double("fit_c1");
        fit_c0=prm.get_double("fit_c0");
    }

    point_load = prm.get_bool("point_load");

    std::vector<std::string> point_load_info;
    if(point_load)
    {
        point_load_info.clear();
        point_load_info = Utilities::split_string_list(prm.get("point_load_loc"));
        for(unsigned int dim_i = 0; dim_i < dim; dim_i++)
        {
            point_load_loc(dim_i) = ::atof(point_load_info[dim_i].c_str());
        }
        point_load_info.clear();
        point_load_info = Utilities::split_string_list(prm.get("point_load_dir"));
        for(unsigned int dim_i = 0; dim_i < dim; dim_i++)
        {
            point_load_dir(dim_i) = ::atof(point_load_info[dim_i].c_str());
        }
        point_load_info.clear();

        point_load_mag = prm.get_double("point_load_mag");
    }

    N_pod=prm.get_integer("N_pod");
    pod_thresh=prm.get_double("pod_thresh");

    //Output parameters for troubleshooting
    std::ofstream prm_out("prm_out.txt");
    prm.print_parameters(prm_out, ParameterHandler::Text);

}

template <int dim>
void AcoustoElastic<dim>::run()
{
    deallog << "Initialising and reading grid... ";

    make_grid();

    read_params();

    if(f_1>f_0)
    {
        sweep=true;
    }
    else
    {
        omega=f_0*2*M_PI;
        sweep=false;
    }


    setup_system();
    assemble_system();

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
                deallog <<"Frequency: " <<omega/2/M_PI<<" Hz ["<<n_sweep+1<<"/"<<n_f<<"]"<<std::endl;

                assemble_impedance();

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
            assemble_impedance();
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
