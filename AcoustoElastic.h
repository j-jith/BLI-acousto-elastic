#ifndef ACOUSTO_ELASTIC_H
#define ACOUSTO_ELASTIC_H

#include "includes.h"
#include "my_point_source_vector.h"
#include "ParameterReader.h"

#define ACOUSTO_ELASTIC_STEP    0
#define BLI_STEP_1              1
#define BLI_STEP_2              2

#define SOLID_DOMAIN            101
#define FLUID_DOMAIN            102
#define SOLID_FIXED_BOUNDARY    201


template <int dim>
class AcoustoElastic
{
   public:
   AcoustoElastic(ParameterHandler &, const unsigned int, bool);
   void run(unsigned int step);

   private:
   enum
   {
      solid_domain_id,
      fluid_domain_id
   };

   const unsigned int degree;

   bool assemble_and_dump;

   static bool
   cell_is_in_solid_domain(const typename hp::DoFHandler<dim>::cell_iterator &cell); 

   static bool
   cell_is_in_fluid_domain(const typename hp::DoFHandler<dim>::cell_iterator &cell); 

   void make_grid();
   void set_active_fe_indices();
   void my_sparsity_pattern();
   void setup_system();
   void read_input(unsigned int);
   void assemble_system(bool step1);
   void assemble_system_aux();
   void read_pod();
   void pod_aux();
   void dump_matrices();
   void solve();
   void solve_reduced();
   void output_results(unsigned int num);
   void output_binary(unsigned int num);
   void output_ktr(unsigned int num);

   ParameterHandler &prm;

   std::vector<Vector<double> > pod_modes;
   unsigned int N_pod;
   double pod_thresh;

   double omega;
   double f_0, f_1;
   unsigned int n_f;
   bool sweep;
   
   double E_s, nu_s, rho_s;
   double lambda_s, mu_s;

   double rho_0, c_0, mu, kappa, cp, gamma;

   double fit_c2, fit_c1, fit_c0;

   Triangulation<dim> triangulation;
   FESystem<dim> elastic_fe;
   FESystem<dim> acoustic_fe;
   hp::FECollection<dim> fe_collection;
   hp::DoFHandler<dim> dof_handler;

   ConstraintMatrix constraints;

   SparsityPattern sparsity_pattern;
   SparseMatrix<double> system_matrix;
   SparseMatrix<double> system_matrix_aux;
   SparseMatrix<double> mass_matrix;
   SparseMatrix<double> damp_matrix;
   SparseMatrix<double> stif_matrix;
   Vector<double> system_rhs;
   Vector<double> solution;
   Vector<double> global_ktr;

   FullMatrix<double> red_mass_matrix;
   FullMatrix<double> red_damp_matrix;
   FullMatrix<double> red_stif_matrix;
   FullMatrix<double> red_matrix_aux;
   Vector<double> red_rhs;
};

#endif
