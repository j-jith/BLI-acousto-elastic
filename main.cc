#include "AcoustoElastic.h"

int main()
{
   Timer timer;
   timer.start ();

   //deallog.attach(std::cout);

   ParameterHandler  prm;
   ParameterReader   param(prm);
   param.read_parameters("params.prm");

   const unsigned int dim = 3;
   unsigned int degree = prm.get_integer("degree_shape");

   bool assemble_and_dump=false;
   AcoustoElastic<dim> ac_el(prm, degree, assemble_and_dump);

   ac_el.run();

   timer.stop ();
   deallog << "\n\n***All done!*** (" << timer.wall_time() << "s)" << std::endl;
}
