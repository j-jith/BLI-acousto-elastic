#include "AcoustoElastic.h"

int main()
{
   Timer timer;
   timer.start ();

   //deallog.attach(std::cout);

   ParameterHandler  prm;
   ParameterReader   param(prm);
   param.read_parameters("params.prm");

   bool assemble_and_dump=false;
   AcoustoElastic<3> ac_el(prm, 2, assemble_and_dump);
   ac_el.run(0);

   timer.stop ();
   deallog << "\n\n***All done!*** (" << timer.wall_time() << "s)" << std::endl;
}
