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

   // set point load
   Point<3> force_loc(0.120,0.0,0.011);
   Point<3> force_dir(0,0,1);
   double force_mag = 1e3;
   ac_el.set_point_load(force_loc, force_dir, force_mag);

   ac_el.run(ACOUSTO_ELASTIC_STEP);

   timer.stop ();
   deallog << "\n\n***All done!*** (" << timer.wall_time() << "s)" << std::endl;
}
