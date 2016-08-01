#include "ParameterReader.h"

// constructor
ParameterReader::ParameterReader(ParameterHandler &paramhandler) : prm(paramhandler) {}

void ParameterReader::declare_parameters()
{
  //prm.enter_subsection ("Physical constants");
  //{
      prm.declare_entry("f_0", "100",
                        Patterns::Double(0),
                        "Initial Frequency");

      prm.declare_entry("f_1", "100",
                        Patterns::Double(0),
                        "Final Frequency");

      prm.declare_entry("n_f", "1",
                        Patterns::Integer(0),
                        "No. of Frequencies");

      prm.declare_entry("E_s", "2e11",
                        Patterns::Double(0),
                        "Young's modulus");

      prm.declare_entry("nu_s", "0.3",
                        Patterns::Double(0),
                        "Poisson's ratio");

      prm.declare_entry("rho_s", "7800",
                        Patterns::Double(0),
                        "Density of solid material");

      prm.declare_entry("rho_0", "286.13",
                        Patterns::Double(0),
                        "Density of fluid material");

      prm.declare_entry("c_0", "283.76",
                        Patterns::Double(0),
                        "Speed of sound in fluid");

      prm.declare_entry("mu", "5.5459e-5",
                        Patterns::Double(0),
                        "Viscosity of fluid");

      prm.declare_entry("kappa", "0.0745",
                        Patterns::Double(0),
                        "Thermal conductivity of fluid");

      prm.declare_entry("cp", "863.6",
                        Patterns::Double(0),
                        "Isobaric heat capacity of fluid");

      prm.declare_entry("gamma", "1.28",
                        Patterns::Double(0),
                        "Ratio of heat capacities of fluid");

      prm.declare_entry("fit_c2", "0.002533",
                        Patterns::Double(0),
                        "2nd order coefficient of polynomial fit of omega*sqrt(omega)");
      prm.declare_entry("fit_c1", "95.49",
                        Patterns::Double(0),
                        "1st order coefficient of polynomial fit of omega*sqrt(omega)");
      prm.declare_entry("fit_c0", "0",
                        Patterns::Double(),
                        "constant coefficient of polynomial fit of omega*sqrt(omega)");
      prm.declare_entry("N_pod", "1",
                        Patterns::Integer(0),
                        "No. of POD input solutions");
      prm.declare_entry("pod_thresh", "1e-12",
                        Patterns::Double(0),
                        "POD threshold for truncation");
    //}
    //prm.leave_subsection ();
}

void ParameterReader::read_parameters (const std::string parameter_file)
{
    declare_parameters();

    prm.read_input (parameter_file);
}


