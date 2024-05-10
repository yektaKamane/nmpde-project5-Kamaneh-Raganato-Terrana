#include "Fisher_Kolmogorov_solver.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  // const unsigned int               mpi_rank =
  //   Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  const unsigned int degree = 1;

  const double T      = 2.0;
  const double deltat = 0.1;
  const unsigned int dim = 2;

  // ParameterHandler parameter_handler;
  // FisherKol::Parameters parameters;
  // parameters.declare_parameters(parameter_handler);
  // parameter_handler.parse_input("../input/test1.prm");
  // parameters.get_parameters(parameter_handler);

  // FisherKol problem("../mesh/mesh-cube-20.msh", degree, T, deltat);
  // FisherKol problem("../mesh/brain-h3.0.msh", degree, T, deltat);
  FisherKol<dim> problem("../mesh/mesh-square-40.msh", degree, T, deltat);

  problem.setup();
  problem.solve();

  return 0;
}