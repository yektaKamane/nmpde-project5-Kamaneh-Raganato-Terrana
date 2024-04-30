#include "Fisher_Kolmogorov_solver.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const unsigned int degree = 1;

  // const double T      = 20.0;
  // const double T      = 4.9;
  const double T      = 0.9;
  const double deltat = 0.1;

  // FisherKol problem("../mesh/brain_3.msh", degree, T, deltat);
  // FisherKol problem("../mesh/brain_benchmark.msh", degree, T, deltat);
  FisherKol problem("../mesh/brain-h3.0_diretta.msh", degree, T, deltat);
  // FisherKol problem("../mesh/mesh-cube-5.msh", degree, T, deltat);
  // FisherKol problem("../mesh/mesh-cube-40.msh", degree, T, deltat);
  // FisherKol problem("../mesh/mesh-cube-10.msh", degree, T, deltat);

  problem.setup();
  problem.solve();

  return 0;
}