#include "FisherKolmogorov.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const unsigned int degree = 1;

  const double T      = .0;
  const double deltat = 0.1;

  FisherKolmogorov problem("../mesh/mesh-cube-20.msh", degree, T, deltat);

  problem.setup();
  problem.solve();

  return 0;
}