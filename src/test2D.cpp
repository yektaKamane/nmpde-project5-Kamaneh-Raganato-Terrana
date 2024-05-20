#include "FisherKolmogorov2D.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const unsigned int degree = 1;

  // const double T      = 20.0;
  const double T      = 1e-3;
  const double deltat = 1e-5;

  // FisherKolmogorov problem("../mesh/mesh-square-h0.100000.msh", degree, T, deltat);
  FisherKolmogorov problem("../mesh/mesh-square-h0.150000.msh", degree, T, deltat);

  problem.setup();
  problem.solve();

  return 0;
}


