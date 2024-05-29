#include <deal.II/base/convergence_table.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "FisherKolmogorov1D.hpp"

// Main function.
int
main(int argc, char * argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  
  const unsigned int N = 200;
  const unsigned int r = 1;
  const double T       = 20.0;
  const double deltat  = 0.1;
  std::vector<double> deltat_vals = {deltat/4, deltat/2, deltat,
                                     deltat*2, deltat*3, deltat*4};

  for (const int deltat : deltat_vals)
  {
    FisherKol problem(N, r, T, deltat);

    problem.setup();
    problem.solve();
    // non so, perché vorremmo salvare tutto, oppure trovare il valore finale di c e vedere quando converge e quando no?
    // poi qualcosa tipo switch case per i tre valori di alpha e d
  }

  return 0;
}



