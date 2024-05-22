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
  const double       T = 20.0;
  const double deltat  = 0.1;

  FisherKol problem(N, r, T, deltat);

  problem.setup();
  problem.solve();

  return 0;
}

// Main function.
// int
// main(int /*argc*/, char * /*argv*/[])
// {
//   ConvergenceTable table;

//   // const std::vector<unsigned int> N_values = {9, 19, 39, 79, 159, 319};
//   const unsigned int N      = 200;
//   const unsigned int degree = 1;

//   std::ofstream convergence_file("convergence.csv");
//   convergence_file << "h,eL2,eH1" << std::endl;

//   for (const unsigned int &N : N_values)
//     {
//       FisherKol problem(N, degree);

//       problem.setup();
//       problem.assemble();
//       problem.solve();
//       problem.output();

//       const double h        = 1.0 / (N + 1.0);
//       const double error_L2 = problem.compute_error(VectorTools::L2_norm);
//       const double error_H1 = problem.compute_error(VectorTools::H1_norm);

//       table.add_value("h", h);
//       table.add_value("L2", error_L2);
//       table.add_value("H1", error_H1);

//       convergence_file << h << "," << error_L2 << "," << error_H1 << std::endl;
//     }

//   table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);

//   table.set_scientific("L2", true);
//   table.set_scientific("H1", true);

//   table.write_text(std::cout);

//   return 0;
// } 

// Main function.
// int
// main(int argc, char *argv[])
// {
//   Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
//   const unsigned int               mpi_rank =
//     Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

//   const unsigned int N      = 200;
//   const unsigned int degree = 1;

//   const double T      = 1e-3;
//   const double deltat = 1e-5;
//   // const double theta  = 0.5;

//   std::vector<double> errors_L2;
//   std::vector<double> errors_H1;

//   for (unsigned int i = 0; i < meshes.size(); ++i)
//     {
//       FisherKol problem(meshes[i], degree, T, deltat);

//       problem.setup();
//       problem.solve();

//       errors_L2.push_back(problem.compute_error(VectorTools::L2_norm));
//       errors_H1.push_back(problem.compute_error(VectorTools::H1_norm));
//     }

//   // Print the errors and estimate the convergence order.
//   if (mpi_rank == 0)
//     {
//       std::cout << "==============================================="
//                 << std::endl;
//       ConvergenceTable table;
//       std::ofstream convergence_file("convergence.csv");
//       convergence_file << "h,eL2,eH1" << std::endl;

//       for (unsigned int i = 0; i < h_vals.size(); ++i)
//         {
//           table.add_value("h", h_vals[i]);
//           table.add_value("L2", errors_L2[i]);
//           table.add_value("H1", errors_H1[i]);
          
//           convergence_file << h_vals[i] << "," << errors_L2[i] << ","
//                            << errors_H1[i] << std::endl;

//           std::cout << std::scientific << "h = " << std::setw(4)
//                     << std::setprecision(2) << h_vals[i];

//           std::cout << std::scientific << " | eL2 = " << errors_L2[i];

//           // Estimate the convergence order.
//           if (i > 0)
//             {
//               const double p =
//                 std::log(errors_L2[i] / errors_L2[i - 1]) /
//                 std::log(h_vals[i] / h_vals[i - 1]);

//               std::cout << " (" << std::fixed << std::setprecision(2)
//                         << std::setw(4) << p << ")";
//             }
//           else
//             std::cout << " (  - )";

//           std::cout << std::scientific << " | eH1 = " << errors_H1[i];

//           // Estimate the convergence order.
//           if (i > 0)
//             {
//               const double p =
//                 std::log(errors_H1[i] / errors_H1[i - 1]) /
//                 std::log(h_vals[i] / h_vals[i - 1]);

//               std::cout << " (" << std::fixed << std::setprecision(2)
//                         << std::setw(4) << p << ")";
//             }
//           else
//             std::cout << " (  - )";

//           std::cout << "\n";
//         }

//       table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);
//       table.set_scientific("L2", true);
//       table.set_scientific("H1", true);
//       table.write_text(std::cout);
//     }

//   return 0;
// }