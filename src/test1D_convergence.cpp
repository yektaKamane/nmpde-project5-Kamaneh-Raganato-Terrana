#include <deal.II/base/convergence_table.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "FisherKolmogorov1D_convergence.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  const unsigned int               mpi_rank =
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

    // Check if the correct number of command-line arguments is provided
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " [parameter_file]" << std::endl;
    return 1;
  }



  // Read the optional parameter file name or use a default value
  std::string parameter_file;
  if (argc >= 3)
  {
    parameter_file = argv[1];
  }
  else
  {
    parameter_file = "../input/test1D_convergence.prm"; // Default parameter file
  }
  
  // const std::vector<int> N_vals = {200,
  //                                  400,
  //                                  800,
  //                                  1600};
  // const std::vector<double> h_vals = {1/200,
  //                                     1/400,
  //                                     1/800,
  //                                     1/1600};
  const std::vector<int> N_vals = {25,
                                   50,
                                   100,
                                   200};
  const std::vector<double> h_vals = {1/25,
                                      1/50,
                                      1/100,
                                      1/200};

  std::vector<double> errors_L2;
  std::vector<double> errors_H1;

  for (unsigned int i = 0; i < N_vals.size(); ++i)
    {
      // Default dimension for this 1D convergence test is 1
      FisherKol<1> problem(N_vals[i], parameter_file);

      problem.setup();
      problem.solve();

      errors_L2.push_back(problem.compute_error(VectorTools::L2_norm));
      errors_H1.push_back(problem.compute_error(VectorTools::H1_norm));
    }

  // Print the errors and estimate the convergence order.
  if (mpi_rank == 0)
    {
      std::cout << "==============================================="
                << std::endl;
      ConvergenceTable table;
      std::ofstream convergence_file("convergence.csv");
      convergence_file << "h,eL2,eH1" << std::endl;

      for (unsigned int i = 0; i < h_vals.size(); ++i)
        {
          table.add_value("h", h_vals[i]);
          table.add_value("L2", errors_L2[i]);
          table.add_value("H1", errors_H1[i]);
          
          convergence_file << h_vals[i] << "," << errors_L2[i] << ","
                           << errors_H1[i] << std::endl;

          std::cout << std::scientific << "h = " << std::setw(4)
                    << std::setprecision(2) << h_vals[i];

          std::cout << std::scientific << " | eL2 = " << errors_L2[i];

          // Estimate the convergence order.
          if (i > 0)
            {
              const double p =
                std::log(errors_L2[i] / errors_L2[i - 1]) /
                std::log(h_vals[i] / h_vals[i - 1]);

              std::cout << " (" << std::fixed << std::setprecision(2)
                        << std::setw(4) << p << ")";
            }
          else
            std::cout << " (  - )";

          std::cout << std::scientific << " | eH1 = " << errors_H1[i];

          // Estimate the convergence order.
          if (i > 0)
            {
              const double p =
                std::log(errors_H1[i] / errors_H1[i - 1]) /
                std::log(h_vals[i] / h_vals[i - 1]);

              std::cout << " (" << std::fixed << std::setprecision(2)
                        << std::setw(4) << p << ")";
            }
          else
            std::cout << " (  - )";

          std::cout << "\n";
        }

      table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);
      table.set_scientific("L2", true);
      table.set_scientific("H1", true);
      table.write_text(std::cout);
    }

  return 0;
}



