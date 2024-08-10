#ifndef FISHER_KOLMOGOROV_HPP
#define FISHER_KOLMOGOROV_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/base/parameter_handler.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// Class representing the non-linear diffusion problem.
template <int dim>
class FisherKol
{
public:

  // Function of the fiber field
  class FunctionN
  {
  public:
    Tensor<2, dim>
    isotropic(const Point<dim> & /*p*/) const
    {
      Tensor<2, dim> values;
      for (unsigned int i = 0; i < dim; ++i)
      {
        values[i][i] = 0.0;
      }
      return values;
    }

  };

  // Function for initial conditions.
  class FunctionU0 : public Function<dim>
  {
  public:

    /* NEW */
    FunctionU0(const FisherKol<dim> &fisher_kol) : fisher_kol(fisher_kol) {}

    virtual double
    value(const Point<dim> & p,
          const unsigned int /*component*/ = 0) const override
    {
      const double radius = fisher_kol.parameters.get_double("radius");
      const double center_x = fisher_kol.parameters.get_double("center_x");
      const double center_y = fisher_kol.parameters.get_double("center_y");
      const double temp = (p[0] - center_x)*(p[0] - center_x) + (p[1] - center_y)*(p[1] - center_y);
      double distance = 0.0;

      if (dim == 2){
        distance = std::sqrt(temp);
      }

      if (dim == 3){
        const double center_z = fisher_kol.parameters.get_double("center_z");
        distance = std::sqrt(temp + (p[2] - center_z)*(p[2] - center_z));
      }

      if (distance <= radius)
        return 0.95;
      else
        return 0.0;
    }
    
    private:
      const FisherKol<dim> &fisher_kol;  // Added member to store reference to FisherKol
  
  };

   // Exact solution
  class ExactSolution : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      double temp_val = std::cos(M_PI * p[0]) * std::cos(M_PI * p[1]);
      if (dim == 2) 
        return (temp_val + 2) * std::exp(-this->get_time());

      if (dim == 3)
        return (temp_val * std::cos(M_PI * p[2])) * std::exp(-this->get_time());

      else return 0.0;
      
    }

    virtual Tensor<1, dim>
    gradient(const Point<dim> &p,
             const unsigned int /*component*/ = 0) const override
    {
      Tensor<1, dim> result;

      if (dim == 2)
      {
        result[0] = -M_PI * std::sin(M_PI * p[0]) * std::cos(M_PI * p[1]) *
                    std::exp(-this->get_time());
        result[1] = -M_PI * std::cos(M_PI * p[0]) * std::sin(M_PI * p[1]) *
                    std::exp(-this->get_time());
      }
      
      if (dim == 3)
      {
        result[0] = -M_PI * std::sin(M_PI * p[0]) * std::cos(M_PI * p[1]) *
                    std::cos(M_PI * p[2]) * std::exp(-this->get_time());
        result[1] = -M_PI * std::sin(M_PI * p[1]) * std::cos(M_PI * p[0]) *
                    std::cos(M_PI * p[2]) * std::exp(-this->get_time());
        result[2] = -M_PI * std::sin(M_PI * p[2]) * std::cos(M_PI * p[0]) *
                    std::cos(M_PI * p[1]) * std::exp(-this->get_time());
      }

      return result;
    }
  };

  // Constructor. We provide the final time, time step Delta t and theta method
  // parameter as constructor arguments.
  FisherKol(const std::string  &mesh_file_name_,
                // const unsigned int &r_,
                // const double       &T_,
                // const double       &deltat_,
                const std::string  &prm_file_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    // , T(T_)
    , mesh_file_name(mesh_file_name_)
    // , r(r_)
    // , deltat(deltat_)
    , prm_file(prm_file_)
    , mesh(MPI_COMM_WORLD)
    /*NEW*/
    , u_0(*this) // Changed: Pass reference of this to u_0
  {
      parameters.declare_entry("coef_alpha", "1.0", Patterns::Double(), "dummy");
      parameters.declare_entry("coef_dext", "1.0", Patterns::Double(), "dummy");
      parameters.declare_entry("coef_daxn", "1.0", Patterns::Double(), "dummy");
      parameters.declare_entry("fib", "0", Patterns::Integer(), "dummy");

      parameters.declare_entry("T", "0", Patterns::Double(), "dummy");
      parameters.declare_entry("deltat", "0", Patterns::Double(), "dummy");
      parameters.declare_entry("degree", "0", Patterns::Integer(), "dummy");

      parameters.declare_entry("radius", "10.0", Patterns::Double(), "dummy");
      parameters.declare_entry("center_x", "0.0", Patterns::Double(), "dummy");
      parameters.declare_entry("center_y", "0.0", Patterns::Double(), "dummy");
      parameters.declare_entry("center_z", "0.0", Patterns::Double(), "dummy");

      parameters.parse_input(prm_file);
  }

  // Initialization.
  void
  setup();

  // Solve the problem.
  void
  solve();

  // Compute the error for convergence analysis.
  double
  compute_error(const VectorTools::NormType &norm_type);

protected:
  // Assemble the tangent problem.
  void
  assemble_system();

  // Solve the linear system associated to the tangent problem.
  void
  solve_linear_system();

  // Solve the problem for one time step using Newton's method.
  void
  solve_newton();

  // Output.
  void
  output(const unsigned int &time_step) const;


  // MPI parallel. /////////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Problem definition. ///////////////////////////////////////////////////////

  FunctionN fiber;

  // Initial conditions.
  FunctionU0 u_0;

  // Exact solution.
  ExactSolution exact_solution;

  // Current time.
  double time;

  // // Final time.
  // const double T;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  const std::string mesh_file_name;

  // Polynomial degree.
  // const unsigned int r;

  // Time step.
  // const double deltat;

  const std::string prm_file;

  ParameterHandler parameters;

  // Mesh.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;

  // DoFs relevant to the current process (including ghost DoFs).
  IndexSet locally_relevant_dofs;

  // Jacobian matrix.
  TrilinosWrappers::SparseMatrix jacobian_matrix;

  // Residual vector.
  TrilinosWrappers::MPI::Vector residual_vector;

  // Increment of the solution between Newton iterations.
  TrilinosWrappers::MPI::Vector delta_owned;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::Vector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::Vector solution;

  // System solution at previous time step.
  TrilinosWrappers::MPI::Vector solution_old;
};

#endif