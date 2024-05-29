#ifndef FISHER_KOLMOGOROV_HPP
#define FISHER_KOLMOGOROV_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// Class representing the non-linear diffusion problem.
template <int dim>
class FisherKol
{
public:
  // Physical dimension (1D, 2D, 3D)
  // static constexpr unsigned int dim = 1; // MODIFIED

  // Function for the alpha coefficient.
  class FunctionAlpha : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 2.0;
    }
  };
  
  // Function of the fiber field
  class FunctionN
  {
  public:
    Tensor<2, dim>
    isotropic(const Point<dim> & p) const
    {
      Tensor<2, dim> values;
      for (unsigned int i = 0; i < dim; ++i)
      {
        values[i][i] = 0.0;
      }
      return values;
    }

  };

  // Function of the matrix D
  class FunctionD
  {
  public:
    Tensor<2, dim> matrix_value(const Point<dim> & /*p*/ /* ,
                   Tensor<2,dim> &values */) const
    {
      Tensor<2, dim> values;
      for (unsigned int i = 0; i < dim; ++i)
      {
        values[i][i] = 0.0002;
      }
      // values[1][1] += 10.0;
      return values;
    }
  };

  // Function for the forcing term.
  class ForcingTerm : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & p,
          const unsigned int /*component*/ = 0) const override
    {
      double temp_val = std::cos(M_PI * p[0]) * std::exp(-this->get_time());

      if (dim == 1)
      {
        temp_val = std::cos(M_PI * p[0]);
        return (M_PI * M_PI - 2) * temp_val + temp_val * temp_val;
      };
      
      if (dim == 2) 
        // return ((2 * M_PI * M_PI - 1) * temp_val - 2) * std::exp(-this->get_time()) +
        //        (temp_val*temp_val + 3*temp_val + 2)* std::exp(-this->get_time() * 2);
        return -2 * M_PI * M_PI * temp_val * std::exp(-this->get_time()) -
               (temp_val * temp_val + 4 - 4 * temp_val) * std::exp(-this->get_time() * 2);

      if (dim == 3)
        return 0.0;

      else return 0.0;
    }
  };

  // Function for Dirichlet boundary conditions.
  // class FunctionG : public Function<dim>
  // {
  // public:
  //   virtual double
  //   value(const Point<dim> & /*p*/,
  //         const unsigned int /*component*/ = 0) const override
  //   {
  //     return 0.0;
  //   }
  // };

  // Function for Neumann boundary condition.
  class FunctionH : public Function<dim>
  {
  public:
    // Constructor.
    FunctionH()
    {}

    // Evaluation.
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 0.0;
    }
  };

  // Function for initial conditions.
  class FunctionU0 : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & p,
          const unsigned int /*component*/ = 0) const override
    {
      // if (p[0] < 0.55 && p[0] > 0.45 && p[1] < 0.55 && p[1] > 0.45 && p[2] < 0.55 && p[2] > 0.45)
      // // if (p[0] < 0.05 && p[0] > -0.05 && p[1] < 0.05 && p[1] > -0.05 && p[2] < 0.05 && p[2] > -0.05)
      // // if(p[0] == 0)
      // {
      //   return 0.1;
      // }
      
      // if (p[0] < 80.0 && p[0] > 70.0 && p[1] < 95.0 && p[1] > 90.0 /*&& p[2] < 50.0 && p[2] > 40.0*/)
      // if (p[0] < 0.55 && p[0] > 0.45 && p[1] < 0.55 && p[1] > 0.45)
      // {
      //   return 0.3;
      // }
      // if (p[0]>=0.45 && p[0]<=0.55)
      //   return 0.1;
      return std::cos(M_PI * p[0]);
  
      // return 0.0;
      // return p[0] * (1 - p[0]) * p[1] * (1 - p[1]);
    }
  };

  // Exact solution.
  class ExactSolution : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      return std::cos(M_PI*p[0]) * std::exp(-this->get_time());
    }

    virtual Tensor<1, dim>
    gradient(const Point<dim> &p,
             const unsigned int /*component*/ = 0) const override
    {
      Tensor<1, dim> result;

      // duex / dx
      result[0] = -M_PI * std::sin(M_PI * p[0]) * std::exp(-this->get_time());

      return result;
    }
  };

  // Constructor. We provide the final time, time step Delta t and theta method
  // parameter as constructor arguments.
  FisherKol(const unsigned int N_,
                const unsigned int &r_,
                const double       &T_,
                const double       &deltat_,
                const std::string  &prm_file_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , T(T_)
    , N(N_)
    , r(r_)
    , deltat(deltat_)
    , prm_file(prm_file_)
    , mesh(MPI_COMM_WORLD)
  {
      parameters.declare_entry("coef_alpha", "1.0", Patterns::Double(), "dummy");
      parameters.declare_entry("coef_dext", "1.0", Patterns::Double(), "dummy");
      parameters.declare_entry("coef_daxn", "0.0", Patterns::Double(), "dummy");
      parameters.declare_entry("fib", "0", Patterns::Integer(), "dummy");

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

  // alpha coefficient.
  FunctionAlpha alpha;

  // ...
  FunctionN fiber;

  // matrix D.
  FunctionD D;

  // Forcing term.
  ForcingTerm forcing_term;

  // Dirichlet boundary conditions.
  // FunctionG function_g;

  // Neumann boundary condition.
  FunctionH function_h;

  // Initial conditions.
  FunctionU0 u_0;

  // Exact solution.
  ExactSolution exact_solution;

  // Current time.
  double time;

  // Final time.
  const double T;

  // Number of elements.
  const unsigned int N;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  // const std::string mesh_file_name;

  // Polynomial degree.
  const unsigned int r;

  // Time step.
  const double deltat;

  const std::string prm_file;

  ParameterHandler parameters;

  // Mesh.
  parallel::fullydistributed::Triangulation<dim> mesh;
  // Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // Quadrature formula used on boundary lines.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_boundary;

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