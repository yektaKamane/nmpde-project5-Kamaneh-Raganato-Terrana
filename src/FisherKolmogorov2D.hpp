#ifndef FK_2D_HPP
#define FK_2D_HPP

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

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/base/tensor_function.h> // AGGIUNTO

#include <fstream>
#include <iostream>

using namespace dealii;

// Class representing the Fisher-Kolmogorov problem.
class FisherKolmogorov
{
public:
  // Physical dimension (1D, 2D, 3D).
  static constexpr unsigned int dim = 2;

  // Function for the alpha coefficient, that represents the growth of c,
  // the relative concentration of misfolded protein.
  class FunctionAlpha : public Function<dim>
  {
  public:
    // Constructor.
    FunctionAlpha()
    {}

    // Evaluation.
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 1.0;
    }
  };

  // Function for the D coefficient, that represents the spreading of c,
  // the relative concentration of misfolded protein.
  // class FunctionD : public Function<dim>
  // {
  // public:
  //   // Constructor.
  //   FunctionD()
  //   {}

  //   // Evaluation.
  //   virtual double
  //   value(const Point<dim> & /*p*/,
  //         const unsigned int /*component*/ = 0) const override
  //   {
  //     return 0.0001;
  //   }
  // };

  // https://www.dealii.org/current/doxygen/deal.II/step_21.html
  // https://www.dealii.org/current/doxygen/deal.II/step_61.html
  // unit_symmetric_tensor() https://www.dealii.org/current/doxygen/deal.II/classSymmetricTensor.html#ae3f87e0ac60ba27e61941451b7185626
  // outer_product() https://www.dealii.org/current/doxygen/deal.II/symmetric__tensor_8h.html#a548c077df22f64aae4000c382c253dce
  class TensorFunctionD : public TensorFunction<2, dim>
  {
  public:
    // Constructor.
    TensorFunctionD()
      : TensorFunction<2, dim>()
    {}

    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<Tensor<2, dim>> &  values) const override
    {
      AssertDimension(points.size(), values.size());

      for (unsigned int p = 0; p < points.size(); ++p)
        {
          values[p].clear();

          // Extracellular diffusion contribution.
          values[p] = d_ext * unit_symmetric_tensor<dim>();

          // Axonal transport contribution.
          Tensor<1, dim> n;
          for (unsigned int d = 0; d < dim; ++d)
            n[d] = 1.0;
          values[p] += d_axn * outer_product(n, n);
        }
    }

  protected:
    // Extracellular diffusion coefficient, associated with isotropic diffusion
    // of misfolded protein through the extracellular space.
    const double d_ext = 0.0001; // TO FIX

    // Axonal transport coefficient, associated with anisotropic diffusion
    // of misfolded protein along the local axonal direction.
    const double d_axn = 0.0001; // TO FIX
  };

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

  // TENERE o BUTTAR VIA?
  // Function for the forcing term.
  class ForcingTerm : public Function<dim>
  {
  public:
    // Constructor.
    ForcingTerm()
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
  class FunctionC0 : public Function<dim>
  {
  public:
    // Constructor.
    FunctionC0()
    {}

    // Evaluation.
    virtual double
    value(const Point<dim> & p,
          const unsigned int /*component*/ = 0) const override
    {
      if (p[0] < 0.55 && p[0] > 0.45 && p[1] < 0.55 && p[1] > 0.45)
      {
        return 0.8;
      }
      return 0.0;
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
      return (std::cos(M_PI*p[0]) * std::cos(M_PI*p[1]) + 2) * std::exp(-get_time());
    }

    virtual Tensor<1, dim>
    gradient(const Point<dim> &p,
             const unsigned int /*component*/ = 0) const override
    {
      Tensor<1, dim> result;

      // duex / dx
      result[0] = -M_PI * std::sin(M_PI * p[0]) * std::cos(M_PI * p[1]) * std::exp(-get_time());

      // duex / dy
      result[1] = -M_PI * std::cos(M_PI * p[0]) * std::sin(M_PI * p[1]) * std::exp(-get_time());

      return result;
    }
  };

  // Constructor.
  FisherKolmogorov(const std::string  &mesh_file_name_,
                   const unsigned int &r_,
                   const double       &T_,
                   const double       &deltat_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , T(T_)
    , mesh_file_name(mesh_file_name_)
    , r(r_)
    , deltat(deltat_)
    , mesh(MPI_COMM_WORLD)
  {}

  // Initialization.
  void
  setup();

  // Solve the problem.
  void
  solve();

  // Compute the error.
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

  // D coefficient.
  TensorFunctionD D;

  // Neumann boundary condition.
  FunctionH function_h;

  // Forcing term.
  // ForcingTerm forcing_term;

  // Initial conditions.
  FunctionC0 c_0;

  // Exact solution.
  ExactSolution exact_solution;

  // Current time.
  double time;

  // Final time.
  const double T;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  const std::string mesh_file_name;

  // Polynomial degree.
  const unsigned int r;

  // Time step.
  const double deltat;

  // Mesh.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // Quadrature formula on boundary lines.
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