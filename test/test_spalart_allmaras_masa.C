//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-


#include "grins_config.h"

#ifdef GRINS_HAVE_MASA

#include <iostream>

// GRINS
#include "grins/mesh_builder.h"
#include "grins/simulation.h"
#include "grins/simulation_builder.h"
#include "grins/multiphysics_sys.h"
#include "grins/parabolic_profile.h"

//libMesh
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/mesh_function.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/gmv_io.h"
#include "libmesh/exact_solution.h"

//MASA
#include "masa.h"

// GRVY
#ifdef GRINS_HAVE_GRVY
#include "grvy.h"
#endif

libMesh::Real
initial_values( const libMesh::Point& p, const libMesh::Parameters &params,
		const std::string& system_name, const std::string& unknown_name );

libMesh::Number
exact_solution( const libMesh::Point& p,
		const libMesh::Parameters&,   // parameters, not needed
		const std::string&,  // sys_name, not needed
		const std::string&); // unk_name, not needed);

libMesh::Gradient
exact_derivative( const libMesh::Point& p,
		  const libMesh::Parameters&,   // parameters, not needed
		  const std::string&,  // sys_name, not needed
		  const std::string&); // unk_name, not needed);

class MasaBCFactory : public GRINS::BoundaryConditionsFactory
{
public:

  MasaBCFactory( )
    : GRINS::BoundaryConditionsFactory()
  { return; };

  ~MasaBCFactory(){return;};

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > build_dirichlet( );
};

// Class to construct the Dirichlet boundary object and operator for the MASA specified u, v velocity and pressure components
class MasaBdyFunctionU : public libMesh::FunctionBase<libMesh::Number>
{
public:
  MasaBdyFunctionU ()
  { this->_initialized = true; }

  virtual libMesh::Number operator() (const libMesh::Point&, const libMesh::Real = 0)
  { libmesh_not_implemented(); }

  virtual void operator() (const libMesh::Point& p,
                           const libMesh::Real t,
                           libMesh::DenseVector<libMesh::Number>& output)
  {
    output.resize(2);
    output.zero();

    // The x-component
    output(0) = MASA::masa_eval_exact_u<libMesh::Real>  ( p(0), p(1) );
    // The y-component
    output(1) = MASA::masa_eval_exact_v<libMesh::Real>  ( p(0), p(1) );
    // The pressure
    //output(2) = MASA::masa_eval_exact_p<libMesh::Real>  ( p(0), p(1) );

  }

  virtual libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > clone() const
  { return libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > (new MasaBdyFunctionU()); }

};

// Class to construct the Dirichlet boundary object and operator for the MASA specified nu
class MasaBdyFunctionNu : public libMesh::FunctionBase<libMesh::Number>
{
public:
  MasaBdyFunctionNu ()
  { this->_initialized = true; }

  virtual libMesh::Number operator() (const libMesh::Point&, const libMesh::Real = 0)
  { libmesh_not_implemented(); }

  virtual void operator() (const libMesh::Point& p,
                           const libMesh::Real t,
                           libMesh::DenseVector<libMesh::Number>& output)
  {
    output.resize(1);
    output.zero();

    // The turbulent viscosity
    output(0) = MASA::masa_eval_exact_nu<libMesh::Real>  ( p(0), p(1) );

  }

  virtual libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > clone() const
  { return libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> > (new MasaBdyFunctionNu()); }

};

int main(int argc, char* argv[])
{

#ifdef GRINS_USE_GRVY_TIMERS
  GRVY::GRVY_Timer_Class grvy_timer;
  grvy_timer.Init("GRINS Timer");
#endif

  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify libMesh input file." << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  // libMesh input file should be first argument
  std::string libMesh_input_filename = argv[1];

  // Create our GetPot object.
  GetPot libMesh_inputfile( libMesh_input_filename );

#ifdef GRINS_USE_GRVY_TIMERS
  grvy_timer.BeginTimer("Initialize Solver");
#endif

  // Initialize libMesh library.
  libMesh::LibMeshInit libmesh_init(argc, argv);

  GRINS::SimulationBuilder sim_builder;

  std::tr1::shared_ptr<MasaBCFactory> bc_factory( new MasaBCFactory );

  sim_builder.attach_bc_factory(bc_factory);

  GRINS::Simulation grins( libMesh_inputfile,
			   sim_builder,
                           libmesh_init.comm() );

  // Assign initial conditions
  std::string system_name = libMesh_inputfile( "screen-options/system_name", "GRINS" );
  std::tr1::shared_ptr<libMesh::EquationSystems> es = grins.get_equation_system();
  const libMesh::System& system = es->get_system(system_name);
  libMesh::Parameters &params = es->parameters;

  system.project_solution( initial_values, NULL, params );

#ifdef GRINS_USE_GRVY_TIMERS
  grvy_timer.EndTimer("Initialize Solver");

  // Attach GRVY timer to solver
  grins.attach_grvy_timer( &grvy_timer );
#endif

  // Solve
  grins.run();

  // Create Exact solution object and attach exact solution quantities
  libMesh::ExactSolution exact_sol(*es);

  exact_sol.attach_exact_value(&exact_solution);
  //exact_sol.attach_exact_deriv(&exact_derivative);

  // Compute error and get it in various norms
  exact_sol.compute_error("GRINS", "u");

  double l2error_u = exact_sol.l2_error("GRINS", "u");
  double h1error_u = 0.0; //exact_sol.h1_error("GRINS", "u");

  const double errortol = 1.0e-10;

  int return_flag = 0;

  if( l2error_u > errortol || h1error_u > errortol )
    {
      return_flag = 1;

      std::cout << "Tolerance exceeded for velocity in SA shear test." << std::endl
		<< "l2 error = " << l2error_u << std::endl
		<< "h1 error = " << h1error_u << std::endl;
    }

  // Compute error and get it in various norms
  exact_sol.compute_error("GRINS", "v");

  double l2error_v = exact_sol.l2_error("GRINS", "v");
  double h1error_v = 0.0; //exact_sol.h1_error("GRINS", "v");

  if( l2error_v > errortol || h1error_v > errortol )
    {
      return_flag = 1;

      std::cout << "Tolerance exceeded for velocity in SA shear test." << std::endl
		<< "l2 error = " << l2error_v << std::endl
		<< "h1 error = " << h1error_v << std::endl;
    }

  // Compute error and get it in various norms
  exact_sol.compute_error("GRINS", "p");

  double l2error_p = exact_sol.l2_error("GRINS", "p");
  double h1error_p = 0.0; //exact_sol.h1_error("GRINS", "p");

  if( l2error_p > errortol || h1error_p > errortol )
    {
      return_flag = 1;

      std::cout << "Tolerance exceeded for pressure in SA shear test." << std::endl
		<< "l2 error = " << l2error_p << std::endl
		<< "h1 error = " << h1error_p << std::endl;
    }

  // Compute error and get it in various norms
  //std::cout<<"Not computing nu error"<<std::endl;
  exact_sol.compute_error("GRINS", "nu");

  double l2error_nu = exact_sol.l2_error("GRINS", "nu");
  double h1error_nu = 0.0; //exact_sol.h1_error("GRINS", "nu");

  if( l2error_nu > errortol || h1error_nu > errortol )
    {
      return_flag = 1;

      std::cout << "Tolerance exceeded for turbulent viscosity in SA shear test." << std::endl
  		<< "l2 error = " << l2error_nu << std::endl
  		<< "h1 error = " << h1error_nu << std::endl;
    }

  return return_flag;
}

std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > MasaBCFactory::build_dirichlet( )
{
  std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> > uvp_func( new MasaBdyFunctionU );
  std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Number> > nu_func( new MasaBdyFunctionNu );

  GRINS::DBCContainer cont_U;
  cont_U.add_var_name( "u" );
  cont_U.add_var_name( "v" );
  //cont_U.add_var_name( "p" );
  cont_U.add_bc_id( 0 );
  cont_U.add_bc_id( 1 );
  cont_U.add_bc_id( 2 );
  cont_U.add_bc_id( 3 );

  cont_U.set_func( uvp_func );

  GRINS::DBCContainer cont_nu;
  cont_nu.add_var_name( "nu" );
  cont_nu.add_bc_id( 0 );
  cont_nu.add_bc_id( 1 );
  cont_nu.add_bc_id( 2 );
  cont_nu.add_bc_id( 3 );

  cont_nu.set_func( nu_func );

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > mymap;

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::incompressible_navier_stokes,  cont_U) );

  //std::cout<<"Not inserting Spalart Allmaras boundary condition map."<<std::endl;
  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::spalart_allmaras,  cont_nu) );

  return mymap;
}

libMesh::Real
initial_values( const libMesh::Point& p, const libMesh::Parameters &params,
		const std::string& , const std::string& unknown_name )
{
  libMesh::Real value = 0.0;

  const double x = p(0);
  const double y = p(1);

  if( unknown_name == "u" ) value = MASA::masa_eval_exact_u<libMesh::Real>  ( x, y );
  if( unknown_name == "v" ) value = MASA::masa_eval_exact_v<libMesh::Real>  ( x, y );
  if( unknown_name == "p" ) value = MASA::masa_eval_exact_p<libMesh::Real>  ( x, y );
  if( unknown_name == "nu" ) value = MASA::masa_eval_exact_nu<libMesh::Real>  ( x, y );

  return value;
}

libMesh::Number
exact_solution( const libMesh::Point& p,
		const libMesh::Parameters& /*params*/,   // parameters, not needed
		const std::string& /*sys_name*/,  // sys_name, not needed
		const std::string& var )  // unk_name, not needed);
{
  // Is the MASA solution we want still loaded and set to the params we want
  //MASA::masa_display_param<libMesh::Real>();

  const double x = p(0);
  const double y = p(1);

  libMesh::Number f = 0.0;
  // Hardcoded to velocity in input file.
  if( var == "u" ) f = MASA::masa_eval_exact_u<libMesh::Real>  ( x, y );
  if( var == "v" ) f = MASA::masa_eval_exact_v<libMesh::Real>  ( x, y );
  if( var == "p" ) f = MASA::masa_eval_exact_p<libMesh::Real>  ( x, y );
  if( var == "nu" ) f = MASA::masa_eval_exact_nu<libMesh::Real>  ( x, y );

  return f;
}

libMesh::Gradient
exact_derivative( const libMesh::Point& p,
		  const libMesh::Parameters& /*params*/,   // parameters, not needed
		  const std::string& /*sys_name*/,  // sys_name, not needed
		  const std::string& var)  // unk_name, not needed);
{
  const double x = p(0);
  const double y = p(1);

  libMesh::Gradient g;

  // Hardcoded to velocity in input file.
  if( var == "u" )
  {
    g(0) = MASA::masa_eval_grad_u<libMesh::Real>  ( x, y, 0 );
    g(1) = MASA::masa_eval_grad_u<libMesh::Real>  ( x, y, 1 );
  }

  if( var == "v" )
  {
    g(0) = MASA::masa_eval_grad_v<libMesh::Real>  ( x, y, 0 );
    g(1) = MASA::masa_eval_grad_v<libMesh::Real>  ( x, y, 1 );
  }

  if( var == "p" )
  {
    g(0) = MASA::masa_eval_grad_p<libMesh::Real>  ( x, y, 0 );
    g(1) = MASA::masa_eval_grad_p<libMesh::Real>  ( x, y, 1 );
  }

  // MASA does not compute gradients of nu at this point, but we need to add it in
  // if( var == "nu" )
  // {
  //   g(0) = MASA::masa_eval_grad_nu<libMesh::Real>  ( x, y, 0 );
  //   g(1) = MASA::masa_eval_grad_nu<libMesh::Real>  ( x, y, 1 );
  // }

  return g;
}

#endif // GRINS_HAVE_MASA
