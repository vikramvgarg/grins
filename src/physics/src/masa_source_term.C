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

// This class
#include "grins/masa_source_term.h"

// GRINS
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{
  MasaSourceTerm::MasaSourceTerm( const std::string& physics_name, const GetPot& input )
    : SourceTermBase(physics_name,input),
      _flow_vars(input,incompressible_navier_stokes),
      _turbulence_vars(input, spalart_allmaras),
      solution_name(input("Physics/"+masa_source_term+"/solution_name", "fans_sa_transient_free_shear"))
  {
    // initialize the problem with the solution the user asked for
    MASA::masa_init<libMesh::Real>("sa_example",solution_name);

    // call the sanity check routine
    // (tests that all variables have been initialized)
    MASA::masa_sanity_check<libMesh::Real>();

    // Set parameters
    MASA::masa_set_param<libMesh::Real>("mu", 2.0);
    MASA::masa_set_param<libMesh::Real>("u_0", 1.0);
    MASA::masa_set_param<libMesh::Real>("u_x", 0.0);
    MASA::masa_set_param<libMesh::Real>("u_y", 0.2);
    MASA::masa_set_param<libMesh::Real>("v_0", 1.0);
    MASA::masa_set_param<libMesh::Real>("v_x", 0.2);
    MASA::masa_set_param<libMesh::Real>("v_y", 0.0);
    MASA::masa_set_param<libMesh::Real>("rho_x", 0.0);
    MASA::masa_set_param<libMesh::Real>("rho_y", 0.0);
    MASA::masa_set_param<libMesh::Real>("p_0", 1.0);
    MASA::masa_set_param<libMesh::Real>("p_x", 0.0);
    MASA::masa_set_param<libMesh::Real>("p_y", 0.0);
    MASA::masa_set_param<libMesh::Real>("nu_sa_x", 0.1);
    MASA::masa_set_param<libMesh::Real>("nu_sa_y", 0.0);
    MASA::masa_set_param<libMesh::Real>("nu_sa_0", 1.0);

    MASA::masa_display_param<libMesh::Real>();
  }

  MasaSourceTerm::~MasaSourceTerm()
  {
    return;
  }

  void MasaSourceTerm::element_time_derivative( bool /*compute_jacobian*/,
                                                  AssemblyContext& context,
                                                  CachedValues& /*cache*/ )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(this->_flow_vars.u_var()).size();
    const unsigned int n_nu_dofs = context.get_dof_indices(this->_turbulence_vars.nu_var()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u_var())->get_JxW();

    // The velocity and turbulent viscosity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& phi_u =
      context.get_element_fe(this->_flow_vars.u_var())->get_phi();
    const std::vector<std::vector<libMesh::Real> >& phi_nu =
      context.get_element_fe(this->_turbulence_vars.nu_var())->get_phi();

    const std::vector<libMesh::Point>& qp_loc = context.get_element_fe(this->_flow_vars.u_var())->get_xyz();

    // Getting time, but we are using steady solutions for now
    libMesh::Real t = context.get_time();

    // Only have 2d turb solutions
    libmesh_assert(context.get_dim()==2);

    // Get residuals
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u_var()); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v_var()); // R_{v}
    libMesh::DenseSubVector<libMesh::Number> &Fnu = context.get_elem_residual(this->_turbulence_vars.nu_var()); // R_{nu}

    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {

      for (unsigned int i=0; i != n_u_dofs; i++)
      {
	// third arg here is time, but we're using as steady solution is
	// steady so t is irrelevant
	Fu(i) += (MASA::masa_eval_source_rho_u<libMesh::Real>  ( qp_loc[qp](0), qp_loc[qp](1), 0.0 ) )*phi_u[i][qp]*JxW[qp];
	Fv(i) += (MASA::masa_eval_source_rho_v<libMesh::Real>  ( qp_loc[qp](0), qp_loc[qp](1), 0.0 ) )*phi_u[i][qp]*JxW[qp];

	//std::cout<<"Fu at ("<<qp_loc[qp](0)<<","<<qp_loc[qp](1)<<") is: "<<Fu(i)<<std::endl;
	//std::cout<<"Fv at ("<<qp_loc[qp](0)<<","<<qp_loc[qp](1)<<") is: "<<Fv(i)<<std::endl;
      }

      for (unsigned int i=0; i != n_nu_dofs; i++)
      {
	// third arg here is time, but we're using as steady solution is
	// steady so t is irrelevant
	Fnu(i) += (MASA::masa_eval_source_nu<libMesh::Real>  ( qp_loc[qp](0), qp_loc[qp](1), 0.0 ) )*phi_nu[i][qp]*JxW[qp];

	//std::cout<<"Fnu at ("<<qp_loc[qp](0)<<","<<qp_loc[qp](1)<<") is: "<<Fnu(i)<<std::endl;
      }
    }

    return;
  }

  // void MasaSourceTerm::element_constraint( bool compute_jacobian,
  // 					   AssemblyContext& context,
  // 					   CachedValues& /*cache*/ )
  // {
  //   // The number of local degrees of freedom in each variable.
  //   const unsigned int n_p_dofs = context.get_dof_indices(this->_flow_vars.p_var()).size();

  //   // We get some references to cell-specific data that
  //   // will be used to assemble the linear system.

  //   // Element Jacobian * quadrature weights for interior integration.
  //   const std::vector<libMesh::Real> &JxW =
  //     context.get_element_fe(this->_flow_vars.u_var())->get_JxW();

  //   // The pressure shape functions at interior quadrature points.
  //   const std::vector<std::vector<libMesh::Real> >& phi_p =
  //     context.get_element_fe(this->_flow_vars.p_var())->get_phi();

  //   const std::vector<libMesh::Point>& qp_loc =
  //     context.get_element_fe(this->_flow_vars.u_var())->get_xyz();

  //   // Getting time, but we are using steady solutions for now
  //   libMesh::Real t = context.get_time();

  //   // Only have 2d turb solutions
  //   libmesh_assert(context.get_dim()==2);

  //   libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(this->_flow_vars.p_var()); // R_{p}

  //   // Add the constraint given by the continuity equation.
  //   unsigned int n_qpoints = context.get_element_qrule().n_points();
  //   for (unsigned int qp=0; qp != n_qpoints; qp++)
  //     {
  //       // Now a loop over the pressure degrees of freedom.  This
  //       // computes the contributions of the continuity equation.
  //       for (unsigned int i=0; i != n_p_dofs; i++)
  //         {
  //           Fp(i) += (MASA::masa_eval_source_rho<libMesh::Real>  ( qp_loc[qp](0), qp_loc[qp](1), 0.0 ) )*phi_p[i][qp]*JxW[qp];

  // 	  } // end loop over p dofs


  //     } // end of the quadrature point (qp) loop

  //   return;
  // }

} // end namespace GRINS

#endif // GRINS_HAVE_MASA
