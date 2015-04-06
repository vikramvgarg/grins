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
      solution_name(input("Physics/"+spalart_allmaras+"/Parameters/solution_name","fans_sa_transient_shear"))
  {
    // initialize the problem with the solution the user asked for
    err = MASA::masa_init<Scalar>("sa example",solution_name);

    // call the sanity check routine
    // (tests that all variables have been initialized)
    err = MASA::masa_sanity_check<Scalar>();
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

    const std::vector<libMesh::Point>& x_qp = context.get_element_fe(this->_flow_vars.u_var())->get_xyz();

    libMesh::Real t = context.get_time();

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

      for (unsigned int i=0; i != n_dofs; i++)
      {
	F_u(i) += (MASA::masa_eval_source_rho_u<Scalar>  (x_qp[0],x_qp[1],t))*phi_u[i][qp]*JxW[qp];
	F_v(i) += (MASA::masa_eval_source_rho_v<Scalar>  (x_qp[0],x_qp[1],t))*phi_u[i][qp]*JxW[qp];
	F_nu(i) += (MASA::masa_eval_source_nu<Scalar>  (x_qp[0],x_qp[1],t))*phi_nu[i][qp]*JxW[qp];
      }
    }

    return;
  }

} // end namespace GRINS
