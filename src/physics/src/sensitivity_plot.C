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
#include "grins/sensitivity_plot.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/assembly_context.h"
#include "grins/generic_ic_handler.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{

  SensitivityPlot::IncompressibleNavierStokesBase(const std::string& physics_name, const GetPot& input )
    : IncompressibleNavierStokesBase(physics_name,
                                         incompressible_navier_stokes, /* "core" Physics name */
                                         input)
  {
    this->read_input_options(input);

    return;
  }

  SensitivityPlot::~Sensitivity()
  {
    return;
  }

  void SensitivityPlot::read_input_options( const GetPot& input )
  {
    return;
  }

  void SensitivityPlot::init_context( AssemblyContext& context );
  {
    // Get a reference to the MultiphysicsSystem using the context
    const MultiphysicsSystem& mphysics_sys = dynamic_cast<const MultiphysicsSystem&>(context.get_system());

    // Get a reference the INSBase physics which the MultiphysicsSystem owns
    const IncompressibleNavierStokesBase<Mu>& ins_physics = dynamic_cast<const IncompressibleNavierStokesBase<Mu>& >(*mphysics_sys.get_physics(incompressible_navier_stokes));

    // Get a pointer to the adjoint solution belonging to the ins_physics
    NumericVector<Number> &adjoint_solution =
        const_cast<System &>(ins_physics).get_adjoint_solution(0);

    // Add this adjoint solution to the vectors that diff context should localize
    context.add_localized_vector(adjoint_solution, ins_physics);

    // Call the base class init_context
    IncompressibleNavierStokesBase::int_context();
  }

  void SensitivityPlot::element_time_derivative( bool compute_jacobian,
						 AssemblyContext& context,
						 CachedValues& /*cache*/ )
  {
    // Get a reference to the MultiphysicsSystem using the context
    const MultiphysicsSystem& mphysics_sys = dynamic_cast<const MultiphysicsSystem&>(context.get_system());

    // Get a reference the INSBase physics which the MultiphysicsSystem owns
    const IncompressibleNavierStokesBase<Mu>& ins_physics = dynamic_cast<const IncompressibleNavierStokesBase<Mu>& >(*mphysics_sys.get_physics(incompressible_navier_stokes));

    // Get a pointer to the adjoint solution belonging to the ins_physics
    NumericVector<Number> &adjoint_solution =
        const_cast<System &>(ins_physics).get_adjoint_solution(0);

     // The number of local degrees of freedom in each variable.
    const unsigned int n_S_dofs = context.get_dof_indices(this->_flow_vars.u_var()).size();

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(this->_flow_vars.u_var())->get_JxW();

    // The Sensitivity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& S_phi =
      context.get_element_fe(this->_flow_vars.u_var())->get_phi();

    const std::vector<libMesh::Point>& S_qpoint =
      context.get_element_fe(this->_flow_vars.u_var())->get_xyz();

    libMesh::DenseSubMatrix<libMesh::Number> &KSS = context.get_elem_jacobian(this->_flow_vars.u_var(), this->_flow_vars.u_var()); // R_{S},{S}

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // Get the previous adjoint solution values at all the qps

    std::vector<libMesh::Gradient> grad_u_adjoint (n_qpoints, 0);

    context.interior_gradient<libMesh::Gradient>(0, adjoint_solution, grad_u_adjoint);

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // Compute the solution & its gradient at the old Newton iterate.
	libMesh::Number S;

	S = context.interior_value(this->_sensitivity_vars.S_var(), qp);

	libMesh::Number p, u, v;
        p = context.interior_value(this->_flow_vars.p_var(), qp);
        u = context.interior_value(this->_flow_vars.u_var(), qp);
        v = context.interior_value(this->_flow_vars.v_var(), qp);

	libMesh::Gradient grad_u, grad_v, grad_w;
        grad_u = context.interior_gradient(this->_flow_vars.u_var(), qp);
        grad_v = context.interior_gradient(this->_flow_vars.v_var(), qp);

	libMesh::Real jac = JxW[qp];

	for (unsigned int i=0; i != n_u_dofs; i++)
          {
            Fu(i) += jac *
	      (grad_u_adjoint[qp](0)*S_phi[i][qp] - S*S_phi[i][qp]);

	    if (compute_jacobian)
              {
		libmesh_not_implemented();
	      }
	  } // End dof loop
      } // End qp loop

    return;
  }
