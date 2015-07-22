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
#include "grins/drag.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/assembly_context.h"

#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/spalart_allmaras_viscosity.h"
#include "grins/drag_macro.h"

#include "grins/inc_navier_stokes_base.h"
// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/fe_base.h"
#include "libmesh/equation_systems.h"

namespace GRINS
{
  template<class Mu>
  Drag<Mu>::Drag( const std::string& qoi_name , const GetPot& input)
    : QoIBase(qoi_name)
      //_mu(input)
  {
    return;
  }

  template<class Mu>
  Drag<Mu>::~Drag()
  {
    return;
  }

  template<class Mu>
  QoIBase* Drag<Mu>::clone() const
  {
    Drag<Mu>& my_ref = const_cast<Drag<Mu>&>(*this);
    return new Drag ( my_ref );
  }

  template<class Mu>
  void Drag<Mu>::init( const GetPot& input, MultiphysicsSystem& system )
  {
    // Extract subdomain on which to compute to qoi
    //int num_ids = input.vector_variable_size( "QoI/Vorticity/enabled_subdomains" );

    //if( num_ids == 0 )
    //{
    //std::cerr << "Error: Must specify at least one subdomain id on which to compute vorticity." << std::endl;
    //libmesh_error();
    //}

    //for( int i = 0; i < num_ids; i++ )
    //{
    //libMesh::subdomain_id_type s_id = input( "QoI/Vorticity/enabled_subdomains", -1, i );
    //	_subdomain_ids.insert( s_id );
    //}

    // Should we define the adjoint boundary condition here ?

    // Grab velocity variable indices
    std::string u_var_name = input("Physics/VariableNames/u_velocity", u_var_name_default);
    std::string v_var_name = input("Physics/VariableNames/v_velocity", v_var_name_default);
    std::string p_var_name = input("Physics/VariableNames/pressure", p_var_name_default);

    this->_u_var = system.variable_number(u_var_name);
    this->_v_var = system.variable_number(v_var_name);
    this->_p_var = system.variable_number(p_var_name);

    // Initialize the viscosity object
    //this->_mu.init(&(equation_system->get_system(0)));

    return;
  }

  template<class Mu>
  void Drag<Mu>::init_context( AssemblyContext& context )
  {
    libMesh::FEBase* u_fe = NULL;
    libMesh::FEBase* v_fe = NULL;

    context.get_element_fe<libMesh::Real>(this->_u_var, u_fe);
    context.get_element_fe<libMesh::Real>(this->_v_var, v_fe);

    u_fe->get_dphi();
    u_fe->get_JxW();

    v_fe->get_dphi();

    return;
  }

  template<class Mu>
  void Drag<Mu>::element_qoi( AssemblyContext& context,
                               const unsigned int qoi_index )
  {
    //if( _subdomain_ids.find( (&context.get_elem())->subdomain_id() ) != _subdomain_ids.end() )
    //{
    // A reference to the system
    const libMesh::System & sys = context.get_system();

    // Create a boundary info object for checking boundary ids later
    const libMesh::BoundaryInfo& boundary_info = *sys.get_mesh().boundary_info;

    // Get a reference to the current element
    libMesh::Elem &elem = context.get_elem();

    // Check if this element is on the boundary
    if(elem.on_boundary())
      {
	// The number of sides on the element
	unsigned int n_sides = elem.n_sides();

	// Loop over the sides that are on the boundary
	for(unsigned char s=0; s != n_sides; ++s)
	  {
	    //Is one of the sides on the boundary we care about
	    if(boundary_info.has_boundary_id(&elem,s,0))
	      {
		// We now know that we are on boundary element with a boundary
		// that corresponds to the one for which we want the drag
		// and can proceed with obtaining the contribution from this element
		libMesh::FEBase* element_fe;
		context.get_element_fe<libMesh::Real>(this->_u_var, element_fe);
		const std::vector<libMesh::Real> &JxW = element_fe->get_JxW();

		// We also need basis functions and derivatives for this QoI
		const std::vector<std::vector<libMesh::Real>>& u_phi =
		  context.get_element_fe(_u_var)->get_phi();
		const std::vector<std::vector<libMesh::RealGradient> >& du_phi =
		  context.get_element_fe(_u_var)->get_dphi();

		// Local DOF count and quadrature point count
		const unsigned int n_u_dofs = context.get_dof_indices(this->_u_var).size();
		unsigned int n_qpoints = context.get_element_qrule().n_points();

		/*! \todo Need to generalize this to the multiple QoI case */
		libMesh::Number& qoi = context.get_qois()[qoi_index];

		for( unsigned int qp = 0; qp != n_qpoints; qp++ )
		  {
		    libMesh::Real u = 0.;
		    libMesh::Real v = 0.;
		    context.interior_value ( this->_u_var, qp, u);
		    context.interior_value ( this->_v_var, qp, v);

		    libMesh::Gradient grad_u = 0.;
		    libMesh::Gradient grad_v = 0.;
		    context.interior_gradient( this->_u_var, qp, grad_u );
		    context.interior_gradient( this->_v_var, qp, grad_v );

		    libMesh::Real p = 0.;
		    context.interior_value (this->_p_var, qp, p);

		    // Get a reference to the MultiphysicsSystem using the context
		    const MultiphysicsSystem& mphysics_sys = dynamic_cast<const MultiphysicsSystem&>(context.get_system());
		    // Get a reference the INSBase physics which the MultiphysicsSystem owns
		    const IncompressibleNavierStokesBase<Mu>& ins_physics = dynamic_cast<const IncompressibleNavierStokesBase<Mu>& >(*mphysics_sys.get_physics(incompressible_navier_stokes));
		// Use the get_viscosity_value function to get the viscosity at this qp
		libMesh::Real _mu_qp = ins_physics.get_viscosity_value(context, qp);

		    for( unsigned int i = 0; i != n_u_dofs; i++ )
		      {
			// Need to access nu from INS
			qoi += ( _mu_qp*( 2*grad_u(0)*du_phi[i][qp](0) ) + _mu_qp*( (grad_u(1) + grad_v(0))*du_phi[i][qp](1) ) + ( u*grad_u(0)*u_phi[i][qp] + v*grad_u(1)*u_phi[i][qp] ) - ( p*du_phi[i][qp](0) ) ) * JxW[qp];
		      } // End dof (basis function) loop

		  } // End qp loop
	      } // End if has boundary id
	  } // End loop over sides
      } // End if elem on boundary
    //}

    return;
  }

  template<class Mu>
  void Drag<Mu>::element_qoi_derivative( AssemblyContext& context,
                                                  const unsigned int qoi_index )
  {
    // Drag QoI RHS defined completely by adjoint dirichlet bc so no need to do anything
    return;
  }

  // template<class Mu>
  // void Drag<Mu>::register_parameter
  //   ( const std::string & param_name,
  //     libMesh::ParameterMultiPointer<libMesh::Number> & param_pointer )
  //   const
  // {
  //   ParameterUser::register_parameter(param_name, param_pointer);
  //   this->_mu.register_parameter(param_name, param_pointer);
  // }

} //namespace GRINS

// Instantiate
INSTANTIATE_DRAG_SUBCLASS(Drag);
