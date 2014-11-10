//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "grins/full_electroosmosis.h"
#include "grins/assembly_context.h"
#include "grins/physics.h"

// libMesh
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"
#include "libmesh/getpot.h"

namespace GRINS
{

  FullElectroosmosis::FullElectroosmosis( const std::string& physics_name, const GetPot& input )
    : Physics(physics_name, input)
  {
    this->_Debye_Huckel = input("Physics/FullElectroosmosis/Debye_Huckel", 1.0); 
    this->_penalty = input("Physics/FullElectroosmosis/penalty", 1.e+8); 

    this->_E_FE_family = libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+full_electroosmosis+"/E_FE_family", "LAGRANGE") );
   
    this->_T_FE_family = libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+full_electroosmosis+"/T_FE_family", "LAGRANGE") );

    this->_V_FE_family = libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+full_electroosmosis+"/V_FE_family", "LAGRANGE") );

    this->_P_FE_family = libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/"+full_electroosmosis+"/P_FE_family", "LAGRANGE") );

    
    this->_E_order = libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+full_electroosmosis+"/E_order", "SECOND") );

    this->_T_order = libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+full_electroosmosis+"/T_order", "SECOND") );

    this->_V_order = libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+full_electroosmosis+"/V_order", "SECOND") );

    this->_P_order = libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/"+full_electroosmosis+"/P_order", "SECOND") );

    typedef std::set<libMesh::subdomain_id_type>::iterator _enabled_subdomains_iterator;
    
    _enabled_subdomains_iterator _enabled_subdomains_it = _enabled_subdomains.begin();

    const _enabled_subdomains_iterator _enabled_subdomains_end = _enabled_subdomains.end();

    for(; _enabled_subdomains_it !=  _enabled_subdomains.end() ;++_enabled_subdomains_it)
      {
	std::cout<<"Full Electrosomosis enable on subdomain: "<<*_enabled_subdomains_it;
      }

 
    return;
  }

  FullElectroosmosis::~FullElectroosmosis()
  {
    return;
  }

  void FullElectroosmosis::init_variables( libMesh::FEMSystem* system )
  {
    // Create a set of active subdomains
    std::set<subdomain_id_type> active_subdomains;
    // Add the subdomain id for the Full EOF model space
    active_subdomains.insert(0);

    // Get libMesh to assign an index for each variable    
    // E is only active on subdomain 0
    _E_var = system->add_variable( "E", this->_E_order, _E_FE_family,&active_subdomains);

    // T, u, v, p are active on both subdomain 0 and 1
    active_subdomains.insert(1);
    _T_var = system->add_variable( "T", this->_T_order, _T_FE_family,&active_subdomains);

    _u_var = system->add_variable( "u", this->_V_order, _V_FE_family,&active_subdomains);

    _v_var = system->add_variable( "v", this->_V_order, _V_FE_family,&active_subdomains);

    _p_var = system->add_variable( "p", this->_P_order, _P_FE_family,&active_subdomains);
    
    return;
  }

  void FullElectroosmosis::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march temperature forward in time
    system->time_evolving(_E_var);
    system->time_evolving(_T_var);
    system->time_evolving(_u_var);
    system->time_evolving(_v_var);
    system->time_evolving(_p_var);

    return;
  }

  void FullElectroosmosis::init_context( AssemblyContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.get_element_fe(_E_var)->get_JxW();
    context.get_element_fe(_E_var)->get_phi();
    context.get_element_fe(_E_var)->get_dphi();    

    context.get_side_fe(_E_var)->get_JxW();
    context.get_side_fe(_E_var)->get_phi();
    context.get_side_fe(_E_var)->get_dphi();    

    context.get_element_fe(_T_var)->get_JxW();
    context.get_element_fe(_T_var)->get_phi();
    context.get_element_fe(_T_var)->get_dphi();
    context.get_element_fe(_T_var)->get_xyz();

    context.get_side_fe(_T_var)->get_JxW();
    context.get_side_fe(_T_var)->get_phi();
    context.get_side_fe(_T_var)->get_dphi();
    context.get_side_fe(_T_var)->get_xyz();

    context.get_element_fe(_u_var)->get_JxW();
    context.get_element_fe(_u_var)->get_phi();
    context.get_element_fe(_u_var)->get_dphi();    

    context.get_side_fe(_u_var)->get_JxW();
    context.get_side_fe(_u_var)->get_phi();
    context.get_side_fe(_u_var)->get_dphi();    

    context.get_element_fe(_p_var)->get_JxW();
    context.get_element_fe(_p_var)->get_phi();
    context.get_element_fe(_p_var)->get_dphi();    

    context.get_side_fe(_p_var)->get_JxW();
    context.get_side_fe(_p_var)->get_phi();
    context.get_side_fe(_p_var)->get_dphi();

    return;
  }

  void FullElectroosmosis::element_time_derivative( bool compute_jacobian,
						AssemblyContext& context,
						CachedValues& /*cache*/ )
  {
    // The number of local degrees of freedom in each variable.
    //const unsigned int n_E_dofs = context.dof_indices_var(_E_var).size();
    const unsigned int n_T_dofs = context.get_dof_indices(_T_var).size();
    const unsigned int n_u_dofs = context.get_dof_indices(_u_var).size();
    const unsigned int n_p_dofs = context.get_dof_indices(_p_var).size();

    // Check number of dofs is same for _u_var, _v_var and _T_var.
    //libmesh_assert (n_u_dofs == context.get_dof_indices(_E_var).size());    
    libmesh_assert (n_u_dofs == context.get_dof_indices(_v_var).size());
    libmesh_assert (n_u_dofs == context.get_dof_indices(_T_var).size());    

    // We get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(_T_var)->get_JxW();

    // The temperature shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(_T_var)->get_phi();
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(_T_var)->get_dphi();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.get_element_fe(_p_var)->get_phi();

    const std::vector<libMesh::Point>& q_points = 
      context.get_element_fe(_T_var)->get_xyz();

    // The subvectors and submatrices we need to fill:
    //
    // K_{\alpha \beta} = R_{\alpha},{\beta} = \partial{ R_{\alpha} } / \partial{ {\beta} } (where R denotes residual)
    // e.g., for \alpha = T and \beta = v we get: K_{Tu} = R_{T},{u}
    //

    libMesh::DenseSubMatrix<Number> &KEE = context.get_elem_jacobian(_E_var, _E_var); // R_{E},{E}

    libMesh::DenseSubMatrix<Number> &KTT = context.get_elem_jacobian(_T_var, _T_var); // R_{T},{T}
    
    libMesh::DenseSubMatrix<libMesh::Number> &KuT = context.get_elem_jacobian(_u_var, _T_var); // R_{u},{V}
    libMesh::DenseSubMatrix<libMesh::Number> &Kuu = context.get_elem_jacobian(_u_var, _u_var); // R_{u},{u}
    libMesh::DenseSubMatrix<libMesh::Number> &Kup = context.get_elem_jacobian(_u_var, _p_var); // R_{u},{p}

    libMesh::DenseSubMatrix<libMesh::Number> &KvT = context.get_elem_jacobian(_v_var, _T_var); // R_{v},{V}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvv = context.get_elem_jacobian(_v_var, _v_var); // R_{v},{v}
    libMesh::DenseSubMatrix<libMesh::Number> &Kvp = context.get_elem_jacobian(_v_var, _p_var); // R_{v},{p}
    
    libMesh::DenseSubVector<Number> &FE = context.get_elem_residual(_E_var); // R_{E}
    libMesh::DenseSubVector<Number> &FT = context.get_elem_residual(_T_var); // R_{T}
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(_u_var); // R_{u}
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(_v_var); // R_{v}
    
    // Now we will build the element Jacobian and residual.
    // Constructing the residual requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	// Compute the solution & its gradient at the old Newton iterate.
	libMesh::Number E;
	E = context.interior_value(_E_var, qp);

	libMesh::Gradient grad_E;
	grad_E = context.interior_gradient(_E_var, qp);

	libMesh::Number T;
	T = context.interior_value(_T_var, qp);

	libMesh::Gradient grad_T;
	grad_T = context.interior_gradient(_T_var, qp);

	libMesh::Number u, v, p;
	u = context.interior_value(_u_var, qp);
	v = context.interior_value(_v_var, qp);
	p = context.interior_value(_p_var, qp);
	
	libMesh::NumberVectorValue U(u, v);

	libMesh::Gradient grad_u, grad_v;
	grad_u = context.interior_gradient(_u_var, qp);
	grad_v = context.interior_gradient(_v_var, qp);
		
	libMesh::Number T_x, T_y, u_x, u_y, v_x, v_y;
	T_x = grad_T(0);
	T_y = grad_T(1);

	u_x = grad_u(0);
	u_y = grad_u(1);

	v_x = grad_v(0);
	v_y = grad_v(1);

	// First, an i-loop over the  degrees of freedom.
	for (unsigned int i=0; i != n_T_dofs; i++)
	  {
	    FE(i) += JxW[qp] *
	      (-T_gradphi[i][qp]*grad_E // diffusion term
	       + (pow(_Debye_Huckel, 2.)*sinh(E)*T_phi[i][qp])); // nonlinear term

	    FT(i) += JxW[qp] *
	      (-T_gradphi[i][qp]*grad_T);  // diffusion term
	    
	    Fu(i) += JxW[qp] * 
	      (p*T_gradphi[i][qp](0)               // pressure term
	       - (grad_u*T_gradphi[i][qp])  // diffusion term    
	       - (sinh(E)*T_x*T_phi[i][qp]) );  // Nonlinear Transport term
	    
	    Fv(i) += JxW[qp]*
	      (p*T_gradphi[i][qp](1)                // pressure term
	       - (grad_v*T_gradphi[i][qp]) // diffusion term
	       - (sinh(E)*T_y*T_phi[i][qp]) );  // Nonlinear Transport term

	    if (compute_jacobian)
	      {
		for (unsigned int j=0; j != n_T_dofs; j++)
		  {		    
		    KEE(i,j) += JxW[qp] *
		      ( T_gradphi[i][qp]*T_gradphi[j][qp] // diffusion term
			- (pow(_Debye_Huckel, 2.)*cosh(E)*T_phi[i][qp]*T_phi[j][qp])); // nonlinear term

		    KTT(i,j) += JxW[qp] *
		      ( T_gradphi[i][qp]*T_gradphi[j][qp]); // diffusion term
		    
		    Kuu(i,j) += JxW[qp]*
		    ( -(T_gradphi[i][qp]*T_gradphi[j][qp]) // diffusion term 
		      - (cosh(E)*T_x*T_phi[i][qp]*T_phi[j][qp]) // nonlinear transport term
		      - (sinh(E)*T_gradphi[j][qp](0)*T_phi[i][qp])); // nonlinear transport term

		    Kvv(i,j) += JxW[qp]*
		      (-(T_gradphi[i][qp]*T_gradphi[j][qp])  // diffusion term
		       - (cosh(E)*T_y*T_phi[i][qp]*T_phi[j][qp]) // nonlinear transport term
		       - (sinh(E)*T_gradphi[j][qp](1)*T_phi[i][qp])); // nonlinear transport term
		  }

	      for (unsigned int j=0; j<n_p_dofs; j++)
		{
		  Kup(i,j) += JxW[qp]*(p_phi[j][qp]*T_gradphi[i][qp](0));
		  Kvp(i,j) += JxW[qp]*(p_phi[j][qp]*T_gradphi[i][qp](1));

		} // end of the inner dof (j) loop

	      } // end - if (compute_jacobian && context.elem_solution_derivative)

	  } // end of the outer dof (i) loop
      } // end of the quadrature point (qp) loop
  
    return;
  }

    void FullElectroosmosis::element_constraint( bool compute_jacobian,
				   AssemblyContext& context,
				   CachedValues& /*cache*/ )
  {
    // The number of local degrees of freedom in each variable.
    const unsigned int n_u_dofs = context.get_dof_indices(_u_var).size();
    const unsigned int n_p_dofs = context.get_dof_indices(_p_var).size();

    // We get some references to cell-specific data that
    // will be used to assemble the linear system.

    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(_u_var)->get_JxW();

    // The velocity shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(_u_var)->get_dphi();

    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.get_element_fe(_p_var)->get_phi();

    // The subvectors and submatrices we need to fill:
    //
    // Kpu, Kpv, Kpw, Fp
        
    libMesh::DenseSubMatrix<libMesh::Number> &Kpu = context.get_elem_jacobian(_p_var, _u_var); // R_{p},{u}
    libMesh::DenseSubMatrix<libMesh::Number> &Kpv = context.get_elem_jacobian(_p_var, _v_var); // R_{p},{v}
    

    libMesh::DenseSubVector<libMesh::Number> &Fp = context.get_elem_residual(_p_var); // R_{p}

    // Add the constraint given by the continuity equation.
    unsigned int n_qpoints = context.get_element_qrule().n_points();
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
	// Compute the velocity gradient at the old Newton iterate.
	libMesh::Gradient grad_u, grad_v, grad_w;
	grad_u = context.interior_gradient(_u_var, qp);
	grad_v = context.interior_gradient(_v_var, qp);	

	// Now a loop over the pressure degrees of freedom.  This
	// computes the contributions of the continuity equation.
	for (unsigned int i=0; i != n_p_dofs; i++)
	  {
	    Fp(i) += JxW[qp] * p_phi[i][qp] *
	      (grad_u(0) + grad_v(1));	    

	    if (compute_jacobian)
	      {
		for (unsigned int j=0; j != n_u_dofs; j++)
		  {
		    Kpu(i,j) += JxW[qp]*p_phi[i][qp]*u_gradphi[j][qp](0);
		    Kpv(i,j) += JxW[qp]*p_phi[i][qp]*u_gradphi[j][qp](1);		    
		  } // end of the inner dof (j) loop

	      } // end - if (compute_jacobian && context.elem_solution_derivative)

	  } // end of the outer dof (i) loop
      } // end of the quadrature point (qp) loop

    return;
  }

void FullElectroosmosis::side_constraint (bool compute_jacobian,
					     AssemblyContext& context, CachedValues& /*cache*/)
{ 
  // The number of local degrees of freedom in each variable
  //const unsigned int n_E_dofs = context.get_dof_indices(_E_var).size();
  const unsigned int n_T_dofs = context.get_dof_indices(_T_var).size();
  const unsigned int n_u_dofs = context.get_dof_indices(_u_var).size();
  const unsigned int n_p_dofs = context.get_dof_indices(_p_var).size();
  
  // Check number of dofs is same for _u_var, _v_var and _T_var.
  //libmesh_assert (n_u_dofs == context.get_dof_indices(_E_var).size());
  libmesh_assert (n_u_dofs == context.get_dof_indices(_v_var).size());
  libmesh_assert (n_u_dofs == context.get_dof_indices(_T_var).size());    

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weights for interior integration

    const std::vector<Real> &JxW_side = context.get_side_fe(_T_var)->get_JxW();

  const std::vector<std::vector<Real> > &phi_side = context.get_side_fe(_T_var)->get_phi();

  // The velocity shape function gradients at boundary quadrature points, needed for boundary condition coupling

  const std::vector<std::vector<RealGradient> > &dphi_side =
    context.get_side_fe(_T_var)->get_dphi();
  
  const std::vector<Point > &qside_point = context.get_side_fe(_T_var)->get_xyz();

  const std::vector<Point> &face_normals = context.get_side_fe(_T_var)->get_normals();
  
  // The subvectors and submatrices we need to fill:
  DenseSubMatrix<Number> &KEE = context.get_elem_jacobian(_E_var,_E_var);

  DenseSubMatrix<Number> &KTT = context.get_elem_jacobian(_T_var,_T_var);
    
  DenseSubMatrix<Number> &Kuu = context.get_elem_jacobian(_u_var,_u_var);
  
  DenseSubMatrix<Number> &Kvv = context.get_elem_jacobian(_v_var,_v_var);  
  
  DenseSubVector<Number> &FE = context.get_elem_residual(_E_var);
  DenseSubVector<Number> &FT = context.get_elem_residual(_T_var);
  DenseSubVector<Number> &Fu = context.get_elem_residual(_u_var);
  DenseSubVector<Number> &Fv = context.get_elem_residual(_v_var);

  unsigned int n_qpoints = context.get_side_qrule().n_points();
  
  Real TOL = 1.e-10;

  // Get all the boundary ids associated with this side
  std::vector<BoundaryID> ids = context.side_boundary_ids();

  for( std::vector<BoundaryID>::const_iterator it = ids.begin();
	 it != ids.end(); it++ )
    {
      libmesh_assert (*it != libMesh::BoundaryInfo::invalid_id);

      // Get the bdry id for this bdry
      short int bc_id = *it;
      
      libMesh::Real pot_set_0 = 0.;
      libMesh::Real pot_set_1 = 1.;
      libMesh::Real pot_set_2 = 1.;
      
      libMesh::Real E_pot = 1.;

      for (unsigned int qp=0; qp != n_qpoints; qp++)
	{
	  // Compute the solution at the old Newton iterate
	  libMesh::Number E = context.side_value(_E_var, qp);

	  libMesh::Gradient grad_T = context.side_gradient(_T_var, qp);
	  
	  libMesh::Number T = context.side_value(_T_var, qp), u = context.side_value(_u_var, qp), v = context.side_value(_v_var, qp), p = context.side_value(_p_var, qp);
	  
	  // Variable to hold Dirchlet data for potential
	  
	  Real T_dirichlet ;

      // x and y co-ordinates of the current qp
      Real x = qside_point[qp](0);
      Real y = qside_point[qp](1);

      switch (bc_id)
	{

	case 0 : // Station 0
	  {	    	   
	    T_dirichlet = pot_set_0;
	    
	    // Penalize potential
	    for (unsigned int i=0; i != n_u_dofs; i++)
	      {
		FT(i) += JxW_side[qp] * _penalty * ( T - T_dirichlet) * phi_side[i][qp];
		
		if (compute_jacobian)
		  for (unsigned int j=0; j != n_u_dofs; ++j)
		    KTT(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp];
	      }

	    for (unsigned int i=0; i != n_u_dofs; i++)
	      {
	    	Fu(i) += JxW_side[qp] * _penalty * (u - 0.0) * phi_side[i][qp] ; // Enforces Dirichlet Condition for p
	    	
	    	if (compute_jacobian)
	    	  for (unsigned int j=0; j != n_u_dofs; ++j)
	    	    Kuu(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp];		      		      		      	    	      		      		      	    	    
	      }
	    	    
	    break;	    
	  }

	case 1 : // Station 1
	  {
	    T_dirichlet = pot_set_1;
	    
	    // Penalize potential
	    for (unsigned int i=0; i != n_u_dofs; i++)
	      {
		FT(i) += JxW_side[qp] * _penalty * ( T - T_dirichlet) * phi_side[i][qp];
		
		if (compute_jacobian)
		  for (unsigned int j=0; j != n_u_dofs; ++j)
		    KTT(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp];
	      }

	    for (unsigned int i=0; i != n_u_dofs; i++)
	      {
	    	Fv(i) += JxW_side[qp] * _penalty * (v - 0.0) * phi_side[i][qp] ; // Enforces Dirichlet Condition for p
	    	
	    	if (compute_jacobian)
	    	  for (unsigned int j=0; j != n_u_dofs; ++j)
	    	    Kvv(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp];		      		      		      	    	      		      		      	    	    
	      }

	    break;

	  }

	case 2: // Station 2
	  {
	    T_dirichlet = pot_set_2;
	    
	      // Penalize potential
	      for (unsigned int i=0; i != n_u_dofs; i++)
		{
		  FT(i) += JxW_side[qp] * _penalty * ( T - T_dirichlet) * phi_side[i][qp];		  
		  
		  if (compute_jacobian)		    
		  for (unsigned int j=0; j != n_u_dofs; ++j)
		    KTT(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp];
		}
	      
	      for (unsigned int i=0; i != n_u_dofs; i++)
		{
		  Fv(i) += JxW_side[qp] * _penalty * (v - 0.0) * phi_side[i][qp] ; // Enforces Dirichlet Condition for p
		  
		  if (compute_jacobian)
		    for (unsigned int j=0; j != n_u_dofs; ++j)
		      Kvv(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp];		      		      		      	    	      		      		      	    	    
		}
	      
	      break;
	    }

	case 3: // Wall Section 1
	  {	    	   
	    // Penalize non-slip velocities
	    for (unsigned int i=0; i != n_u_dofs; i++)
	      {
		FE(i) += JxW_side[qp] * _penalty * ( E - E_pot ) * phi_side[i][qp]; // Enforces No Penetration for u						 		    
		Fu(i) += JxW_side[qp] * _penalty * ( u - 0.0 ) * phi_side[i][qp]; // Enforces No Penetration for u									    		
		Fv(i) += JxW_side[qp] * _penalty * ( v - 0.0 ) * phi_side[i][qp]; // Enforces No Penetration for u									    		

		if (compute_jacobian)
		  for (unsigned int j=0; j != n_u_dofs; ++j)
		    {
		      KEE(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp]; // Enforces No Penetration for u		      

		      Kuu(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp]; // Enforces No Penetration for u		      
		      		      		      		      
		      Kvv(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp]; // Enforces Slip Condition for v		      
		    }
						    
	      }
	    break;
	  }

	  case 4: // Wall Section 2
	  {	    	   
	    // Penalize non-slip velocities
	    for (unsigned int i=0; i != n_u_dofs; i++)
	      {
		FE(i) += JxW_side[qp] * _penalty * ( E - E_pot ) * phi_side[i][qp]; // Enforces No Penetration for u						 		    
		Fu(i) += JxW_side[qp] * _penalty * ( u - 0.0 ) * phi_side[i][qp]; // Enforces Slip Condition for u

		Fv(i) += JxW_side[qp] * _penalty * ( v - 0.0 ) * phi_side[i][qp]; // Enforces No Penetration for v
							    		
		if (compute_jacobian)
		  for (unsigned int j=0; j != n_u_dofs; ++j)
		    {
		      KEE(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp]; // Enforces No Penetration for u
		      
		      Kuu(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp]; // Enforces Slip Condition for u		      
		      		      		      
		      Kvv(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp]; // Enforces No Penetration for v		      
		    }
						    
	      }
	    break;
	  }

	  case 5: // Wall Section 3, note the change of sign on the _penalty
	  {	    	   
	    // Penalize non-slip velocities
	    for (unsigned int i=0; i != n_u_dofs; i++)
	      {
		FE(i) += JxW_side[qp] * _penalty * ( E - E_pot ) * phi_side[i][qp]; // Enforces No Penetration for u						 		    
		Fu(i) += JxW_side[qp] * _penalty * (u - 0.0 ) * phi_side[i][qp]; // Enforces Slip Condition for u

		Fv(i) += JxW_side[qp] * _penalty * ( v - 0.0 ) * phi_side[i][qp]; // Enforces No Penetration for v
							    		
		if (compute_jacobian)
		  for (unsigned int j=0; j != n_u_dofs; ++j)
		    {
		      KEE(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp]; // Enforces No Penetration for u

		      Kuu(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp]; // Enforces Slip Condition for u		      
		      		      		      
		      Kvv(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp]; // Enforces No Penetration for v		      
		    }
						    
	      }
	    break;
	  }

	case 6: // Wall Section 4
	  {	    	   
	    // Penalize non-slip velocities
	    for (unsigned int i=0; i != n_u_dofs; i++)
	      {
		FE(i) += JxW_side[qp] * _penalty * ( E - E_pot ) * phi_side[i][qp]; // Enforces No Penetration for u						 		    
		Fu(i) += JxW_side[qp] * _penalty * ( u - 0.0 ) * phi_side[i][qp]; // Enforces Slip Condition for u 

		Fv(i) += JxW_side[qp] * _penalty * ( v  - 0.0) * phi_side[i][qp]; // Enforces No Penetration for v
							    		
		if (compute_jacobian)
		  for (unsigned int j=0; j != n_u_dofs; ++j)
		    {
		      KEE(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp]; // Enforces No Penetration for u

		      Kuu(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp]; // Enforces Slip Condition for u 		      
		      		      		      
		      Kvv(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp]; // Enforces No Penetration for v
		    }
						    
	      }
	    break;
	  }

	  case 7: // Wall Section 5, note the change of sign on the _penalty
	  {	    	   
	    // Penalize non-slip velocities
	    for (unsigned int i=0; i != n_u_dofs; i++)
	      {
		FE(i) += JxW_side[qp] * _penalty * ( E - E_pot ) * phi_side[i][qp]; // Enforces No Penetration for u						 		    
		Fu(i) += JxW_side[qp] * _penalty * ( u - 0.0 ) * phi_side[i][qp]; // Enforces No Penetration for u

		Fv(i) += JxW_side[qp] * _penalty * ( v - 0.0 ) * phi_side[i][qp]; // Enforces Slip Condition for v
							    		
		if (compute_jacobian)
		  for (unsigned int j=0; j != n_u_dofs; ++j)
		    {
		      KEE(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp]; // Enforces No Penetration for u

		      Kuu(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp]; // Enforces No Penetration for u		      
		      		      		      		      
		      Kvv(i,j) += JxW_side[qp] * _penalty * phi_side[j][qp] * phi_side[i][qp]; // Enforces Slip Condition for v		      
		    }
						    
	      }
	    break;
	  }
	
	default :
	    {
	      std::cout<<"Velocity Dirchlet boundary condition not set ! Boundary id: "<<bc_id<<std::endl<<std::endl;
	      libmesh_error();
	    }
	    
	
	} // end switch    
	
    } // end of the quadrature point qp-loop
      
    } // Loop over all boundary ids associated with this side 
 }


  void FullElectroosmosis::mass_residual( bool compute_jacobian,
				      AssemblyContext& context,
				      CachedValues& /*cache*/ )
  {
    libmesh_not_implemented();

    // // First we get some references to cell-specific data that
    // // will be used to assemble the linear system.

    // // Element Jacobian * quadrature weights for interior integration
    // const std::vector<Real> &JxW = 
    //   context.get_element_fe(_T_var)->get_JxW();

    // // The shape functions at interior quadrature points.
    // const std::vector<std::vector<Real> >& phi = 
    //   context.get_element_fe(_T_var)->get_phi();

    // // The number of local degrees of freedom in each variable
    // const unsigned int n_T_dofs = context.get_dof_indices(_T_var).size();

    // // The subvectors and submatrices we need to fill:
    // DenseSubVector<Real> &F = context.get_elem_residual(_T_var);

    // DenseSubMatrix<Real> &M = context.get_elem_jacobian(_T_var)(_T_var);

    // unsigned int n_qpoints = context.element_qrule->n_points();

    // for (unsigned int qp = 0; qp != n_qpoints; ++qp)
    //   {
    // 	// For the mass residual, we need to be a little careful.
    // 	// The time integrator is handling the time-discretization
    // 	// for us so we need to supply M(u_fixed)*u for the residual.
    // 	// u_fixed will be given by the fixed_interior_* functions
    // 	// while u will be given by the interior_* functions.
    // 	Real T_dot = context.interior_value(_T_var, qp);

    // 	for (unsigned int i = 0; i != n_T_dofs; ++i)
    // 	  {
    // 	    F(i) += JxW[qp]*(_rho*_Cp*T_dot*phi[i][qp] );

    // 	    if( compute_jacobian )
    //           {
    //             for (unsigned int j=0; j != n_T_dofs; j++)
    //               {
    // 		    // We're assuming rho, cp are constant w.r.t. T here.
    //                 M(i,j) += JxW[qp]*_rho*_Cp*phi[j][qp]*phi[i][qp] ;
    //               }
    //           }// End of check on Jacobian
          
    // 	  } // End of element dof loop
      
    //   } // End of the quadrature point loop

    // return;
  }

} // namespace GRINS
