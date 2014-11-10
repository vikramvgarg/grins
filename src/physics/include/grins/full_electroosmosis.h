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
#ifndef GRINS_FULL_ELECTROOSMOSIS_H
#define GRINS_FULL_ELECTROOSMOSIS_H

#include "grins/physics.h"
#include "grins/inc_navier_stokes_base.h"

namespace GRINS
{

  class FullElectroosmosis : public Physics
  {

  public:

    FullElectroosmosis( const std::string& physics_name, const GetPot& input );

    ~FullElectroosmosis();

        //! Initialize variables for this physics.
    virtual void init_variables( libMesh::FEMSystem* system );

    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext& context );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    //! Time dependent part(s) of physics for element interiors
    virtual void element_time_derivative( bool compute_jacobian,
					  AssemblyContext& context,
					  CachedValues& cache );
    
    // Constraint part(s)
    virtual void element_constraint( bool compute_jacobian,
				     AssemblyContext& context,
				     CachedValues& cache );

    virtual void side_constraint( bool compute_jacobian,
				     AssemblyContext& context,
				     CachedValues& cache );

    // Mass matrix part(s)
    virtual void mass_residual( bool compute_jacobian,
				AssemblyContext& context,
				CachedValues& cache );

    protected:

    // The slip parameter lambda u.t = -lambda grad_T.t 
    libMesh::Real _Debye_Huckel;
    
    // The penalty parameter to enforce boundary conditions
    libMesh::Real _penalty;

    //! Indices for each variable;
    VariableIndex _E_var; /* Index for potential field */
    VariableIndex _T_var; /* Index for potential field */
    VariableIndex _u_var; /* Index for x-velocity field */
    VariableIndex _v_var; /* Index for y-velocity field */
    VariableIndex _p_var; /* Index for pressure field */

    std::string _E_var_name, _T_var_name, _u_var_name, _v_var_name, _p_var_name;

    libMeshEnums::Order _E_order, _T_order, _V_order, _P_order;

    libMeshEnums::FEFamily _E_FE_family, _T_FE_family, _V_FE_family, _P_FE_family;

  private:
    
    FullElectroosmosis();

  };

} // end namespace block

#endif // GRINS_STOKES_H
