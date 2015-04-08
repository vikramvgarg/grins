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

#ifndef GRINS_MASA_SOURCE_TERM_H
#define GRINS_MASA_SOURCE_TERM_H

// GRINS
#include "grins/source_term_base.h"
#include "primitive_flow_fe_variables.h"
#include "grins/turbulence_fe_variables.h"

// MASA
#ifdef GRINS_HAVE_MASA
#include "masa.h"

namespace GRINS
{

  class MasaSourceTerm : public SourceTermBase
  {
  public:

    MasaSourceTerm( const std::string& physics_name, const GetPot& input );

    virtual ~MasaSourceTerm();

    virtual void element_time_derivative( bool compute_jacobian,
					  AssemblyContext& context,
					  CachedValues& cache );

  protected:

    // The flow variables
    PrimitiveFlowFEVariables _flow_vars;

    // These are defined for each physics
    TurbulenceFEVariables _turbulence_vars;

    // A string to hold the particular manufactured solution we want to verify against
    std::string solution_name;

  private:

    MasaSourceTerm();

  };

} // end namespace GRINS

#endif // End if GRINS_HAVE_MASA

#endif // GRINS_MASA_SOURCE_TERM_H
