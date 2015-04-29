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
#include "grins/forcing_function_evaluation.h"

// GRINS
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  ForcingFunctionEvaluation::ForcingFunctionEvaluation (const std::string& physics_name, const GetPot& input) : Physics(physics_name, input)
  {
    //Constructor stuff

    return;
  }

  ForcingFunctionEvaluation::~ForcingFunctionEvaluation()
  {
    return;
  }

  void ForcingFunctionEvaluation::init_variables( libMesh::FEMSystem* system)
  {
    return;
  }

  libMesh::Real ForcingFunctionEvaluation::compute_masa_forcing( libMesh::Real x, libMesh::Real y, std::string& unknown_name)
  {
    libMesh::Real f = 0.0;

    if( unknown_name == "u" ) value = MASA::masa_eval_source_rho_u<libMesh::Real>  ( x, y );
    if( unknown_name == "v" ) value = MASA::masa_eval_source_rho_v<libMesh::Real>  ( x, y );
    if( unknown_name == "p" ) value = MASA::masa_eval_source_rho<libMesh::Real>  ( x, y );
    if( unknown_name == "nu" ) value = MASA::masa_eval_source_nu<libMesh::Real>  ( x, y );

    return f;
  }

} // end namespace GRINS

#endif // End if GRINS_HAVE_MASA
