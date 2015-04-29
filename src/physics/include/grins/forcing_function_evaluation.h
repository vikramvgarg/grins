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


#ifndef GRINS_FORCING_FUNCTION_EVALUATION_H
#define GRINS_FORCING_FUNCTION_EVALUATION_H

// GRINS
#include "grins/physics.h"

// MASA
#ifdef GRINS_HAVE_MASA
#include "masa.h"

namespace GRINS
{
  class ForcingFunctionEvaluation : public Physics
   {
    public:

    ForcingFunctionEvaluation( const std::string& physics_name, const GetPot& input );

    ~ForcingFunctionEvaluation();

    libMesh::Real compute_masa_forcing( libMesh::Real x, libMesh::Real y, std::string& unknown_name);

    virtual void init_variables( libMesh::FEMSystem* system );

    };

} // end namespace GRINS

#endif // End if GRINS_HAVE_MASA

#endif // GRINS_FORCING_FUNCTION_EVALUATION_H
