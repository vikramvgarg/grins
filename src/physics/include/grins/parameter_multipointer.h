//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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


// The ParameterMultiPointer Class offers the capability to organize and modify different 
// instances of a variable with the same name and meaning (for e.g. density 'rho') existing
// as member variables of different Physics ('rho' for example is declared in Navier-Stokes 
// as well as Heat Transfer)

// A ParameterHandler object will instruct each instantiated Physics to register its instance
// of a variable with the corresponding ParameterMultiPointer object. For e.g. a ParameterHandler 
// will tell the NavierStokes and HeatTransfer Physics to register their instance of 'rho'
// with a single RhoMultiPointer object. The ParameterHandler can then ask RhoMultiPointer 
// to modify all the different instances of 'rho' at once.

#ifndef GRINS_PARAMETER_MULTIPOINTER_H
#define GRINS_PARAMETER_MULTIPOINTER_H

// GRINS include files

// libMesh include files


namespace GRINS
{
 public:
  /**
   * Default constructor: Needs a ParameterVector to hold the various instances
   * of the parameter we are pointing to
   */
  ParameterMultiPointer (libMesh::ParameterVector)
