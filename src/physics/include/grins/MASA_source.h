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

// Class definition to support manufactured solution verification via
// the MASA library

#ifndef __masa_source_h__
#define __masa_source_h__

#ifdef GRINS_HAVE_MASA
#  include "libmesh/point.h"
#  include "grins/physics.h"
#  include "masa.h"
#endif

namespace GRINS {

  //! Interface for MASA manufactured solutions
  /*! 
    This physics provides access to MASA's library of manufactured solutions. To use this
    capability, the user needs to include "MASA source" to the list of instantiated physics
    ,specify the source_name string to the particular solution they want to use and specify
    the relevant parameters
  */

  
