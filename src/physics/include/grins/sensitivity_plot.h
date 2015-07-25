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


#ifndef GRINS_SENSITIVITY_PLOT_H
#define GRINS_SENSITIVITY_PLOT_H

//GRINS
#include "grins/inc_navier_stokes_base.h"

namespace GRINS
{

  //! Physics class for Sensitivity Plot
  /*!
    This class simply solves the system int S phi_j = int partialR\partialp(u_h, phi_j;z)
    for S = sum S_i phi_i
  */

  class SensitivityPlot : public IncompressibleNavierStokesBase
  {
  public:

    SensitivityPlot(const std::string& physics_name, const GetPot& input);

    ~SensitivityPlot();

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    // Context initialization
    // We also need to init the contexts for the adjoint variables (see init_context
    // of adjoint example 5 in libMesh)
    virtual void init_context( AssemblyContext& context );

    // Time dependent part(s)
    virtual void element_time_derivative( bool compute_jacobian,
					  AssemblyContext& context,
					  CachedValues& cache );

  private:

    SensitivityPlot();

  };

} // End namespace block

#endif // GRINS_SENSITIVITY_PLOT_H
