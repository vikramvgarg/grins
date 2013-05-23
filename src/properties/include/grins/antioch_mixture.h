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

#ifndef GRINS_ANTIOCH_MIXTURE_H
#define GRINS_ANTIOCH_MIXTURE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// libMesh
#include "libmesh/libmesh_common.h"

// Antioch
#include "antioch/vector_utils.h"
#include "antioch/chemical_mixture.h"

// Boost
#include <boost/scoped_ptr.hpp>

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  class AntiochMixture
  {
  public:

    AntiochMixture( const GetPot& input );
    ~AntiochMixture();

    libMesh::Real M( unsigned int species ) const;

    libMesh::Real M_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    libMesh::Real R( unsigned int species ) const;

    libMesh::Real R_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    libMesh::Real X( unsigned int species, libMesh::Real M, libMesh::Real mass_fraction ) const;

    void X( libMesh::Real M, const std::vector<libMesh::Real>& mass_fractions, 
	    std::vector<libMesh::Real>& mole_fractions ) const;

    unsigned int species_index( const std::string& species_name ) const;

    std::string species_name( unsigned int species_index ) const;

  protected:

    boost::scoped_ptr<Antioch::ChemicalMixture<double> > _antioch_gas;

  private:

    AntiochMixture();

  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  libMesh::Real AntiochMixture::M( unsigned int species ) const 
  {
    return _antioch_gas->M(species);
  }

  inline
  libMesh::Real AntiochMixture::M_mix( const std::vector<libMesh::Real>& mass_fractions ) const 
  {
    return _antioch_gas->M(mass_fractions);
  }

  inline
  libMesh::Real AntiochMixture::R( unsigned int species ) const 
  {
    return _antioch_gas->R(species);
  }

  inline
  libMesh::Real AntiochMixture::R_mix( const std::vector<libMesh::Real>& mass_fractions ) const 
  {
    return _antioch_gas->R(mass_fractions);
  }

  inline
  libMesh::Real AntiochMixture::X( unsigned int species, const libMesh::Real M,
                                   const libMesh::Real mass_fraction ) const
  {
    return _antioch_gas->X(species,M,mass_fraction);
  }

  inline
  void AntiochMixture::X( libMesh::Real M,
                          const std::vector<libMesh::Real>& mass_fractions, 
                          std::vector<libMesh::Real>& mole_fractions ) const
  {
    _antioch_gas->X(M,mass_fractions,mole_fractions);
    return;
  }

  inline
  unsigned int AntiochMixture::species_index( const std::string& species_name ) const
  {
    return _antioch_gas->active_species_name_map().find(species_name)->second;
  }

  inline
  std::string AntiochMixture::species_name( unsigned int /*species_index*/ ) const
  {
    libmesh_not_implemented();
    return "dummy";
  }
  
} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_MIXTURE_H
