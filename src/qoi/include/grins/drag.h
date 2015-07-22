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


#ifndef GRINS_DRAG_H
#define GRINS_DRAG_H

// GRINS
#include "grins/qoi_base.h"
#include "grins/variable_name_defaults.h"

namespace GRINS
{
  //! Drag QoI
  /*!
    This class implements a drag QoI that can be used to both compute
    QoI values and drive QoI-based adaptive refinement. Currently, this QoI
    is only implemented in 2D and will error if it detects a three-dimensional
    problem.
   */
  template<class Viscosity>
  class Drag : public QoIBase
  {
  public:

    //! Constructor
    /*! Constructor takes GetPot object to read any input options associated
        with this QoI */
    Drag( const std::string& qoi_name, const GetPot& input );

    virtual ~Drag();

    //! Required to provide clone (deep-copy) for adding QoI object to libMesh objects.
    virtual QoIBase* clone() const;

    virtual bool assemble_on_interior() const;

    virtual bool assemble_on_sides() const;

    //! Initialize local variables
    /*! Any local variables that need information from libMesh get initialized
        here. For example, variable indices. */
    virtual void init( const GetPot& input, MultiphysicsSystem& system );

    virtual void init_context( AssemblyContext& context );

    //! Compute the qoi value.
    /*! Currently, only implemented for 2D. Assumes that the vorticity will be
        computed over area of input subdomain id. Drag computed as an interior
	integral, see Bangerth pg */
    virtual void element_qoi( AssemblyContext& context,
                              const unsigned int qoi_index );

    //! Compute the qoi derivative with respect to the solution.
    /*! Drag is a boundary flux, its contribution to the adjoint RHS
      is via a Dirichlet boundary condition. */

    // Registers all parameters in this QoI
    /* virtual void register_parameter */
    /*   ( const std::string & param_name, */
    /*     libMesh::ParameterMultiPointer<libMesh::Number> & param_pointer ) */
    /* const; */

  protected:

    //! u-velocity component variable index
    VariableIndex _u_var;

    //! v-velocity component variable index
    VariableIndex _v_var;

    //! pressure variable index
    VariableIndex _p_var;

    //! Viscosity object
    //Viscosity _mu;

    //! List of sumdomain ids for which we want to compute this QoI
    //std::set<libMesh::subdomain_id_type> _subdomain_ids;

  private:
    //! User never call default constructor.
    Drag();

  };

  template<class Mu>
    inline
  bool Drag<Mu>::assemble_on_interior() const
  {
    return true;
  }

  template<class Mu>
    inline
  bool Drag<Mu>::assemble_on_sides() const
  {
    return false;
  }
}
#endif //GRINS_DRAG_H
