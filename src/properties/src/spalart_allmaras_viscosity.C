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


// This class
#include "grins/spalart_allmaras_viscosity.h"
#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"

//GRINS
#include "grins/grins_physics_names.h"
#include "grins/turbulent_viscosity_macro.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/parsed_function.h"

namespace GRINS
{
  template<class Mu>
  SpalartAllmarasViscosity<Mu>::SpalartAllmarasViscosity( const GetPot& input ):
    ParameterUser("SpalartAllmarasViscosity"),
    _mu(input),
    _turbulence_vars(input, spalart_allmaras),
    _sa_params(input)
  {
    if( !input.have_variable("Materials/Viscosity/mu") )
      {
        std::cerr<<"No viscosity has been specified."<<std::endl;

        libmesh_error();
      }
    return;
  }

  template<class Mu>
  void SpalartAllmarasViscosity<Mu>::register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
    const
  {
    ParameterUser::register_parameter(param_name, param_pointer);
    this->_mu.register_parameter(param_name, param_pointer);
    this->_sa_params.register_parameter(param_name, param_pointer);
  }

  template<class Mu>
  void SpalartAllmarasViscosity<Mu>::init( libMesh::FEMSystem* system )
  {
    this->_turbulence_vars.init(system);
  }

  template<class Mu>
  libMesh::Real SpalartAllmarasViscosity<Mu>::operator()(AssemblyContext& context, unsigned int qp) const
  {
    // The physical viscosity
    libMesh::Real mu_physical = this->_mu(context, qp);

    // The unscaled turbulent viscosity (the nu the SA physics solves for)
    libMesh::Real nu = context.interior_value(this->_turbulence_vars.nu_var(),qp);

    // Assert that mu_value is greater than 0
    if(nu < 0.0)
      {
        libmesh_warning("Negative turbulent viscosity encountered !");

        // We are using a negative S-A model, so will set eddy viscosity to zero
        // if the turbulent viscosity nu < 0.0
        nu = 0.0;
      }

    // Step 1
    libMesh::Real chi = nu/mu_physical;

    // Step 2
    libMesh::Real fv1 = _sa_params.fv1(chi);

    // Step 3
    libMesh::Real mu_turbulent = nu*fv1;

    // Compute the value of the total viscosity and return it
    libMesh::Number mu_value = mu_turbulent + mu_physical; // Turbulent viscosity + physical viscosity

    return mu_value;
  }

} // namespace GRINS

INSTANTIATE_TURBULENT_VISCOSITY_SUBCLASS(SpalartAllmarasViscosity);
