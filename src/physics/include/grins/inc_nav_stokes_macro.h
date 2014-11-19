#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/constant_conductivity.h"
#include "grins/parsed_conductivity.h"

#ifndef GRINS_INC_NAV_STOKES_MACRO_H
#define GRINS_INC_NAV_STOKES_MACRO_H

#define INSTANTIATE_INC_NS_SUBCLASS(class_name) \
  template class GRINS::class_name<GRINS::ConstantViscosity, GRINS::ConstantConductivity>; \
  template class GRINS::class_name<GRINS::ConstantViscosity, GRINS::ParsedConductivity>; \
  template class GRINS::class_name<GRINS::ParsedViscosity, GRINS::ConstantConductivity>; \
  template class GRINS::class_name<GRINS::ParsedViscosity, GRINS::ParsedConductivity>

#endif // GRINS_INC_NAV_STOKES_MACRO_H
