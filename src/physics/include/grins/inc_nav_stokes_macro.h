#ifndef GRINS_INC_NAV_STOKES_MACRO_H
#define GRINS_INC_NAV_STOKES_MACRO_H

#define INSTANTIATE_INC_NS_SUBCLASS(class_name) \
template class GRINS::class_name<GRINS::ConstantViscosity>; \
template class GRINS::class_name<GRINS::ConstantConductivity>

template class GRINS::class_name<GRINS::ConstantViscosity>; \
template class GRINS::class_name<GRINS::ParsedConductivity>

template class GRINS::class_name<GRINS::ParsedViscosity>; \
template class GRINS::class_name<GRINS::ConstantConductivity>

template class GRINS::class_name<GRINS::ParsedViscosity>; \
template class GRINS::class_name<GRINS::ParsedConductivity>

#endif // GRINS_INC_NAV_STOKES_MACRO_H
