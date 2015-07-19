#ifndef GRINS_DRAG_MACRO_H
#define GRINS_DRAG_MACRO_H

#define INSTANTIATE_DRAG_SUBCLASS(class_name) \
template class GRINS::class_name<GRINS::ConstantViscosity>; \
 template class GRINS::class_name<GRINS::ParsedViscosity>; \
template class GRINS::class_name<GRINS::SpalartAllmarasViscosity<GRINS::ConstantViscosity> >

#endif // GRINS_DRAG_MACRO_H
