// This file is meant for defining shared macros

#ifndef BSACL_MACROS_H_
#define BSACL_MACROS_H_


#ifdef __BSACL_USE_CONST_EXT_POINTERS
#  define __EXT_PTR_CONST const
#else
#  define __EXT_PTR_CONST
#endif


#endif // BSACL_MACROS_H_