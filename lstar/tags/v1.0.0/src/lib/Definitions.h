//
//  Definitions.h
//  UBJDevelopment
//
//  Created by Kyungguk Min on 3/30/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#ifndef UBJDevelopment_Definitions_h
#define UBJDevelopment_Definitions_h

//
// Convenient macros
//

//
// Extern macros
//
#if !defined(UBK_C_EXTERN)
#if defined(__cplusplus)
#define UBK_C_EXTERN extern "C"
#else
#define UBK_C_EXTERN extern
#endif
#endif

//
// Inline macros
//
#if !defined(UBK_INLINE)
#if defined(__GNUC__)
#define UBK_INLINE static __inline__ __attribute__((always_inline))
#elif defined(__MWERKS__) || defined(__cplusplus)
#define UBK_INLINE static inline
#elif defined(_MSC_VER)
#define UBK_INLINE static __inline
#elif TARGET_OS_WIN32
#define UBK_INLINE static __inline__
#endif
#endif

//
// Logical type. Here BOOL which is char type replaces the role of c++'s bool type.
//
typedef char BOOL;
#ifndef YES
#define YES 1
#endif

#ifndef NO
#define NO 0
#endif

#endif
