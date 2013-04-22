//
//  ubk_dlm_c.h
//  UBJDevelopment
//
//  Created by Kyungguk Min on 4/11/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#ifndef UBJDevelopment_ubk_dlm_c_h
#define UBJDevelopment_ubk_dlm_c_h

#include "UBKLstarxx.h"
#include "idl_export.h"

//
// IDL DLM Loader
//
UBK_C_EXTERN int IDL_Load( void );

//
//
//
UBK_C_EXTERN void ubk_cotrans(int argc, IDL_VPTR *argv,char *argk);

//
//
//
UBK_C_EXTERN void ubk_ts_field(int argc, IDL_VPTR *argv,char *argk);

//
//
//
UBK_C_EXTERN void ubk_field_line(int argc, IDL_VPTR *argv,char *argk);

//
//
//
UBK_C_EXTERN void ubk_lstar(int argc, IDL_VPTR *argv,char *argk);

#endif
