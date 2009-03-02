/* macros for ana_spec_inv */
/* function ana_spec_inv defined in invlib.h */

#ifndef SPECINV_H

/* macros for analytical spectral inversion, match these to instructions in header file */
/* macros for function bitmaps */
#define ASI_FXN_PL  (1)
#define ASI_FXN_EXP (2)
#define ASI_FXN_RM  (4)
#define ASI_FXN_PLE  (8)
#define ASI_FXN_RM2  (16)
#define ASI_MAX_POW2 (4)
#define ASI_FXN_ALL (ASI_FXN_PL | ASI_FXN_EXP | ASI_FXN_RM | ASI_FXN_PLE | ASI_FXN_RM2)

/* maximum number of free parameters for any fit */
#define ASI_MAX_NQ (10)

/* returns offset for first support_data for i'th fit */
#define ASI_SD_START(i) ((ASI_MAX_NQ+1)*i)

#define SPECINV_H 1
#endif
