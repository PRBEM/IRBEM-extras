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

#define ASI_DE_INCLUDED (0)
#define ASI_DE_TRAPZ (1)
#define ASI_DE_PLATEAU (2)

/* maximum number of free parameters for any fit */
#define ASI_MAX_NQ (10)

/* returns offset for first support_data for i'th fit, i is zero-based */
#define ASI_SD_START(i,NY) ((2+(ASI_MAX_NQ)+(NY))*(i))

#define SPECINV_H 1
#endif
