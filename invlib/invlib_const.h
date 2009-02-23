/*
  Spectral inversion library - macros/constants.

 */

#ifndef INVLIB_CONST_H

/* macros for return codes for all invlib functions */
#define INVLIB_SUCCESS         (1)
#define INVLIB_ERROR            (0)
#define INVLIB_ERR_NULL         (-101)
#define INVLIB_ERR_DATAEMPTY    (-102)
#define INVLIB_ERR_DATANAN      (-103)
#define INVLIB_ERR_NOCOUNTS     (-104)
#define INVLIB_ERR_NOFXN        (-201)
#define INVLIB_ERR_INVALIDFXN   (-202)
#define INVLIB_ERR_RME0         (-301)
#define INVLIB_ERR_PLE          (-302)
#define INVLIB_ERR_INVALIDMIN   (-401)
#define INVLIB_ERR_INVALIDITER  (-402)
#define INVLIB_ERR_VERBOSE      (-501)
#define INVLIB_ERR_OUTFILE      (-502)
#define INVLIB_ERR_IALPHA0      (-601)
#define INVLIB_ERR_INVALIDMETH  (-602)

#define INVLIB_CONST_H 1
#endif

