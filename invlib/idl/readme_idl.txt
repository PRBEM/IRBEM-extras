;******************************************************************************
File readme_idl.txt

;******************************************************************************
A brief HOWTO on making an IDL - callable shareable object for Paul O'Briens
invlib C library.

Author: Reiner Friedel, LANL

Date: 12 February 2009  
;******************************************************************************

Summary:

As part of the invlib distribution you should see the following IDL file:

make_dll_invlib_obrien.pro

which contains two routines:

make_dll_invlib_obrien
test_invlib_obrien

make_dll_invlib_obrien uses native IDL routines to build the shareable object
file, while test_invlib_obrien is an IDL implementation of Paul's specinv_test
routine which uses the IDL call_external method to access ana_spec_inv.


;******************************************************************************

Notes:

In calling ana_spec_inv from IDL, the passing of paramters needs to be treated
as strict (IDL is notoriously lax in this regard). 

You need to specify the type and size of variables EXACTLY matching the
variables declared in C.

test_invlib_obrien was implemented on a 64bit linux PC. Please check that on
your architecture the varibale declarions in IDL still match the C conventions
default sizes for C ints, long ints etc can vary from OS architecturevto OS
architecture... 

See comments in make_dll_invlib_obrien.pro for details.

Special note for passing strings:

C treats srings a null terminated byte arrays. To pass strings, declare them
in the following way in IDL:

str_variable = [byte('puke.dat'), 0b]
 

;******************************************************************************

Sample run output:

The idl routines in make_dll_invlib_obrien.pro need to be run from the
directory of the invlib distribution  (or use keyword CODE_BASE to specify the
full path of the invlib location).

 

IDL> .r make_dll_invlib_obrien
% Compiled module: MAKE_DLL_INVLIB_OBRIEN.
% Compiled module: TEST_INVLIB_OBRIEN.
IDL> make_dll_invlib_obrien
% MAKE_DLL_INVLIB_OBRIEN: Shared library created:
/a/toaster-g3/vol/home/u/friedel/idl/call_external/C/invlib_obrien/invlib_idl.so
% MAKE_DLL_INVLIB_OBRIEN: See info in dir compile_temp for build details
% MAKE_DLL_INVLIB_OBRIEN: Use test_invlib_obrien for IDL version of specinv_test
IDL> test_invlib_obrien
specinf test: NC= 5, NE=100, NEout=50
Selecting Gaussian Relative Error penalty function for y[0]=2.09324e+06
Selecting Gaussian Relative Error penalty function for y[1]=463346
Selecting Gaussian Relative Error penalty function for y[2]=130927
Selecting Gaussian Relative Error penalty function for y[3]=24851.7
Selecting Gaussian Relative Error penalty function for y[4]=3733.96
Trying Power-Law
optimize invoked with minimizer vector_bfgs2, MaxIter=100
optimize: 1/100: 2.32315558193513e+03
optimize: completed after 2/100: iteration is not making progress towards solution
Fit results, q:
1.000000
3.000000

Trying Exponential
optimize invoked with minimizer vector_bfgs2, MaxIter=100
optimize: completed after 1/100: iteration is not making progress towards solution
Fit results, q:
1.000000
-1.000000

RESULT          LONG      =            1
% TEST_INVLIB_OBRIEN: Return Code:           1
IDL>