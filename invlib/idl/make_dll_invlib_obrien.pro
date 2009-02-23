;******************************************************************************
;* PROCEDURE:   make_dll_invlib_obrien
;* 
;* DESCRIPTION:  
;*	This procedure is used to make the idl-callablel sharebable
;*      object for the Paul O'Biren's invlib c-library. 
;*      It uses the IDL procedure MAKE_DLL.
;*      run from location of invlib files. 
;*
;* INPUTS:
;*      none
;*
;* KEYWORDS:
;*      CODE_DIR -  set to location of invlib codes. 
;*                  Default current directory.
;*
;* MODIFICATION HISTORY:       
;*      written Feb 2009, Reiner Friedel
;******************************************************************************
PRO make_dll_invlib_obrien, CODE_DIR = code_dir, VERBOSE=verbose

IF keyword_set(CODE_DIR) THEN code_dir = CODE_DIR ELSE $
  cd, current = code_dir
code_dir = code_dir+'/'

OutputFile = 'invlib_idl'

InputFiles = file_basename(file_search(code_dir+'/*.c'), '.c')

compile_dir = code_dir+'compile_temp'
file_mkdir, compile_dir

export_rtns = ['ana_spec_inv']

MAKE_DLL, InputFiles,  OutputFile, export_rtns, $
          INPUT_DIR = code_dir, DLL_PATH = shlib, $
          COMPILE_DIR = compile_dir, OUTPUT_DIR = code_dir, $
          VERBOSE = verbose, SHOW_ALL_OUTPUT = verbose, /nocleanup, $
          EXTRA_LFLAGS='-lgsl -lgslcblas'

message, 'Shared library created:', /info
print, shlib

message, 'See info in dir compile_temp for build details', /info
message, 'Use test_invlib_obrien for IDL version of specinv_test', /info

END 


;******************************************************************************
;* PROCEDURE:   test_invlib_obrien
;* 
;* DESCRIPTION:  
;*	This procedure implements specinv_test.c in IDL
;*      using invlib_obrien.so with call_external in IDL.
;*
;* KEYWORDS:
;*      CODE_DIR -  set to location of invlib codes. 
;*                  Default current directory.
;* MODIFICATION HISTORY:       
;*      written Feb 2009, Reiner Friedel
;******************************************************************************
PRO test_invlib_obrien, CODE_DIR = code_dir


IF keyword_set(CODE_DIR) THEN code_dir = CODE_DIR ELSE $
  cd, current = code_dir
code_dir = code_dir+'/'

;read c, dc, Egrid, H, b from ascii file specinv_test.in1
fln = code_dir+'specinv_test.in1'

N_C = 0l & N_E = 0l & NE_out = 0l &  instr = '' & dmy = 0.0d

openr, u, fln, /get_lun

readf, u, N_C, N_E, NE_out

c = dblarr(N_C)
dc = dblarr(N_C)
Egrid = dblarr(N_E)
H = dblarr(N_C*N_E)
b = dblarr(N_C)

Eout = dblarr(NE_out)
flux = dblarr(NE_out)
dlogflux = dblarr(NE_out)

readf, u, instr
FOR i = 0, N_C-1 DO BEGIN & readf, u, dmy & c[i] = dmy & ENDFOR 
readf, u, instr
FOR i = 0, N_C-1 DO BEGIN & readf, u, dmy & dc[i] = dmy & ENDFOR 
readf, u, instr
FOR i = 0, N_E-1 DO BEGIN & readf, u, dmy & Egrid[i] = dmy & ENDFOR 
readf, u, instr
FOR i = 0, N_E*N_C-1 DO BEGIN & readf, u, dmy & H[i] = dmy & ENDFOR 
readf, u, instr
FOR i = 0, N_C-1 DO BEGIN & readf, u, dmy & b[i] = dmy & ENDFOR 
readf, u, instr
FOR i = 0, NE_out-1 DO BEGIN & readf, u, dmy & Eout[i] = dmy & ENDFOR 

close, u &  free_lun, u

print, n_c, n_e, NE_out, $
  format = "('specinf test: NC=',i2,', NE=',i3,', NEout=',i2)"

;setup parametrs for call to ana_spec_inv

int_prms = lon64arr(7)
real_prms = dblarr(3)

int_prms[0] = N_C
int_prms[1] = N_E
int_prms[2] = NE_out
int_prms[3] = 1+2                ;analtyical functions bitmap
int_prms[4] = 0                  ;minimizer, 0=BFGS, 3=NM
int_prms[5] = 100                ;maximumn # of iterations
int_prms[6] = 1                  ;1 = verbose to standard out

real_prms[0] = 0.511             ;electron rest energy, MeV
real_prms[1] = 100               ;E_break for proton PLE, MeV
real_prms[2] = 345               ;E0 for proton PLE, MeV

filename = 'puke.dat'
filename = [byte('puke.dat'), 0b]

lib_name = code_dir+'invlib_idl.so'

;NOTE!!!
;YOU HAVE to declare the calling paramters that get passed throguh
;call_extranlal exactly as what they are expected as in the routine
;called.
;IDL will segfault ungracefully IF there is a paramater mismatch.
;Here's the def of the parmaters on the C side:
; 
; int ana_spec_inv(const double *y, const double *dy, const double*Egrid, 
;                  const double *H, const double *b,
;                const long int *int_params, const double *real_params,
;                char *filename, double *Eout, double *flux, double *dlogflux);
;
; Also note that depending on your OS achitecture, the default size of
; the C variables may differ. E.G, under 64 bit linux a C long is 8
; bytes which is an IDL lon64, while on 32 bits it's 4 butes and this
; a IDL long....
; Also note that the ana_spec_inv call had to be altered to handle the
; "outFilePtr" argument. Originally this was passed in a a file handle
; pointer - for which IDL has no type equivalent that woukld work.
; This was changed to a filename, and the opening of the fiel is
; handled in ana_spec_inv (not gracefully, at each call the fiel gets
; opened and closed...). So hwen calling the rotuine from IDL, it
; would be good to avoud the Verbose = 3 setting for now. 

result = call_external(lib_name, 'ana_spec_inv', $
              c, dc, Egrid, H, b,$                 ;input prms
              int_prms, real_prms, filename,$      ;control prms
              Eout, flux, dlogflux, $              ;output prms
              /AUTO_GLUE, /UNLOAD, RETURN_TYPE = 3)

help, result
message, 'Return Code:'+string(result), /info

END 
