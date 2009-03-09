SHELL=/bin/sh
NULL=true

# Lines starting with #. will be printed as helps
#. make inversion library (invlib.so or invlib.dll)
#. (depends on libnnlib.a in ../nnlib)
#. make help - display this help
#. make invlib.so - make unix/cygwin shared object
#. make invlib.dll - make windows shared object (from MSYS/MINGW only)
#. make specinv_test.exe - make unix/cygwin executable specinv_test.exe
#. make omni2uni_test.exe - make unix/cygwin executable omni2uni_test.exe
#. make dist - make a tarball for distribution
#. make clean - remove compiled files, invlib.tar files, and emacs ~ files

#------------------------------------------------------------------------------
# macros.
#------------------------------------------------------------------------------

# note: command line macros override what's defined here
MAKE=make
SO=so
CC=gcc

ifeq ($(HOST),coto-ago1)
  IFLAGS= -I./ -I../nnlib/
  LFLAGS= -L./
  CFLAGS= -Wall $(IFLAGS) $(LFLAGS) -DOLDGSL
  OMPCC = gcc
else
  IFLAGS=-I/usr/local/include -I./ -I../nnlib/
  LFLAGS=-L/usr/local/lib -L./
  CFLAGS= -Wall $(IFLAGS) $(LFLAGS)
  ifeq ($(HOST),fellowship)
	OMPCC= gcc42
  else
	# cygwin
	OMPCC= gcc-4
  endif
endif

ifeq ($(OSTYPE),msys)
  STATIC= -static
else
 STATIC=
endif


NNLIBCFILES=../nnlib/matrix_inv.c ../nnlib/nn.c ../nnlib/nnio.c
INVLIBCFILES=invutil.c specinv.c omni2uni.c o2u_TEM1.c optim.c  $(NNLIBCFILES)

help:
	@grep '^.\.' Makefile

.c.o:
	$(CC) $(CFLAGS) -c $<

specinv_test.exe:  specinv_test.c  $(INVLIBCFILES)
	$(CC) $(CFLAGS) -o specinv_test.exe specinv_test.c $(INVLIBCFILES) $(STATIC) -lgsl -lgslcblas -lm

omni2uni_test.exe:  omni2uni_test.c  $(INVLIBCFILES)
	$(CC) $(CFLAGS) -o omni2uni_test.exe omni2uni_test.c $(INVLIBCFILES) $(STATIC) -lgsl -lgslcblas -lm

invlib.dll: $(INVLIBCFILES)
	$(CC) $(CFLAGS) -shared -o invlib.dll $(INVLIBCFILES) -Wl,--add-stdcall-alias -fno-second-underscore -w  $(STATIC) -lgsl -lgslcblas -lm 

invlib.so: $(INVLIBCFILES)
	$(CC) $(CFLAGS) -shared -o invlib.so $(INVLIBCFILES) $(STATIC)  -lgsl -lgslcblas -lm

clean:
	\rm -f *.o *.so *.dll *~ invlib.tar* *_test.exe *.a
