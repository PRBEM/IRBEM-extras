SHELL=/bin/sh
NULL=true

# Lines starting with #. will be printed as helps
#. make inversion library (invlib.so or invlib.dll)
#. (depends on libnnlib.a in ../nnlib)
#. make help - display this help
#. make invlib.so - make unix/cygwin shared object
#. make invlib.dll - make windows shared object (from MSYS/MINGW only)
#. make libinvlib.a - make precompiled library
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
IFLAGS=-I/usr/local/include -I./ -I../nnlib/
LFLAGS=-L/usr/local/lib -L./ -L../nnlib/
CFLAGS= -Wall $(IFLAGS) $(LFLAGS)
NNLIB=../nnlib/libnnlib.a
INVLIB=libinvlib.a

INVLIBFILES=invutil.o specinv.o omni2uni.o o2u_TEM1.o optim.o 
LIBS=-lgsl -lgslcblas -lnnlib

help:
	@grep '^.\.' Makefile

.c.o:
	$(CC) $(CFLAGS) -c $<

libinvlib.a: $(INVLIBFILES)
	ar rcs libinvlib.so $(INVLIBFILES)

specinv_test.exe:  specinv_test.o $(INVLIB)
	$(CC) $(CFLAGS) -o specinv_test.exe specinv_test.o $(LIBS) -linvlib -static

omni2uni_test.exe:  omni2uni_test.o $(INVLIB)
	$(CC) $(CFLAGS) -o omni2uni_test.exe omni2uni_test.o $(LIBS) -linvlib -static

invlib.dll: $(INVLIBFILES)
	$(CC) $(CFLAGS) -o invlib.dll $(INVLIBFILES) -Wl,--add-stdcall-alias -fno-second-underscore -w -shared $(LIBS) -static

invlib.so: $(INVLIBFILES)
	$(CC) $(CFLAGS) -o invlib.so $(INVLIBFILES) -Wl,--add-stdcall-alias -fno-second-underscore -shared $(LIBS) -static

clean:
	\rm -f *.o *.so *.dll *~ invlib.tar* *_test.exe
