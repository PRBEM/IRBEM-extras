SHELL=/bin/sh
NULL=true

# Lines starting with #. will be printed as helps
#. make neural network library (nnlib.so or nnlib.dll)
#. make help - display this help
#. make nnlib.so - make unix/cygwin shared object
#. make nnlib.dll - make windows shared object (from MSYS/MINGW only)
#. make nntrain.exe - make executable nntrain.exe (from MSYS/MINGW only)
#. make nntrain.x - make executable nntrain.x (unix/cygwin)
#. make nntrain.mpix - make MPI executable nntrain.mpix (unix)
#. make nntrain.ompx - make executable with OpenMP enabled (unix)
#. make dist - make a tarball for distribution
#. make clean - remove compiled files, nnlib.tar files, and emacs ~ files

#------------------------------------------------------------------------------
# macros.
#------------------------------------------------------------------------------

# note: command line macros override what's defined here
MAKE=make
SO=so
CC=gcc
MPICC=mpicc
ifeq ($(HOST),coto-ago1)
  IFLAGS= -I./
  LFLAGS= -L./
  CFLAGS= $(IFLAGS) $(LFLAGS) -DOLDGSL
  OMPCC = gcc
else
  IFLAGS=-I/usr/local/include -I./
  LFLAGS=-L/usr/local/lib -L./
  CFLAGS= -Wall $(IFLAGS) $(LFLAGS)
  ifeq ($(HOST),fellowship.aero.org)
	OMPCC= gcc42
#  else
#	# cygwin
#	OMPCC= gcc-4
  endif
endif



NNLIBCFILES=nn.c nnio.c matrix_inv.c
LIBS=-lgsl -lgslcblas -lm

TARCFILES=nn.c nnio.c matrix_inv.c nntrain.c
TARHFILES=nn.h nnlib.h subscripts.h matrix_inv.h
TARMISC=Makefile nnlib.dll nntrain.exe nnlib.m doc/nnlib.pdf
TARFILES=$(TARCFILES) $(TARHFILES) $(TARMISC)

help:
	@grep '^.\.' Makefile

.c.o:
	$(CC) $(CFLAGS) -c $<

nntrain.exe:  nntrain.c $(NNLIBCFILES) 
	$(CC) -s -O3 $(CFLAGS) -o nntrain.exe nntrain.c $(NNLIBCFILES) -static $(LIBS)

# nntrain.x on unix
nntrain.x: nntrain.o $(NNLIBCFILES) 
	$(CC) -s -O3 $(CFLAGS) -o nntrain.x nntrain.c $(NNLIBCFILES) $(LIBS)

# nntrain.mpix -- enable mpi parallel features
nntrain.mpix: nntrain.c $(NNLIBCFILES)
	$(MPICC) -s -O3 $(CFLAGS) -DUSEMPI -o nntrain.mpix nntrain.c $(NNLIBCFILES) $(LIBS)

# nntrain.ompx -- enable OpenMP parallel features
nntrain.ompx: nntrain.c $(NNLIBCFILES)
	$(OMPCC) -s -O3 $(CFLAGS) -fopenmp -DUSEOMP -o nntrain.ompx nntrain.c $(NNLIBCFILES) $(LIBS)

nnlib.dll: $(NNLIBCFILES)
	$(CC) $(CFLAGS) -shared -o nnlib.dll $(NNLIBCFILES) -Wl,--add-stdcall-alias -fno-second-underscore -w -static $(LIBS)

nnlib.so: $(NNLIBCFILES)
	$(CC) $(CFLAGS) -fPIC -shared -o nnlib.so $(NNLIBCFILES) $(LIBS)

dist:
	\rm -f nnlib.tar.gz nnlib.tar ; \
	\tar cvf nnlib.tar $(TARFILES)   ; \
	\gzip -f nnlib.tar 

clean:
	\rm -f *.o *.so *.dll *~ nnlib.tar* nntrain.exe nntrain.x nntrain.mpix nntrain.ompx *.a
