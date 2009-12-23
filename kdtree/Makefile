SHELL=/bin/sh
NULL=true

# Lines starting with #. will be printed as helps
#. make kdtree library (kdtree.so or kdtree.dll)
#. make help - display this help
#. make kdtree.dll - make windows shared object (from MSYS/MINGW only)
#. make kdtree.so - make unix shared object
#. make kdtree-omp - make unix shared object w/ openmp support (replaced kdtree.so)
#. make kdtree_oct.oct - make oct file
#. make dist - make a tarball for distribution
#. make clean - remove compiled files, .tar files, and emacs ~ files

#------------------------------------------------------------------------------
# macros.
#------------------------------------------------------------------------------

# note: command line macros override what's defined here
MAKE=make
SO=so
CC=gcc
IFLAGS=
LFLAGS=
CFLAGS= -Wall $(IFLAGS) $(LFLAGS) -fPIC

LIBFILES=kdtree.o
LIBS=

CFILES=kdtree.c
HFILES=kdtree.h
MFILES=
OCTFILES=kdtree_oct.cc
TARMISC=Makefile
TARFILES=$(CFILES) $(HFILES) $(MFILES) $(TARMISC)

ifneq (,$(findstring 64,$(MACHTYPE)))
# 64-bit case
  OCTDEFS=-DLINUX64
  # enable use of OpenMP in octfile
  export XTRA_CFLAGS=-fopenmp -DUSEOMP
  export XTRA_CXXFLAGS=-fopenmp -DUSEOMP
  export LDFLAGS=-fopenmp -DUSEOMP
else
# 32-bit case
  OCTDEFS=
endif

help:
	@grep '^.\.' Makefile

.c.o:
	$(CC) $(CFLAGS) -c $<

kdtree.dll: $(LIBFILES)
	$(CC) $(CFLAGS) -o kdtree.dll $(LIBFILES) -Wl,--add-stdcall-alias -w -shared $(LIBS) -static

kdtree.so: $(LIBFILES) $(HFILES)
#	$(CC) $(CFLAGS) -o kdtree.so $(LIBFILES) -Wl,--add-stdcall-alias -shared $(LIBS) -static
	$(CC) $(CFLAGS) -shared -o kdtree.so $(LIBFILES) $(LIBS)

kdtree-omp: $(CFILES) $(HFILES)
	$(CC) $(CFLAGS) -fopenmp -DUSEOMP -shared -o kdtree.so $(CFILES) $(LIBS)

kdtree_oct.oct: $(CFILES) $(HFILES) $(OCTFILES)
	\rm -f *.o ; \
	mkoctfile  $(OCTDEFS) -o kdtree_oct.oct $(OCTFILES) $(CFILES) ; \
	\rm -f *.o

dist: 
	\rm -f kdtree.tar.gz kdtree.tar ; \
	\tar cvf kdtree.tar $(TARFILES)   ; \
	\gzip -f kdtree.tar 

clean:
	\rm -f *.o *.so *.dll *~ kdtree.tar*
