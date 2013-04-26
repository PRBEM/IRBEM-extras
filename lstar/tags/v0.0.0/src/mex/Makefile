################################################################
#
# $Author$
# $LastChangedDate$
# $Revision$
# $Id$
#
################################################################


# ###################USER CHANGE#########################
PROJECT	= lstar
PROJECT_ROOT_DIR=../..


MEX=
MEXOPTS=
INSTALLDIR=.
DEBUG=
# #######################################################

PROJECT_SRC_DIR=$(PROJECT_ROOT_DIR)/src
LIB_SRC_DIR	=$(PROJECT_SRC_DIR)/lib
MEX_SRC_DIR	=$(PROJECT_SRC_DIR)/mex
MEX_TARGET=mex_$(PROJECT)
LIB_TARGET	=$(PROJECT)xx

LIB_NAME=lib$(LIB_TARGET).a


MEX_TSFIELD_SOURCES=UBKTSFieldcxx.cpp
MEX_TSFIELD_OBJECTS=$(MEX_TSFIELD_SOURCES:.cpp=.o)
MEX_TSFIELD_TARGETS=UBKTSFieldcxx

MEX_COTRANS_SOURCES=UBKGeopackCotranscxx.cpp
MEX_COTRANS_OBJECTS=$(MEX_COTRANS_SOURCES:.cpp=.o)
MEX_COTRANS_TARGETS=UBKGeopackCotranscxx

MEX_FIELDLINE_SOURCES=UBKFieldLinecxx.cpp
MEX_FIELDLINE_OBJECTS=$(MEX_FIELDLINE_SOURCES:.cpp=.o)
MEX_FIELDLINE_TARGETS=UBKFieldLinecxx

MEX_LSTAR_SOURCES=UBKLstarcxx.cpp
MEX_LSTAR_OBJECTS=$(MEX_LSTAR_SOURCES:.cpp=.o)
MEX_LSTAR_TARGETS=UBKLstarcxx

TARGETS=$(MEX_FIELDLINE_TARGETS) $(MEX_TSFIELD_TARGETS) \
		$(MEX_LSTAR_TARGETS) $(MEX_COTRANS_TARGETS)

all: $(LIB_TARGET) $(TARGETS)

$(MEX_TSFIELD_TARGETS): $(MEX_TSFIELD_OBJECTS)
	$(MEX) -f $(MEXOPTS) -cxx -outdir $(INSTALLDIR) $^ $(LIB_SRC_DIR)/$(LIB_NAME)

$(MEX_COTRANS_TARGETS): $(MEX_COTRANS_OBJECTS)
	$(MEX) -f $(MEXOPTS) -cxx -outdir $(INSTALLDIR) $^ $(LIB_SRC_DIR)/$(LIB_NAME)

$(MEX_FIELDLINE_TARGETS): $(MEX_FIELDLINE_OBJECTS)
	$(MEX) -f $(MEXOPTS) -cxx -outdir $(INSTALLDIR) $^ $(LIB_SRC_DIR)/$(LIB_NAME)

$(MEX_LSTAR_TARGETS): $(MEX_LSTAR_OBJECTS)
	$(MEX) -f $(MEXOPTS) -cxx -outdir $(INSTALLDIR) $^ $(LIB_SRC_DIR)/$(LIB_NAME)

.cpp.o:
	$(MEX) $(DEBUG) -I$(LIB_SRC_DIR) -f $(MEXOPTS) -cxx -c $<

$(LIB_TARGET):
	$(MAKE) PROJECT=$(PROJECT) -C $(LIB_SRC_DIR) $@

install: all
	@date
	mkdir -p $(INSTALL_DIR)/mex
	cp *.mex* $(INSTALL_DIR)/mex/
	cp -r +ubk $(INSTALL_DIR)/mex/

clean:
	rm -rf *.mexa64 *.mexmaci64 *.o
