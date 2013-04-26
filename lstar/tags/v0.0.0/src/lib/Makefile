################################################################
#
# $Author$
# $LastChangedDate$
# $Revision$
# $Id$
#
################################################################


################################################################
CC		=
CXX		=
LD		=$(CXX)
AR		=
CFLAGS	=-Wall -O3
CXXFLAGS=$(CFLAGS)

PROJECT	= lstar
PROJECT_ROOT_DIR=../..
################################################################

PROJECT_SRC_DIR=$(PROJECT_ROOT_DIR)/src
LIB_SRC_DIR	=$(PROJECT_SRC_DIR)/lib
LIB_TARGET=$(PROJECT)xx

HEADERS = AffineTransform.h Contour.h Coordinator.h Definitions.h \
		  FieldLine.h FieldLineTable.h FieldModel.h Geopack.h Key.h \
		  MagneticFlux.h ODESolver.h Particle.h ParticleTemplate.h \
		  Point.h T02.h T89.h T96.h TS05.h TSExternalField.h \
		  TSFieldComponent.h TSFieldModel.h Thread.h UBKLstarxx.h

SOURCES	= AffineTransform.cpp Contour.cpp Coordinator.cpp FieldLine.cpp \
		  FieldLineTable.cpp FieldModel.cpp Geopack.cpp Key.cpp \
		  MagneticFlux.cpp ODESolver.cpp Particle.cpp ParticleTemplate.cpp \
		  Point.cpp T02.cpp T89.cpp T96.cpp TS05.cpp TSExternalField.cpp \
		  TSFieldModel.cpp Thread.cpp

OBJECTS =$(SOURCES:.cpp=.o)

LIB_NAME=lib$(LIB_TARGET).a

all: $(LIB_TARGET)

$(LIB_TARGET): $(LIB_NAME)
	@date

$(LIB_NAME): $(OBJECTS)
	$(AR) rcs $@ $^

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $^

install: all
	@date
	mkdir -p $(INSTALL_DIR)/include
	cp $(HEADERS) $(INSTALL_DIR)/include/
	mkdir -p $(INSTALL_DIR)/lib
	cp $(LIB_NAME) $(INSTALL_DIR)/lib/

.PHONY: clean $(LIB_TARGET) install

clean:
	rm -rf $(OBJECTS) $(LIB_NAME)
