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


CFLAGS=-g -DDEBUG
LDFLAGS=-lm -pthread
# #######################################################

PROJECT_SRC_DIR=$(PROJECT_ROOT_DIR)/src
LIB_SRC_DIR	=$(PROJECT_SRC_DIR)/lib
EXAMPLE_SRC_DIR	=$(PROJECT_SRC_DIR)/examples
EXAMPLE_TARGET=example_$(PROJECT)
LIB_TARGET	=$(PROJECT)xx

LIB_NAME=lib$(LIB_TARGET).a
LIB_PATH=$(LIB_SRC_DIR)/$(LIB_NAME)


SOURCES=main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
TARGET=$(EXAMPLE_TARGET)


all: $(TARGET)

$(TARGET): $(OBJECTS) $(LIB_PATH)
	$(CXX) $(LDFLAGS) -o $@ $^

.cpp.o:
	$(CXX) -I$(LIB_SRC_DIR) $(CFLAGS) -c $<

$(LIB_PATH): $(LIB_TARGET)
	$(MAKE) PROJECT=$(PROJECT) -C $(LIB_SRC_DIR) $^

install: all
	@date
	mkdir -p $(INSTALL_DIR)/examples
	cp $(EXAMPLE_TARGET) $(INSTALL_DIR)/examples/

.PHONY: $(LIB_TARGET)

clean:
	rm -rf $(OBJECTS) $(EXAMPLE_TARGET)
