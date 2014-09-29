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

PROJECT	=lstar
PROJECT_ROOT_DIR=..
################################################################

LIB_TARGET	=$(PROJECT)xx
DLM_TARGET	=idl_$(PROJECT)
MEX_TARGET	=mex_$(PROJECT)
EXAMPLE_TARGET=example_$(PROJECT)

PROJECT_SRC_DIR=$(PROJECT_ROOT_DIR)/src
LIB_SRC_DIR	=$(PROJECT_SRC_DIR)/lib
DLM_SRC_DIR	=$(PROJECT_SRC_DIR)/dlm
MEX_SRC_DIR	=$(PROJECT_SRC_DIR)/mex
EXAMPLE_SRC_DIR	=$(PROJECT_SRC_DIR)/examples


all: $(PROJECT)

$(PROJECT): $(LIB_TARGET) $(DLM_TARGET) $(MEX_TARGET) $(EXAMPLE_TARGET)
	@date

$(LIB_TARGET):
	$(MAKE) PROJECT=$(PROJECT) -C $(LIB_SRC_DIR) $@

$(DLM_TARGET):
	$(MAKE) PROJECT=$(PROJECT) -C $(DLM_SRC_DIR) all

$(MEX_TARGET):
	$(MAKE) PROJECT=$(PROJECT) -C $(MEX_SRC_DIR) all

$(EXAMPLE_TARGET):
	$(MAKE) PROJECT=$(PROJECT) -C $(EXAMPLE_SRC_DIR) all

install: all
	@date
	$(MAKE) PROJECT=$(PROJECT) INSTALL_DIR=$(INSTALL_DIR) -C $(LIB_SRC_DIR) install
	$(MAKE) PROJECT=$(PROJECT) INSTALL_DIR=$(INSTALL_DIR) -C $(DLM_SRC_DIR) install
	$(MAKE) PROJECT=$(PROJECT) INSTALL_DIR=$(INSTALL_DIR) -C $(MEX_SRC_DIR) install
	$(MAKE) PROJECT=$(PROJECT) INSTALL_DIR=$(INSTALL_DIR) -C $(EXAMPLE_SRC_DIR) install

.PHONY: clean $(PROJECT) $(LIB_TARGET) $(DLM_TARGET) $(MEX_TARGETS) install $(EXAMPLE_TARGET)

clean:
	$(MAKE) PROJECT=$(PROJECT) -C $(DLM_SRC_DIR) clean
	$(MAKE) PROJECT=$(PROJECT) -C $(LIB_SRC_DIR) clean
	$(MAKE) PROJECT=$(PROJECT) -C $(MEX_SRC_DIR) clean
	$(MAKE) PROJECT=$(PROJECT) -C $(EXAMPLE_SRC_DIR) clean
