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
DOXYGEN =doxygen

PROJECT	=lstar
PROJECT_ROOT_DIR=.
INSTALL_DIR=`pwd`/local
################################################################

PROJECT_SRC_DIR=$(PROJECT_ROOT_DIR)/src

LIB_TARGET	=$(PROJECT)xx
DLM_TARGET	=idl_$(PROJECT)
MEX_TARGET	=mex_$(PROJECT)
EXAMPLE_TARGET=example_$(PROJECT)

LIB_SRC_DIR	=$(PROJECT_SRC_DIR)/lib
DLM_SRC_DIR	=$(PROJECT_SRC_DIR)/dlm
MEX_SRC_DIR	=$(PROJECT_SRC_DIR)/mex
EXAMPLE_SRC_DIR	=$(PROJECT_SRC_DIR)/examples

LIB_INSTALL_DIR	=$(PROJECT_INSTALL_DIR)/lib
DLM_INSTALL_DIR	=$(PROJECT_INSTALL_DIR)/dlm
MEX_INSTALL_DIR	=$(PROJECT_INSTALL_DIR)/mex
EXAMPLE_INSTALL_DIR	=$(PROJECT_INSTALL_DIR)/examples


all: $(PROJECT)

$(PROJECT):
	@date
	$(MAKE) PROJECT=$(PROJECT) -C $(PROJECT_SRC_DIR) $@

install: all
	@date
	$(MAKE) PROJECT=$(PROJECT) INSTALL_DIR=$(INSTALL_DIR) -C $(PROJECT_SRC_DIR) install

doc: lstar.dox.conf
	doxygen $^

.PHONY: clean $(PROJECT) install lstar.dox.conf

clean:
	@date
	$(MAKE) PROJECT=$(PROJECT) -C $(PROJECT_SRC_DIR) clean
