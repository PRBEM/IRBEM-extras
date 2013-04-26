################################################################
#
# $Author$
# $LastChangedDate$
# $Revision$
# $Id$
#
################################################################


################################################################
PROJECT	= lstar
PROJECT_ROOT_DIR=../..

PROJECT_SRC_DIR=$(PROJECT_ROOT_DIR)/src
LIB_SRC_DIR	=$(PROJECT_SRC_DIR)/lib
DLM_SRC_DIR	=$(PROJECT_SRC_DIR)/dlm
DLM_TARGET  =idl_$(PROJECT)
LIB_TARGET	=$(PROJECT)xx

LIB_PATH = $(LIB_SRC_DIR)/lib$(LIB_TARGET).a

IDL_DIR =
CC		=
CXX		=
DEBUG   =
C_FLAGS		= -O3 -I$(IDL_DIR)/external/include -c -Wall $(DEBUG) $(CFLAGS) -I$(LIB_SRC_DIR)
LD 		= $(CXX)
SHELL		= /bin/sh
X_CFLAGS	=
X_LD_FLAGS	=
X_LD_POST	= -lpthread -lm -O3 $(LIB_PATH)
SO_EXT		=so
################################################################

SOURCES=ubk_dlm_c.cpp

OBJECTS=$(SOURCES:.cpp=.o)

# The following is the default entry point. This section will determine 
# what system we are on and call make again with the correct flags set.

all : 
	@echo "OS type detected: "`uname`
	@case `uname` in \
	   "SunOS") if [ `/bin/uname -p` != i386 ]; then \
			echo "" ;\
			echo "We can build either a 32-bit or a 64-bit" ;\
			echo "version of this DLM. Please use the command" ;\
			echo "          make sun32" ;\
			echo "for a 32-bit version, or" ;\
			echo "          make sun64" ;\
			echo "for the 64-bit version. It is important that" ;\
			echo "you build the correct version for the IDL" ;\
			echo "binary you are using. Otherwise the DLM will" ;\
			echo "fail to load." ;\
		    else \
			make $(DLM_TARGET) \
			"X_CFLAGS=-K pic" \
			"X_LD_FLAGS=-B symbolic -z nodefs -G" ; \
		    fi \
		;; \
	   "AIX") make $(DLM_TARGET) \
			"LD=cc" \
			"SO_EXT=a" \
			"X_CFLAGS=-qro -qroconst" \
			"X_LD_FLAGS= -bM:SRE -bnoentry -btextro -bE:$(DLM_TARGET).export -bI:$(IDL_DIR)/external/idl.export" ;;\
	   "Darwin" ) make $(DLM_TARGET) \
			"X_CFLAGS= -fPIC" \
			"X_LD_FLAGS= -flat_namespace -undefined suppress -bundle" ;; \
	   "HP-UX") make $(DLM_TARGET)  \
			"SO_EXT= sl" \
			"X_CFLAGS=+z +e -Aa -D_HPUX_SOURCE" \
			"X_LD_FLAGS= -b -Brestricted +s" ;; \
	   "IRIX" )  make $(DLM_TARGET) \
			"X_CFLAGS=-n32 -mips3" \
		       "X_LD_FLAGS=-Bsymbolic -n32 -mips3 -shared -all -rdata_shared" ;;\
	   "IRIX64" )  make $(DLM_TARGET) \
		       "X_CFLAGS=-n32 -mips3" \
		       "X_LD_FLAGS=-Bsymbolic -n32 -mips3 -shared -all -rdata_shared" ;;\
	   "OSF1" )  make $(DLM_TARGET) \
			"X_CFLAGS=-float" \
			"X_LD_FLAGS=-expect_unresolved '*' -shared -all" ;;\
	   "Linux" ) make $(DLM_TARGET) \
			"X_CFLAGS= -fPIC" \
			"X_LD_FLAGS= -shared -Bsymbolic --warn-once -noinhibit-exec " ; \
		;; \
	   *) echo "This system is not supported" ;; \
       esac



# We can't tell if you want a 32-bit or a 64-bit version of the DLM
# on a Sparc/Solaris system because we could build either. Hence, we
# force you to specify the target on the make command line, leading you to
# here
sun32 :
	make $(DLM_TARGET) "X_CFLAGS=-K pic" "X_LD_FLAGS=-B symbolic -z nodefs -G"

sun64 :
	make $(DLM_TARGET) "X_CFLAGS=-K pic -xtarget=ultra -xarch=v9" \
		"X_LD_FLAGS=-B symbolic -z nodefs -G"



# The following is the entry point for the recursive build

$(DLM_TARGET) : $(DLM_TARGET).$(SO_EXT)
	@date

$(DLM_TARGET).$(SO_EXT) : $(OBJECTS) $(LIB_PATH)
	$(LD) $(X_LD_FLAGS) -o $(DLM_TARGET).$(SO_EXT) $(X_LD_POST) $^

$(LIB_PATH): $(LIB_TARGET)
	$(MAKE) PROJECT=$(PROJECT) -C $(LIB_SRC_DIR) $^

.cpp.o :
	$(CXX) $(C_FLAGS) $(X_CFLAGS) $^

install: all
	@date
	mkdir -p $(INSTALL_DIR)/dlm
	cp $(DLM_TARGET).$(SO_EXT) $(INSTALL_DIR)/dlm/
	cp $(DLM_TARGET).dlm $(INSTALL_DIR)/dlm/
	cp -r ubk/*.pro $(INSTALL_DIR)/dlm/

.PHONY: clean $(LIB_TARGET) install

clean : 
	rm -f $(OBJECTS) $(DLM_TARGET).so $(DLM_TARGET).sl $(DLM_TARGET).a \
		so_locations
