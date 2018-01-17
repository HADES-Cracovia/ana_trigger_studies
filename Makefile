APP_NAME      := trigger_studies
SOURCE_FILES  := analysis.cc core.cc 
#INSTALL_DIR   := ..

USES_RFIO     := no
USES_ORACLE   := no
USES_GFORTRAN := yes

include $(HADDIR)/hades.def.mk

LIB_DIRS += $(PLUTODIR) ${HOME}/usr/lib64
INC_DIRS += ${HOME}/usr/include

.PHONY:  default
#default: clean build install
default: build install

CPP_FLAGS += -std=c++11

include $(HADDIR)/hades.app.mk
