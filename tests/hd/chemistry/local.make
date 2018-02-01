.PHONY: clean_krome

<<<<<<< HEAD
KROME_DIR := $(HOME)/codes/krome
KROME_OPTIONS := $(HOME)/codes/krome_models/OptionsUCL_CIE_Ripamonti.txt
PROJECT_NAME := test
KROME_BUILD_DIR := $(KROME_DIR)/build_$(PROJECT_NAME)
KROME_LIBS := lapack

INC_DIRS += $(KROME_BUILD_DIR)
LIBS += $(KROME_LIBS)

KOBJS := krome_all.o krome_user_commons.o opkda1.o opkda2.o opkdmain.o
KROME_OBJS := $(addprefix $(KROME_BUILD_DIR)/, $(KOBJS))
F90FLAGS += $(KROME_OBJS)

$(KROME_OBJS):
	cd $(KROME_DIR) && ./krome -options=$(KROME_OPTIONS) -project=$(PROJECT_NAME) -compact
=======
KROME_DIR := $(HOME)/git/krome
PROJECT_NAME := hdchem
KROME_BUILD_DIR := $(KROME_DIR)/build_$(PROJECT_NAME)

INC_DIRS += $(KROME_BUILD_DIR)

KOBJS := krome_all.o krome_user_commons.o opkda1.o opkda2.o opkdmain.o
KROME_OBJS := $(addprefix $(KROME_BUILD_DIR)/, $(KOBJS))

$(KROME_OBJS):
	cd $(KROME_DIR) && ./krome -test=hello -project=$(PROJECT_NAME) -compact
>>>>>>> 5dfed322b1e6e5f233c6fbf78d40c402ae617142
	$(MAKE) -C $(KROME_BUILD_DIR) gfortran

clean_krome:
	$(MAKE) -C $(KROME_BUILD_DIR) clean

# Only build KROME_OBJS when they do not exist
amrvac.o: | $(KROME_OBJS)
<<<<<<< HEAD

## Use the version below if you already ahve all the Krome files
## in KROME_BUILD_DIR and only have to make it.
# .PHONY: clean_krome
#
# KROME_BUILD_DIR := build_krome_test
# KROME_LIBS := lapack
#
# INC_DIRS += $(KROME_BUILD_DIR)
# LIBS += $(KROME_LIBS)
#
# KOBJS := krome_all.o krome_user_commons.o opkda1.o opkda2.o opkdmain.o
# KROME_OBJS := $(addprefix $(KROME_BUILD_DIR)/, $(KOBJS))
# F90FLAGS += $(KROME_OBJS)
#
# $(KROME_OBJS):
# 	$(MAKE) -C $(KROME_BUILD_DIR) gfortran
#
# clean_krome:
# 	$(MAKE) -C $(KROME_BUILD_DIR) clean
#
# # Only build KROME_OBJS when they do not exist
# amrvac.o: | $(KROME_OBJS)
=======
>>>>>>> 5dfed322b1e6e5f233c6fbf78d40c402ae617142
