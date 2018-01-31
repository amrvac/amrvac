.PHONY: clean_krome

KROME_DIR := $(HOME)/git/krome
PROJECT_NAME := hdchem
KROME_BUILD_DIR := $(KROME_DIR)/build_$(PROJECT_NAME)

INC_DIRS += $(KROME_BUILD_DIR)

KOBJS := krome_all.o krome_user_commons.o opkda1.o opkda2.o opkdmain.o
KROME_OBJS := $(addprefix $(KROME_BUILD_DIR)/, $(KOBJS))

$(KROME_OBJS):
	cd $(KROME_DIR) && ./krome -test=hello -project=$(PROJECT_NAME) -compact
	$(MAKE) -C $(KROME_BUILD_DIR) gfortran

clean_krome:
	$(MAKE) -C $(KROME_BUILD_DIR) clean

# Only build KROME_OBJS when they do not exist
amrvac.o: | $(KROME_OBJS)
