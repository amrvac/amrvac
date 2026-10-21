TEST_DIRS := Advect_BC_from_file_2D AlfvenShock_MHD2D Shu_Osher_1D_HD

# Directories with altered names (useful for cleaning)
CLEAN_TEST := $(TEST_DIRS:%=clean-%)

.PHONY: all clean $(TEST_DIRS) $(CLEAN_TEST)

all: $(TEST_DIRS)

$(TEST_DIRS):
	@$(MAKE) -C $@ -f test.make

clean:	$(CLEAN_TEST)

$(CLEAN_TEST):
	$(MAKE) -C $(@:clean-%=%) clean
