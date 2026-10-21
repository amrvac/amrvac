TEST_DIRS := Shocks_1.75D/Jiang Shocks_1.75D/Fast Shocks_1.75D/Slow Shocks_1.75D/Intermediate Waves_fine_1.75D/Alfven Waves_fine_1.75D/Fast
#TEST_DIRS := Brio_Wu_1.75D Torrilhon_1.75D Convection_2D

# Directories with altered names (useful for cleaning)
CLEAN_TEST := $(TEST_DIRS:%=clean-%)

.PHONY: all clean $(TEST_DIRS) $(CLEAN_TEST)

all: $(TEST_DIRS)

$(TEST_DIRS):
	@$(MAKE) -C $@ -f test.make

clean:	$(CLEAN_TEST)

$(CLEAN_TEST):
	$(MAKE) -C $(@:clean-%=%) clean
