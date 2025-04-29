.RECIPEPREFIX = >

.PHONY: all

## <-- This part should be removed when migration is complete
src_dir := ../src
vacpp := $(src_dir)/vacpp.pl

sources := $(shell find $(src_dir) -name *.t)
fpp_files := $(sources:../%.t=%.fpp)

all: $(fpp_files)

src/%.fpp: $(src_dir)/%.t
> @mkdir -p $(@D)
> $(vacpp) -d=3 $< > $@
## --> end of temporary part

