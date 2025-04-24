hash := $(shell printf '%s' "$(sort $(enabled))" | md5sum | cut -c1-8)
build_dir = build/$(target)-$(hash)

all: build/latest/arch.mk

build/latest::
> @echo -e "Symlink $(_cyan)build/latest$(_reset) -> $(_blue)$(build_dir)$(_reset)"
> @mkdir -p $(build_dir)
> @rm -f build/latest
> @ln -s $(target)-$(hash) build/latest

build/latest/arch.mk: | build/latest
> @ln -s $(amrvac)/arch/$(arch).mk build/latest/arch.mk

