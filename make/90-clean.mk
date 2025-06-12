.PHONY: clean clean-all

clean:
> @echo -n "Removing $(build_dir)"
> @rm -rf $(build_dir)
> @rm -f $(build)/latest
> @rm -f amrvac
> @rm -f config.mk

clean-all:
> @echo -n "Removing build dir"
> @rm -rf $(build)
> @rm -f amrvac
> @rm -f config.mk
