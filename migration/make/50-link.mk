amrvac: $(build_dir)/obj/amrvac
> @rm amrvac
> @ln -s $< $@

$(build_dir)/obj/amrvac:
> @echo -e "Linking $(_green)$@$(_reset)"
> @$(link) $(link_flags) $^ -o $@

