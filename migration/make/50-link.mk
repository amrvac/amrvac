# TODO: link to amrvac.a for later linkage with mod_usr.o

amrvac: $(build_dir)/obj/amrvac
> @rm -f amrvac
> @ln -s $< $@

$(build_dir)/obj/amrvac:
> @echo -e "Linking $(_green)$(notdir $@)$(_reset)"
> @$(link) $(link_flags) $^ -o $@

