.PHONY: build_krome

build_krome:
	cd $(KROME_DIR)
	$(KROME_DIR)/krome -test=hello -compact
	cd -
	$(MAKE) -C build gfortran

amrvac.o: build/krome_all.o build/krome_user_commons.o  build/opkda1.o  build/opkda2.o  build/opkdmain.o
