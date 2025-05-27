arch := llvm

compile = flang
f90_flags += -ffree-form -fimplicit-none -cpp $(shell mpifort --showme:compile)
link_flags += $(f90_flags) $(shell mpifort --showme:link)

