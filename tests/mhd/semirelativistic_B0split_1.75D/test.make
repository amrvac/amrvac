SETUP_FLAGS := -d=1 -v=3

TESTS := g72_nosplit.log g72_split_uniform.log g72_split_shear.log

include ../../test_rules.make

# Every mode solves the same problem, and each is checked against its own reference
g72_nosplit.log:       gombosi.par g72.par nosplit.par       short_regression.par
g72_split_uniform.log: gombosi.par g72.par split_uniform.par short_regression.par
g72_split_shear.log:   gombosi.par g72.par split_shear.par   short_regression.par
