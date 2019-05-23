SETUP_FLAGS := -d=3
TESTS := ltc_sphere_3d.log \

ltc_sphere_3d.log: ltc_3d.par

include ../../test_rules.make
