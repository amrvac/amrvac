SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := promRT_25D.log

promRT_25D.log: promRT_25D.par

include ../../test_rules.make
