SETUP_FLAGS := -d=3
TESTS := rtc_3d.log \

rtc_3d.log: rtc_3d.par

include ../../test_rules.make
