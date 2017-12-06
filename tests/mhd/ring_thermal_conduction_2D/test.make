SETUP_FLAGS := -d=2
TESTS := rtc_2d.log \

rtc_2d.log: rtc_2d.par

include ../../test_rules.make
