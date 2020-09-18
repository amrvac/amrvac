# this is part of the continuous integration pipeline
# call syntax:
# bash test_runner.sh <name>
#
# Run tests in $1 and check if any fails
# return appropriate exit code (0 if and only if all tests pass)

make -j 2 -s $1 | tee $1.out; test ${PIPESTATUS[0]} -eq 0
cond1=$? # first condition: build correctly with make

if [ $cond1 != 0 ] ; then
    echo "BUILD ERROR: make $1 failed to build the executable"
    exit 1
fi

! grep "**FAILED**" $1.out
cond2=$? # second condition: no error detected while running tests

if [ $cond2 != 0 ] ; then
    echo "TEST ERROR: at least one $1 test failed"
    exit 1
fi
