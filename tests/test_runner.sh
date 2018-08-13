# this is part of the travis ci/cd pipeline
# call syntax:
# bash test_runner.sh hd
#
# Run tests in $1 and check if any fails
# return appropriate exit code (0 if and only if all tests pass)

make -s $1 | tee $1.out
! grep "**FAILED**" $1.out
