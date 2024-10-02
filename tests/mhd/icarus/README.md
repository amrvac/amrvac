02/10/24

Tinatin Baratashvili

## How to run Icarus testcase
----------------

The standard procedure has to be followed to run the testcase of Icarus.

First, retrieve the makefile by
```
setup.pl -d=3
```
and compile it with
```
make -j 4
````

you can run the testcase with
```
mpirun –np 4 ./amrvac -i amrvac.par
```


Note* 4 here is number of CPUs, can be modified depending on the availability of CPUs.


## Necessary files
amrvac.par - parameters for the simulation, remove the following line
```
         typefilelog = 'regression_test'
```
edit the following lines to this in order to perform a run of 24 days and save output every 3h:
```
&savelist
        itsave(1,1)   = 0
        itsave(1,2)   = 0
!        tsavestart(2) = 336.0
        dtsave_dat    = 3.0
/


&stoplist
       time_max      = 576.0d0
/
```
The python script - input_boundary_generate.vtk - generates the input VTK file from the standard WSA output file.


The mangetogram timestamp should be indicated in the amrvac.par file in the icarus_list.


The link to a more detailed documentation will be added soon.
