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


## Standard testcase description

This run will be short, to check that it does not crash and can be finished successfully. I would recommend not modifying this amrvac.par file, instead to work with test_icarus.par file which is also uploaded in this directory. By default this testcase runs a simulation in low resolution, uniform grid, without AMR, with 5 cone CMEs. The simulation lasts 24 days and the output is saved after 14 days.


## Necessary files
The most important files to run Icarus are the

-User file setting up the simulation (mod_usr.t)
-Parameter file that fixes the numerical method, computational domain and Icarus specific settings
-The solar wind boundary (.vtk file)
-The CME parameter file (.in file)

## Additional Scripts
The additional files are uploaded to generate necessary files for Icarus.
-input_boundary_generate.vtk - this python script generates the input VTK file from the standard WSA output file. If you have multiple WSA standard boundary files, it can combine into the timedependent input file. The paths and information should be indicated correctly in the user definition segment of the script.
-convert_output_euhforia_format.py - This python script converts the standard output .csv files which are the satellite to the standard EUHFORIA format with the corresponding satellite names in the filename.



The mangetogram timestamp should be indicated in the amrvac.par file in the icarus_list. When input_boundary_generate.vtk is run, it outputs in the terminal delta_phi value that should be adjusted in &icarus_list for variable delta_phi.



The link to a more detailed documentation will be added soon.
