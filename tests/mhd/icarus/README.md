02/10/24

Tinatin Baratashvili

## How to run Icarus testcase
----------------

To run the Icarus test case, you must follow the standard MPI-AMRVAC procedure. First, you install AMRVAC from https://amrvac.org/md_doc_installation.html. 

In order to get the latest version please checkout amrvac3.2 branch:
```
git checkout amrvac3.2
```
and 
```
git pull
```

Then, make sure you are in ~/icarus folder and retrieve the makefile by
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


- User file setting up the simulation (mod_usr.t)
- Parameter file that fixes the numerical method, computational domain and Icarus specific settings
- The solar wind boundary (.vtk file)
- The CME parameter file (.in file)

# Setting Icarus parameters

Icarus-specific parameter list (&icarus_list) in the `.par` file:

Name | Standard values | Description
---|---|---
`amr_criterion` |`shock`, `tracing` | Default AMR criteria according to Baratashvili et al. 2022
`cme_flag` |0, 1 | 0 means no CME injection, 1 activates CME injection
`num_cmes` | $\ge$ 0 |Number of CMEs to be injected. num_cmes $\le$ number of CMEs in the input `cme_parameter_file`
`relaxation` | 14 | Duration is given in days
`cme_insertion` | 0 | Duration is given in days
`cme_parameter_file` | 'cme_parameters.in' | The file containing CME parameters
`magnetogram_time` |'2015-06-25T01:04:00' | Magnetogram timestamp in 'YYYY-MM-DDTHH:mm:ss' format
`delta_phi` | 1.2 | The correction longitude from WSA file given in radians, output of the VTK file generation python script


# Changing the parameters

Some of the settings that you could change in the `.par` files are:

Name | Description
---|---
`base_filename` | Base file name for output
`tsavestart(2)` | Time to start saving output 3D data [in hours]
`dtsave_dat` | Time between 3D data output [in hours]
`time_max` | Time when simulation ends [in hours]
`time_stepper` | time discretization method, e.g., twostep, threestep, fourstep
`flux_scheme` | spatial discretization method, e.g., tvd, tvdlf, hll
`limiter` | which limiter to use in the spatial discretization, e.g., woodward, minmod, koren
`stretch_dim(1)='uni', stretch_uncentered=.false.` | When uncommented, the grid is stretched radially, when commented it is radially uniform
`refine_max_level` | The maximum number of refinement levels
`block_nx1` | 30 if uniform, 6 if radially stretched
`domain_nx1` | Domain Resolutions: 300 - low, 600 - middle, 1200 - high 
`domain_nx2` | Domain Resolutions: 32 - low, 64 - middle, 126 - high 
`domain_nx3` | Domain Resolutions: 96 - low, 192 - middle, 384 - high 
`xprobmin1, xprobmax1` | Lower and outer boundary values in solar radii 
`omega_frame` | The rotation rate of the Sun [radian per hour]
`boundary_data_file_name' | The input boundary file in the VTK format



A complete list of parameters can be found [par.md](par.md).


## Additional Scripts

The additional files are uploaded to generate necessary files for Icarus.

- input_boundary_generate.vtk - this python script generates the input VTK file from the standard WSA output file. If you have multiple WSA standard boundary files, it can combine into the timedependent input file. The paths and information should be indicated correctly in the user definition segment of the script.
- convert_output_euhforia_format.py - This python script converts the standard output .csv files which are the satellite to the standard EUHFORIA format with the corresponding satellite names in the filename.



The mangetogram timestamp should be indicated in the amrvac.par file in the icarus_list. When input_boundary_generate.vtk is run, it outputs in the terminal delta_phi value that should be adjusted in &icarus_list for variable delta_phi.



The link to a more detailed documentation will be added soon.
