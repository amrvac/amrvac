# Frequently Asked Questions

[TOC]

# List of questions {#faq-list}

Below you can find a list of frequently asked questions and answers.

# Is there a mailinglist? {#faq-mailinglist}

Yes, we have a mailinglist that you can send questions to: <mailto:amrvacusers@ls.kuleuven.be>

You can join the mailing-list by [subscribing](https://ls.kuleuven.be/cgi-bin/wa?SUBED1=AMRVACUSERS&A=1), so that you will be informed about important changes or (bug)fixes. You can also search the [mailing list archive](https://ls.kuleuven.be/cgi-bin/wa?A0=AMRVACUSERS).

# I am a user of old MPI-AMRVAC. How can I switch to MPI-AMRVAC 2.0? {#faq-new-version}

MPI-AMRVAC 2.0 is quite different from its previous version. Names of many files, 
subroutines, and parameters are changed. Many subroutines and modules are reorganized.
We suggest following these steps to modify your old application to run with
MPI-AMRVAC 2.0.

1. rename `amrvacusr.t` to `mod_usr.t` and make it a fortran module `module mod_usr`. If you have some
   global variables declared in `amrvacusrpar.t`, declare them in `module mod_usr` and delete `amrvacusrpar.t`.
2. use a tool **amrvac/tools/upgrade.pl** to automatically find the old names and replace them with new ones.
   Just execute `upgrade.pl` in your application folder which contains `mod_usr.t` and `amrvac.par` files.
3. read [instructions](amrvacusr.md) for `mod_usr.t` and modify it accordingly. Pay attention to use mod_hd or mod_mhd at the beginning
   and create subroutine usr_init, in which your old subroutines for initial condition, special boundary 
   condition, special sources,..., need to be pointed to symbolic subroutine pointers collected in _mod_usr_methods.t_,
   and coordinate system and main physics should be specified.
4. If you use special boundary conditions, delete the useless argument 'iw' in `subroutine specialbound_usr`. 
   Special physical sources such as gravity, radiative cooling, thermal conduction, 
   viscosity are now added in a different way and coded in the folder amrvac/src/physics, 
   please read their documents and change user special sources accordingly.
5. Go through [getting started](getting_started.md) to learn the new way to compile and run your application. Note that the `definitions.h` file
   is not used anymore to turn on special functionalities and can be deleted.

Here is an incomplete summary of changes in the `.par` files:

Old parameter | Changes in AMRVAC 2.0
---|---
wnames | Variable names are now automatically defined by the physics modules, so you don't need to specify wnames
dixB | The number of ghost cells is now set automatically, you can override it with nghostcells
fileheadout (etc.) | `base_filename` (in filelist) now defines the base filename for both output and log files
`typeaxial` | Define your geometry with a call to `set_coordinate_system` in your mod_usr.t file, see @ref axial.md
typeadvance | `time_integrator`
typefull1 | `flux_scheme`
typelimiter1 | limiter
fix_small | check_small_values
strictsmall, strictgetaux | Removed
typedivbfix | Has been moved to `mhd_list`, see @ref par.MD

# Variable names in AMRVAC 2.0 {#faq-varnames}

Here are some of the variable names that were changed. Note that e.g. `v1_` en `m1_` were also identical in older version of AMRVAC, as they refer to the same variable (but in primitive/conservative form).

Old name | New name
---|---
`m1_` | mom(1)
`m2_` | mom(2)
`v1_` | mom(1)
`b1_` | mag(1)
`tr1_` | tracer(1)

The variables `e_` (or equivalently, `p_`) and `rho_` have the same name.

# Is there a way I can define my own parameters somewhere in mod_usr.t and configure them through `amrvac.par` ? {#faq-own-parameters}

Indeed, there is a quick and time-saving way to read your own parameters without having to give an explicit value in the usr file and recompile each time. Instead, add this in your usr file : 

1. at the end of the usr_init subroutine, add this line : call params_read(par_files)
2. and just after the usr_init subroutine, define the params_read subroutine you just called. For instance :
  
```{fortran}
    !> Read parameters from a file
    subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /my_list/ my_parameter_1, my_parameter_2 ! where you tell the code to read your own parameters

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, my_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read
```

with `my_parameter_1` and `my_parameter_2` to be defined at the very beginning of the user file, just before the “contains” statement. Doing so, you can use them anywhere in the usr file and they will have the value you defined in your par file adding the following lines :

```{fortran}
&my_list
	my_parameter_1    = 0.2d0
	my_parameter_2    = 1.d0
/
```
# I'm looking for a way to translate some .blk files (using `convert_type= 'onegrid'`) back into standard output `.dat` that I can use to restart a simulation. My goal is to be able to change my data mid-simulation through another program. {#faq-blk-to-dat}

The idea of the `onegrid` option for conversion is that a hierarchical block AMR grid (as stored in the `.dat` files) can be saved to an equivalent uniform grid representation at a user-chosen level (combine `level_io` with `onegrid`). You then can use any software you like to handle uniform grid data. Converting back to the .dat format is then impossible. However, you can write a user-routine to read in the uniform data, and then use it to restart a similar uniform-grid (no AMR, just domain decomposed) simulation. The `.dat` files are all that is needed to do restarts.

If you want to change the data during the simulation, in principle you do not need another program to do that. For that purpose, we provide the generic subroutines

1. usr_process_grid

2. usr_process_global

or immediately after the advance step:

3. usr_process_adv_grid

4. usr_process_adv_global

or just once, at restart (and after reading in a .dat file), then use

5. usr_before_main_loop

which allows you to modify things during runtime (of course, it should make physical sense). See their interface in _mod_usr_methods.t_.

# Is there a way to save (and visualize) the data in the boundary ghost cells? {#faq-show-ghost-cells}

The logical variable `save_physical_boundary` can be set to true, which enforces the `.dat` files to also contain the ghost cell information for physical boundaries (i.e., those beyond the left or right edge of the domain, in each direction). You can use this file (like any other `.dat` file, to restart, and this helps if you want to use saved boundary info in your boundary value handling. However, all our present conversion options (like e.g. to `.vtu` files) do not store this extra info, and you can therefore not use them for visualizing the ghost cell info. For that, you will need to handle the `.dat` files directly, e.g. using python.

# I get an error: mpi.mod was created by a different version of GNU Fortran {#faq-mpi-mod}

This means that your MPI library was compiled with a different version of
`gfortran` (and `GCC`) than you are using now. If you run on a cluster, contact
their support and ask them to fix it.

If you encounter this problem on your own machine, you can try to:

1. Change compilers: Use the version of gfortran that MPI was compiled with
2. Reinstall your MPI library, or install a different one (e.g. MPICH or OpenMPI), which is then hopefully compiled with the right version of GCC/gfortran.
3. Install an operating system with working package management, e.g. Debian ;)

# Can the MHD module handle Hall-MHD? {#faq-hall-mhd}

The MPI-AMRVAC code has Hall-MHD included, as detailed in some of our publications, e.g. in the method paper 

* `MPI-AMRVAC for Solar and Astrophysics', O. Porth, C. Xia, T. Hendrix, S.P. Moschou, & R. Keppens, 2014, ApJS 214,4 [doi:10.1088/0067-0049/214/1/4](http://dx.doi.org/10.1088/0067-0049/214/1/4)

or also in the Kelvin-Helmholtz related application paper 

* `On the influence of environmental parameters on mixing and reconnection caused by the Kelvin-Helmholtz instability at the magnetopause', M.H.J. Leroy & R. Keppens, 2017, PoP 24, 012906 (15pp), [doi:10.1063/1.4974758](http://dx.doi.org/10.1063/1.4974758)

The implementation details are given in the first reference (Porth et al.), and although it works properly on several tests and applications, we note that the time step constraint of our explicit implementation may become prohibitive for particular applications. We just limit the CFL according to \f$\Delta t < \Delta x/ c_w\f$, with time step and spatial step linked by the speed \f$c_w\f$, but in that speed we set \f$c_w= |v|+\max(c_{fast}, \eta_h B/\rho * k_{max})\f$ and \f$k_{max}\f$ is the maximal wavenumber we can represent (i.e. linked to \f$\Delta x\f$). The dispersive nature of the Hall-MHD system may then make \f$\Delta t\f$ going like \f$\Delta x^2\f$, and this limits the current implementation.

In the new version MPI-AMRVAC 2.0, the Hall effect is included when setting the following in the `mhd_list` namelist part

```{fortran}
&mhd_list
      mhd_Hall=.true.
      mhd_etah=  .... (some positive number, quantifying the hall parameter)
/
```

In the same namelist, the optional logical `mhd_4th_order=.true.`
implies a 4th order evaluation for currents, default it is only second order. In any case, you may also need to activate an additional ghost cell layer (or 2 for 4th-order evaluations), through setting appropriately the parameter `nghostcells`.
