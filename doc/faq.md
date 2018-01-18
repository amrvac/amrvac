# Frequently Asked Questions

# List of questions

Below you can find a list of frequently asked questions and answers.

## Is there a mailinglist?

Yes, we now have a mailinglist, which is used to:

* Keep you updated on important (bug)fixes
* Announce new features or changes
* Post general questions

You can join the mailing-list
by [subscribing](https://ls.kuleuven.be/cgi-bin/wa?SUBED1=AMRVACUSERS&A=1). Then
you will be able to send and receive mails
from <mailto:amrvacusers@ls.kuleuven.be>.

## I am a user of old MPI-AMRVAC. How can I switch to MPI-AMRVAC 2.0?

MPI-AMRVAC 2.0 is quite different from its previous version. Names of many files, 
subroutines, and parameters are changed. Many subroutines and modules are reorganized.
We suggest following these steps to modify your old application to run with
MPI-AMRVAC 2.0.

1. rename amrvacusr.t to mod_usr.t and make it a fortran module 'module mod_usr'. If you have some
   global variables declared in amrvacusrpar.t, declare them in 'module mod_usr' and delete amrvacusrpar.t.
2. use a tool amrvac/tools/upgrade.pl to automatically find the old names and replace them with new ones.
   Just execute upgrade.pl in your application folder which contains mod_usr.t and amrvac.par files.
3. read [instructions](amrvacusr.md) for mod_usr.t and modify it accordingly. Pay attention to use mod_hd or mod_mhd at the beginning
   and create subroutine usr_init, in which your old subroutines for initial condition, special boundary 
   condition, special sources,..., need to be pointed to symbolic subroutine pointers collected in mod_usr_methods.t,
   and coordinate system and main physics should be specified.
4. If you use special boundary conditions, delete the useless argument 'iw' in subroutine specialbound_usr. 
   Special physical sources such as gravity, radiative cooling, thermal conduction, 
   viscosity are now added in a different way and coded in the folder amrvac/src/physics, 
   please read their documents and change user special sources accordingly.
5. Go through [getting started](getting_started.md) to learn the new way to compile and run your application. Note that definitions.h file
   is not used anymore to turn on special functionalities and can be deleted.
