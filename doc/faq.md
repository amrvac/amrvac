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

## How does the git-version differ from the previous svn-version?

We have compiled a [list of changes](gitversion.md). Other than that, we try to
keep the [documentation](contents.md) up to date.

## I just did a 'git pull' and now my setup won't compile anymore, what to do?

It is possible that we needed to change the `makefile`. Since the `makefile` is
in your working folder (see below), it does not match the new requirements and
compilation might fail because of that. In your working folder do

    $AMRVAC_DIR/setup.pl -s

and copy the output. Now remove the old makefile

    rm makefile

paste and execute the setup line. This gets the new makefile from the repository
and configures it accordingly. Now try to

    make clean && make

## Corrupt snapshot files?

**Q:** In some cases, the .dat file is screwed up, why is this happening? like,
ndim is a horribly big number, or the other variables at the end of the .dat
file so I cannot use them for restart or convert them post-simulation this
happens apparently random, I mean, with the same compiler, nr. of processors,
etc. (due to Norbert Magyar via Skype...)

**A:** Most likely, you are using typeparIO=1 which uses parallel MPI-IO for
snapshot read/write. Unfortunately this can fail sporadically on some clusters.
Try using the other options typeparIO=0 (serial MPI-IO) or typeparIO=-1 (serial
FortranIO).

## How should I setup a new simulation folder?

The philosophy is now to keep all the files needed for the simulation together
in one *working folder*. These files are typically the following:

* `amrvacusr.t`
* `amrvacusrpar.t`
* `amrvac.par` (any name is possible)
* `definitions.h` (optional)
* `makefile` (usually not needed)

It is good style to also inclue a README file with some explanations as well as
the options for `$AMRVAC_DIR/setup.pl`. These options are now saved in the
**makefile** however.

You can create this folder anywhere on your machine and type `make` there since
`$AMRVAC_DIR` points to the remaining source files. Note also that all the `*.t`
files in your working folder will get precedence over the ones in the `src/*`.
So if you need to make a modification to any of those, you can copy them in,
modify and not contaminate your repository in the process.

## What is the localmakefile for?

You can add additional source files to your working folder. These are found by
make through the **localmakefile**. For example:

    #########################################################################
    # Here we have some auxilary modules
    #########################################################################
    FOBJECTS += init_particles$F  rand3$F
    OBJECTS += init_particles.o   rand3.o

    #########################################################################
    # Dependencies:
    #########################################################################
    init_particles.o : mod_particles.o rand3.o
    advance.o        : mod_particles.o
    amrvacusr.o      : mod_odeint.o mod_constants.o

This migth be useful as the **amrvacusr.t** can become quite unwieldy.

## What is then the purpose of the src/usr folder?

The src/usr folder is **deprecated** and no longer used.

## Compilation of file ... is very slow, why?

It seems that since some version of ifort 13 it takes much longer to compile
several files. We have observed that for example for the advance.f. The 12…
branch is much faster in compilation time and you could switch to that if you
have it. Other than that, you could lower the optimization, e.g. setup with
`-arch=debug` or `-arch=gfortrandebug` until you make a production setup with
`-arch=default`.

