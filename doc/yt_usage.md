# AMRVAC datasets and YT

[TOC]

# Introduction {#yt_introduction}

The release of AMRVAC 2.2 included minor changes in the datfile format, in order to improve compatibility with [`yt`](https://yt-project.org),
an open-source Python package for analysing and visualising volumetric data.
`yt` supports AMR meshes and can handle discrete or sampled data such as particles.
It's query-based so that it can handle larger-than-RAM datasets and makes uses of [Cython](https://cython.org) for speed-up and C-interoperability.

The main goal of this webpage is to show off some of yt's functionality and possibilities.
For a detailed guide explaining the full capabilities of `yt` (with examples!),
we refer to [the YT Cookbook](https://yt-project.org/docs/dev/cookbook/index.html).

The `yt` frontend for AMRVAC was developed by Clément Robert and Niels Claes (who are also the current maintainers).
The actual pull request for the frontend, including it in the main `yt` repository, was accepted on the 13th of November, 2019.

Because the frontend is still relatively new, we would greatly appreciate any feedback or issues that you encounter.
The most efficient ways to get in touch with us are either
- [reporting issues on github](https://github.com/yt-project/yt/issues)
- [joining yt's Slackspace](https://yt-project.org/doc/help/index.html#go-on-slack-or-irc-to-ask-a-question) (join channels **help** and **amrvac**)

## Installing yt {#yt_installation}

The [yt website](https://yt-project.org) has detailed instructions on how to install the `yt` package.
As of right now the only way to use the AMRVAC frontend is building [from source](https://yt-project.org/docs/dev/installing.html#installing-yt-from-source).
It's recommended to use `-e .` as well, so you can immediately pull updates from the repository (note the dot after `-e`, this is important).
Required dependencies are `Numpy` and `Cython`.

    pip install numpy cython
    git clone https://github.com/yt-project/yt
    cd yt
    git checkout master
    pip install -e .

## Latest news and updates for the AMRVAC frontend {#latest_news}
- **2019-12-05**: Datasets can now be consistently renormalised to physical units using the `mod_usr.t` AMRVAC normalisations.
- **2019-12-19**: Default derived field list is extended with kinetic energy density, sound speed, mach number and thermal pressure.
- **2020-01-21**: Dust fields are now detected by the frontend; density, momentum and velocity fields are set up together with total dust density and dust to gas ratio.
- **2020-06-23**: yt-4 was merged into master (don't forget to update!).

## Features currently in development/planned {#features_todo}
- Support for staggered grids (_planned_)
- Support for stretched grids (_in development_)
- Support for particle data (_planned_)

## Current limitations & unsupported features {#limitations}
- `yt` supports particle data, but this has not (yet) been implemented in the AMRVAC frontend. This might come later.
- Staggered grids (version 2.2 and later): `yt` will log a warning if staggered grids are loaded, but the flag is currently ignored.
  You can still use all of `yt`'s features on those datasets though.
- Stretched grids are currently unsupported. However, since various users have already told us that they would love to use yt with their
  stretched-grid datasets, we have put this quite high on our to-do list. We are working on this.

# Examples {#examples}
Since adding lots of examples and figures on this webpage would be quite difficult to maintain, we prepared a Jupyter notebook
which contains various examples of using yt with actual MPI-AMRVAC datasets:

[MPI-AMRVAC example notebook on Jupyter nbviewer](https://nbviewer.jupyter.org/github/n-claes/amrvac_yt_docs/blob/master/amrvac_yt_docs.ipynb)

This notebook is regularly updated whenever we update the frontend or have a fun idea that we think is worthwile to add as an example.
If you want to experiment with it, you can either download the notebook directly from nbviewer and use your own data, or you
can clone the repository with the notebook [here](https://github.com/n-claes/amrvac_yt_docs).
