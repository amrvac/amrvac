# Looking at line-integrated quantities

# Collapsed output

Often one wishes to look at quantities integrated over one particular
direction (collapsing), yielding e.g. column densities. For this purpose
AMRVAC can output line integrals of the simulation variables (plus user-
specified variables) during runtime.
The current implementation is efficient as we take advantage of the grid
structure; however only collapsing along coordinate directions is currently
possible (which is often sufficient). The user can freely specify the
resolution (collapseLevel) of the resulting output data. No re-gridding of the
domain is required, higher levels on the output grid are filled with flat
interpolation.
Note that collapsing of 1D output is not fully implemented. Collapsed output
is sensitive to the switch _saveprim_ of the _filelist_ section, therefore
primitive variable output at runtime is possible. Also the switch _nwauxio_
affects the collapsed view to output user-defined variables.

# Setup in par file

Collapsed views are the fourth file format with output intervals that can be
specified in the _savelist_ section of the par file. To give an example, here
we ask for collapsed arrays in all three coordinate directions, where the
logical _collapse(1)=T_ specifies to output arrays integrated over the first
coordinate direction. The other directions are analoguous.

     &savelist;
            itsave(1,4)      = 0
            dtsave(4)        = 0.1d0
            collapse(1)      = T
            collapse(2)      = T
            collapse(3)      = T
            collapseLevel    = 3
    /

In the standard configuration, AMRVAC will output _*.vti_ files that can be
visualized with paraview. The output filename is composed of the direction and
level. For example, the first collapsed output name reads
_filenameout_d1_l3_nXXXX.vti_ and analoge for the other two.
Note that the order of the (reduced) dimensions in the resulting output files
is preserved, e.g. the first collapsed array in the example above will hold
the y-direction as first coordinate and the z-direction as second coordinate.
This should be kept in mind for visualization.

# Collapsing of existing output

To collapse existing _*.dat_ files, the simulation can be restarted from a
given output time and the code can be brought to a halt after zero iterations.
This is entirely analoguous to the method for [slicing](slices.md) _*.dat_
files and done in the following way: it is best to create a new _*.par_ file
(e.g. collapse.par) and clear the savelist from any output to filetypes other
than _4_. We use itsave to demand a collapse output for the zero-iteration.

     &savelist;
            itsave(1,4)      = 0
            collapse(1)      = T
            collapse(2)      = T
            collapse(3)      = T
            collapseLevel    = 3
    /

The stoplist should look as follows:

     &stoplist;
            reset_it         = .true.
            it_max            = 0
    /

where we reset the iteration counter (so that _itsave(1,4)=0_ will output
collapsed data) and stop the code immediately after the IO (_it_max=0_).

The code can then be started with

    amrvac -i collapse.par -collapsenext 10 -if datamr/data0010.dat

which will take the output _datamr/data0010.dat_ to create new collapsed view 
with index 10 (-collapse 10). The par-file is the newly created collapse.par 
(-i collapse.par) so that the default used to run the code can be left 
untouched. It is a simple exercise in
shell scripting to run along all output-files in one go. For example with the
BASH:

    for i in {0..10}; do printf -v j "%04d" $i; mpirun -np 4 ./amrvac -i collapse.par -collapsenext $i -if datamr/data$j.dat; done

# ASCII output

As an alternative to the cell-type _*.vti_ file format, you can obtain comma-
separated-value _*.csv_ files of the cell-center variables. This can be useful
especially for 2D simulations (ergo 1D line output) which can then be simply
visualized using e.g. gnuplot. For a quick look, the _*.csv_ files can be
visualized with GNUplot using a command similar to the following:

    gnuplot> p 'data_d1_l3_n0000.csv' u 2:1:9 with image

To activate csv ASCII output, the option

     & filelist
            collapse_type    = 'csv'
    /

needs to be set.
