# Slice output

# Introduction

To alleviate the disk space requirements and overhead of full snapshot output,
it is possible to write hypersurfaces at their own (short) intervals.
Currently, these slices are aligned with the grid, with the benefit that only
a sub-space of the forest has to be traversed to obtain a Morton ordered AMR-
aware subgrid. The slice output is useful especially in expensive 3D
simulations in cartesian geometry, but gives valid output for any dimension
and geometry. The output is composed of the grid cells closest to the
specified subdimensional plane and thus reflects non-interpolated simulation
variables which can be handy for debugging purposes.

## Setup in par file

D-1 dimensional slices are the third file format with output intervals that
can be specified in the _savelist_ section of the par file. To give an
example, here we specify three slices, where the first is perpendicular to the
first coordinate direction _slicedir(1)=1_ and intersects the first axis at
a value of 0.6 _slicecoord(1)=0.6_. The second plane is parametrized
perpendicular to the third coordinate direction _slicedir(2)=3_ and
intersects the third axis at a value of 0.8 _slicecoord(2)=0.8_ and analoge
for the third slice.

     &savelist;
            itsave(1,3)=0
            dtsave(3)=0.1d0
            nslices=3
            slicedir(1)=1
            slicecoord(1)=0.6
            slicedir(2)=3
            slicecoord(2)=0.8
            slicedir(3)=2
            slicecoord(3)=0.7
    /

The total number of slices is specified by _nslices_. 
The data type of the slices is set by _slice_type_ in _filelist_, which by
default is 'vtuCC' in ascii form. if _slice_type='dat'_, the implementation
obtains a properly Morton ordered subdimensional forest with the same levels
as the original simulation, such that the output .dat file can be used for 
restarts in lower dimension. User also has the option 'csv' to obtain 
comma-separated-value *.csv files of the cell-center variables. This can be
useful especially for 2D simulations (1D line output) which can then be 
simply visualized using e.g. gnuplot. For 1D simulations, the code will always 
write a single *.csv file and append the point data together with the output 
time. The file then reads _base_filename-d1-x.600.csv_ in the example above. 
For a quick look, the *.csv files can be imported in Paraview and be 
visualized as points using the filter _Table to Points_. 
The output filename is composed of the direction
and offset values. For example, the first slice output name reads
_filenameout-d1-x.600-nXXXX.vtu_ and analoge for the other two slices.
Note, that the order of the (reduced) dimensions in the resulting output files
is preserved, e.g. the third slice in the example above will hold the
x-direction as first coordinate and the z-direction as second coordinate. 
If you restart from a snapshot and continue the run with slicing, please
choose the next file number _x_ by  _slicenext=x_ in _filelist_ to get continuous 
the file number of slices, otherwise it will start from 0 again overwriting old
files.

## Slicing of existing output

To slice existing *.dat files, the simulation can be restarted from a given
output time and the code can be brought to a halt after zero iterations. This
is done in the following way: its best to create a new *.par file (e.g.
slices.par) and clear the savelist from any output to filetypes other than
_3_. We use itsave to demand a slice output for the zero-iteration.

     &savelist;
            itsave(1,3)=0
            nslices=3
            slicedir(1)=1
            slicecoord(1)=0.6
            slicedir(2)=3
            slicecoord(2)=0.8
            slicedir(3)=2
            slicecoord(3)=0.7
    /

The stoplist should look like the following,

     &stoplist;
            reset_it=.true.
            it_max=0
    /

where we reset the iteration counter (so that _itsave(1,3)=0_ will output
slice data) and stop the code immediately after the IO by set _it_max=0_.

The code can then be started with

    amrvac -i slices.par -slice 10 -if datamr/data0010.dat

which will take the output _datamr/data0010.dat_ 
to create new slices with index 10 (-slice 10). The par-file is
the newly created slices.par (-i slices.par) so that the default used to run
the code can be left untouched. It is a simple exercise in shell scripting to
run along all output-files in one go. For example with the BASH:

    for i in {0..10}; do printf -v j "%04d" $i; ./amrvac -i slices.par -slice $i -if datamr/data$j.dat; done
