# 3D Printing

This is a short tutorial that shows how to prepare MPI-AMRVAC output for 3D printing.
The advice in this tutorial is of a general nature. Please contact your
printing service for specifics.

Requirements: A 3D printer needs a description of a surface area in three
dimensions in the form of an .stl file.

Procedure:

  1. Write your output in .vtu format and load it into the VisIt visualization program. This software package is free and can be obtained from: <https://wci.llnl.gov/simulation/computer-codes/visit/>
  2. In visit, create a surface area that represents the result of your simulation. This can be done through:
    * Add --&gt; Pseudocolor
    * OpAtts --&gt; Slicing --&gt; isosurface
    * Make sure to draw only one isosurface. The default is ten, which will not yield a printable result.
    * Scale the result to a size that fits in the printer. This can be done through OpAtts --&gt; Transforms --&gt; Transform --&gt; Scale   Keep in mind that the cost of a 3D print scales with the volume.
  3. Once you are satisfied with the result, export your surface through
    * File --&gt; Export database
    * Choose STL as the file format and use the option to write binary
    * Keep in mind that each plot will generate its own database. If you have multiple plots in one window, VisIt will write a separate file for each.
  4. The STL file produced by VisIt will not be directly printable. It will contain flaws that need to be repaired.
    * This can be done online at either:  <http://www.netfabb.com/netfabbcloud.php> or <https://cloud.materialise.com/>.
    * Alternatively, you can use a program such as Meshlab, available for free at [http://meshlab.sourceforge.net/](http://meshlab.sourceforge.net/). This will also let you  manipulate the STL file by hand, allowing you to merge multiple STL files or insert additional objects. Once you are done, it would be wise to apply one of the online services to your end result in order to fix any errors that you may have missed.
  5. The process demonstrated in the previous steps will yield a printable file. However, it may still cause problems. For example, if part of your surface is not connected to the rest. If that happens, the printer will create two separate objects, which may not be desirable. It is up to the user to make sure to check this in advance. Also, you may run into trouble if the size of the .stl file is too large. Generally, it is a good idea to keep it below 100 MB. This may require that you use a low resolution version of the .vtu file. This can be specified with the 'level_io_max' parameter in the .par file, before converting to .vtu.

