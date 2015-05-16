; $Id: slicer3.pro,v 1.18 1998/03/26 05:20:41 alan Exp $
;
; Copyright (c) 1993-1998, Research Systems, Inc.  All rights reserved.
;	Unauthorized reproduction prohibited.
;
;+
; NAME:
;       SLICER3
;
; PURPOSE:
;       Widget based application to visualize 3D data.
;       This program superseeds the "SLICER" program.
;
; CATEGORY:
;       Volume display / rendering.
;
; CALLING SEQUENCE:
;
;       SLICER3 [, hData3D]
;
; INPUTS:
;       hData3D:    A pointer to the 3D data, or an array of pointers
;                   to multiple 3D arrays.   If multiple arrays are specified,
;                   they all must have the same X, Y, and Z dimensions.
;                   This parameter is optional.   The default is to use a 3D
;                   array created from BYTARR(2,2,2).   While running SLICER3,
;                   the user may interactively load data via the file menu
;                   (see example).   If data is loaded in this fashion,
;                   any data passed to SLICER3 via a pointer (or pointers) is
;                   deleted, and the pointers become invalid.
;
; KEYWORD PARAMETERS:
;       DETACH:     If set, then the drawing area is placed in a base that is
;                   detached from the control panel.   The drawing area can
;                   only be detached if Slicer3 is not run in modal mode.
;       MODAL:      If set, then Slicer3 will block user interaction with all
;                   other widgets (and block the command line) until the user
;                   quits Slicer3.   If Slicer3 is started from some other
;                   widget-based application, then it is usually advisable
;                   to run Slicer3 in Modal mode.
;       GROUP:      This keyword specifies a widget ID of the group leader.
;                   If the specified widget is destroyed, Slicer3 is also
;                   destroyed.   If Slicer3 is started from a widget
;                   application, then GROUP should ALWAYS be specified.
;                   See example.
;       DATA_NAMES: A string array of names for the data. The names appear
;                   on the droplist widget for the current data. If the
;                   number of elements of DATA_NAMES is less than the
;                   number of elements in hData3D then default names will
;                   be generated for the unnamed data.
;
; COMMON BLOCKS:
;       COMMON colors, r, g, b, cur_red, cur_green, cur_blue
;                   These common variables are used by the "STRETCH",
;                   "LOADCT", and "XLOADCT" commands.
;
; SIDE EFFECTS:
;                   Slicer3 modifies the current color table, as well as
;                   various elements of the plotting system (ie, the "!X",
;                   "!Y", "!Z", and "!P" system variables).
;                   If the "MODAL" keyword is set (usually a good idea),
;                   then SLICER3 will, upon exit, restore these system
;                   variables (and the color tables) to the values they
;                   had when SLICER3 was started.
;                   Slicer3 sets the position for the light source and
;                   enables back-facing polygons to be drawn (see the
;                   IDL "SET_SHADING" command).
;
;                   Slicer3 overwrites the existing contents of the
;                   Z-buffer.   Upon exiting Slicer3, the Z-buffer contents
;                   are the same as what was last displayed by Slicer3.
;
;                   On 24-bit displays, Slicer3 sets the device to
;                   non-decomposed color mode (DEVICE, DECOMPOSED=0).
;
;                   Slicer3 breaks the color table into 6 "bands", based upon
;                   the number of available colors (max_color=!D.N_COLORS
;                   on 8-bit displays, and max_color=256 on 24-bit displays) :
;
;                      Band start index:    Band end index:    Used for:
;                      -----------------    ---------------    ---------
;
;                      0                    nColor-1           X Slices.
;                      nColor               (2*nColor)-1       Y Slices.
;                      2*nColor             (3*nColor)-1       Z Slices.
;                      3*nColor             (4*nColor)-1       Iso-surfaces.
;                      4*nColor             (5*nColor)-1       Projections.
;
;                      Where:
;                              nColor = (max_color - 9) / 5
;
;                      Note that the value of !D.N_Colors can vary from
;                      machine to machine, and from run to run, depending
;                      upon available system resources.   Also, !D.N_Colors
;                      is usually not set by IDL until the first window has
;                      been created (or realized) in that IDL session.
;
;                   Annotation colors are the last "band", and they are
;                   set up as :
;
;                      Color index:     Color:
;                      -------------    ------
;
;                      max_color - 1    White.
;                      max_color - 2    Yellow.
;                      max_color - 3    Cyan.
;                      max_color - 4    Purple.
;                      max_color - 5    Red.
;                      max_color - 6    Green.
;                      max_color - 7    Blue.
;                      max_color - 8    Black.
;
;                   On 24-bit displays, improved performance can often be
;                   gained by running Slicer3 in 8-bit mode.   This can be
;                   accomplished (on some platforms) by entering the following
;                   command at the start of the IDL session (before any
;                   windows are created):
;
;                      Device, Pseudo_Color=8
;
;                   See the documentation for additional information.
;
; RESTRICTIONS:
;       The data used by Slicer3 must meet the following conditions:
;          * The data must have three dimensions.
;          * The minimum size of the data array must be 2x2x2.
;          * If multiple volumes are loaded, they all must have the
;            same dimensions.
;
; PROCEDURE:
;
;       "File" menu:
;
;          "Load":
;                  Select a file containing a 3D array (or arrays) to load
;                  into Slicer3.   The file must have been written in a
;                  certain binary format.   For each data array in the file,
;                  the following values are present:
;
;                     data item                     data type     bytes
;                     --------------------------    ----------    ------
;
;                     Number of dimensions          long          4
;                     in array.   Note that
;                     this is always 3 for
;                     valid Slicer3 data.
;
;                     Size of first dimension.      long          4
;                     Size of second dimension.     long          4
;                     Size of third dimension.      long          4
;
;                        If multiple arrays are present in the file,
;                        they must all have the same dimensions.
;
;                     Data type (1 through 5)       long          4
;                     (see the IDL "SIZE"
;                     function for types).
;
;                     Total number of elements.     long          4
;                     (dimX*dimY*dimZ).
;
;                        Note that the all of the above values are the
;                        exact output of the IDL "SIZE" function.
;
;                     Number of characters          long          4
;                     in data name.
;
;                        Note that the above value is the output from
;                        the IDL "STRLEN" function.
;
;                     Data name.                    byte          strlen()
;
;                     3D data.                      varies        varies
;
;                     Note that the 3D data type and number of bytes
;                     is specified by the "size" information above.
;
;                  Any number of 3D datasets can be concatenated into
;                  a single file of this type (as long as they all have
;                  the same dimensions).
;
;                  (See EXAMPLE, below.)
;
;                  NOTE: Files saved by the "Save Subset" operation
;                  (see below) are suitable for input via the "Load"
;                  operation.
;
;                  Data files that are moved from one platform to
;                  another may not load as expected, due to differing
;                  byte order.   See the "BYTEORDER" and "SWAP_ENDIAN"
;                  IDL commands for details.
;
;          "Save / Save Subset":
;                  Slicer3 must be in "Block" mode for this operation to be
;                  available.   When selected, a subset of the 3D data
;                  enclosed in the current block is written to the chosen
;                  save file.   This subset can then be loaded back into
;                  Slicer3 at any time.   If multiple 3D arrays are
;                  currently available in Slicer3, then multiple subsets
;                  are saved to the file.
;
;          "Save / Save Tiff Image":
;                  When selected, a tiff image of the current Slicer3
;                  contents is saved to the chosen file.   When running in
;                  8-bit mode, a "Class P" palette color Tiff file is created.
;                  In 24-bit mode, a "Class R" (interleaved by image) Tiff
;                  file is created.
;
;          "Quit":
;                  Exits Slicer3.
;
;
;       "Tools" menu:
;
;          "Erase":
;                  Erases the display window and deletes all the objects
;                  in the display list.
;
;          "Delete / ...":
;                  As graphical objects are created, they are added to the
;                  display list.   The "Delete" menu allows the user to
;                  delete a specific object from the list.   When an object
;                  is deleted, the screen is redrawn with the remaining
;                  objects.
;
;          "Colors / Reset Colors":
;                  Selecting this will cause the original color scheme to
;                  be restored.
;
;          "Colors / Differential Shading":
;                  This allows the user to change the percentage of
;                  differential shading applied to the X, Y, and Z slices.
;
;          "Colors / Slice/Block":
;                  This allows the user to use the "XLOADCT" operation
;                  to modify the colors used for slices and blocks.
;                  In some cases, the new colors will not be visible
;                  until the user selects "Done" in the XLOADCT tool.
;
;          "Colors / Surface":
;                  This allows the user to use the "XLOADCT" operation
;                  to modify the colors used for iso-surfaces.
;
;          "Colors / Projection":
;                  This allows the user to use the "XLOADCT" operation
;                  to modify the colors used for projections.
;
;          Note that on some platforms, the selected colors may not
;          become visible until after the "XLOADCT" tool is exited.
;
;          "Options":
;                  This brings up a panel allowing the user to set:
;                     The axis visibility.
;                     The wire-frame cube visibility.
;                     The display window size
;                        (the X and Y dimensions are always the same).
;                  If the user selects "Ok", then the display is redrawn.
;
;
;       "About" menu:
;
;          "About Slicer":
;                  Brings up help information about Slicer3.
;
;
;       "Data:" pull-down menu:
;                  If multiple datasets are currently available in Slicer3,
;                  this menu allows the selection of the current data.
;                  Slices, blocks, iso-surfaces, etc. are created from
;                  the currently selected data.   If only one dataset
;                  is currently loaded, then this menu is inactive.
;
;
;       "Mode:" pull-down menu:
;                  This menu is used to select the current mode of operation.
;
;
;       Main Draw Window:
;                  Interaction in the main draw window is dependent upon
;                  the currently selected mode ("Slice", "Block", "Surface",
;                  etc., see below).   In general, when coordinate input is
;                  required from the user, it is performed by clicking a
;                  mouse button on the "surface" of the wire-frame cube that
;                  surrounds the data.   This 3D location is then used as
;                  the basis for whatever input is needed.   In most cases,
;                  the "front" side of the cube is used.   In a few cases,
;                  the coordinate input is on the "back" side of the cube.
;
;
;       "Slice" mode:
;                  To display a slice, click and drag the left mouse button
;                  on the wire-frame cube.   When the button is released, a
;                  slice through the data will be drawn at that location.
;
;          "Draw" mode:
;                  When in Draw mode, new slices will be merged into
;                  the current Z-buffer contents.
;
;          "Expose" mode:
;                  When in Expose mode, new slices will be drawn in
;                  front of everything else.
;
;          "Orthogonal" mode:
;                  When in Orthogonal mode, use the left mouse button
;                  in the big window to position and draw an orthogonal
;                  slicing plane.   Clicking the right mouse button in
;                  the big window (or any mouse button in the small
;                  window) will toggle the slicing plane orientation.
;             "X":
;                  This sets the orthogonal slicing plane orientation
;                  to be perpendicular to the X axis.
;             "Y":
;                  This sets the orthogonal slicing plane orientation
;                  to be perpendicular to the Y axis.
;             "Z":
;                  This sets the orthogonal slicing plane orientation
;                  to be perpendicular to the Z axis.
;
;          "Oblique" mode:
;                  Clicking any mouse button in the small window will
;                  reset the oblique slicing plane to its default
;                  orientation.
;             "Normal" mode:
;                  When in this mode, click and drag the left mouse
;                  button in the big window to set the surface normal
;                  for the oblique slicing plane.
;             "Center" mode:
;                  When in this mode, click and drag the left mouse
;                  button in the big window to set the center point
;                  for the surface normal.
;             "Display":
;                  Clicking this button will cause an oblique slicing
;                  plane to be drawn.
;
;
;       "Block" mode:
;                  When in Block mode, use the left mouse button in the
;                  big window to set the location for the "purple" corner
;                  of the block.   Use the right mouse button to locate
;                  the opposite "blue" corner of the block.
;
;                  When in Block mode, the "Save Subset" operation under
;                  the main "File" menu is available.
;
;             "Add" mode:
;                  When in this mode, the block will be "added" to the
;                  current Z-buffer contents.
;
;             "Subtract" mode:
;                  When in this mode, the block will be "subtracted"
;                  from the current Z-buffer contents.   Subtract mode
;                  is only effective when the block intersects some
;                  other object in the display (such as an iso-surface).
;
;             "Display":
;                  Clicking this button will cause the block to be drawn.
;
;
;        "Surface" mode:
;                  In iso-surface is like a contour line on a contour
;                  map.   On one side of the line, the elevation is higher
;                  than the contour level, and on the other side of the
;                  line, the elevation is lower than the contour level.
;                  An iso-surface, however, is a 3D surface that passes
;                  through the data such that the data values on one side
;                  of the surface are higher than the threshold value,
;                  and on the other side of the surface, the data values
;                  are lower than the threshold value.
;
;                  When in Surface mode, a logarithmic histogram plot
;                  of the data is displayed in the small draw window.
;                  Click and drag a mouse button on this plot to set
;                  the iso-surface threshold value.   This value is
;                  also shown in the text widget below the plot.
;                  The threshold value may also be set by typing a
;                  new value in this text widget.   The histogram
;                  plot is affected by the current threshold settings.
;                  (See Threshold mode, below).
;
;             "Low":
;                  Selecting this mode will cause the iso-surface polygon
;                  facing to face towards the lower data values.
;                  Usually, this is the mode to use when the iso-surface
;                  is desired to surround high data values.
;
;             "High":
;                  Selecting this mode will cause the iso-surface polygon
;                  facing to face towards the higher data values.
;                  Usually, this is the mode to use when the iso-surface
;                  is desired to surround low data values.
;
;             "Shading" pull-down menu:
;                  Iso-surfaces are normally rendered with light-source
;                  shading.   If multiple datasets are currently loaded,
;                  then this menu allows the selection of a different
;                  3D array for the source of the iso-surface shading
;                  values.   If only one dataset is currently loaded,
;                  then this menu is inactive.
;
;             "Display":
;                  Clicking this button will cause the iso-surface to
;                  be created and drawn.   Iso-surfaces often consist
;                  of tens of thousands of polygons, and can sometimes
;                  take considerable time to create and render.
;
;
;        "Projection" mode:
;                  A "voxel" projection of a 3D array is the projection
;                  of the data values within that array onto a viewing
;                  plane.   This is similar to taking an X-ray image of
;                  a 3D object.
;
;             "Max" mode:
;                  Select this mode for a Maximum intensity projection.
;
;             "Avg" mode:
;                  Select this mode for an Average intensity projection.
;
;             "Low" mode:
;                  Select this mode for a Low resolution projection.
;
;             "Med" mode:
;                  Select this mode for a Medium resolution projection.
;
;             "High" mode:
;                  Select this mode for a High resolution projection.
;
;             "Depth Queue %":
;                  Use the slider to set the depth queue percent.
;                  A value of 50, for example, indicates that the
;                  farthest part of the projection will be 50 % as
;                  bright as the closest part of the projection.
;
;             "Display":
;                  Clicking this button will cause the projection to
;                  be calculated and drawn.   Projections can sometimes
;                  take considerable time to display.   Higher resolution
;                  projections take more computation time.
;
;
;        "Threshold" mode:
;                  When in Threshold mode, a logarithmic histogram plot
;                  of the data is displayed in the small draw window.
;                  Click and drag the left mouse button on this plot to
;                  set the minimum and maximum threshold values.
;                  To expand a narrow range of data values into the
;                  full range of available colors, set the threshold
;                  range before displaying slices, blocks, or projections.
;                  The threshold settings also affect the histogram
;                  plot in "Surface" mode.   The minimum and maximum
;                  threshold values are also shown in the text widgets
;                  below the histogram plot.
;
;                  Click and drag the right mouse button on the histogram
;                  plot to set the transparency threshold.
;                  Portions of any slice, block, or projection that are
;                  less than the transparency value are not drawn (clear).
;                  Iso-surfaces are not affected by the transparency
;                  threshold.   The transparency threshold value is also
;                  shown in a text widget below the histogram plot.
;
;           "Min":
;                  In this text widget, a minimum threshold value can
;                  be entered.
;
;           "Max":
;                  In this text widget, a maximum threshold value can
;                  be entered.
;
;           "Transp.":
;                  In this text widget, a transparency threshold value
;                  can be entered.
;
;
;        "Profile" mode:
;                  In Profile mode, a plot is displayed showing the
;                  data values along a line.   This line is also shown
;                  superimposed on the data in the main draw window.
;                  The bottom of the plot corresponds to the "purple"
;                  end of the line, and the top of the plot corresponds
;                  to the "blue" end of the line.
;
;           "Orthogonal" mode:
;                  Click and drag the left mouse button to position the
;                  profile line, based upon a point on the "front"
;                  faces of the wire-frame cube.   Click and drag the
;                  right mouse button to position the profile line,
;                  based upon a point on the "back" faces of the
;                  wire-frame cube.   As the profile line is moved,
;                  The profile plot is dynamically updated.
;
;           "Oblique" mode:
;                  Click and drag the left mouse button to position the
;                  "purple" end of the profile line on one of the "front"
;                  faces of the wire-frame cube.   Click and drag the
;                  right mouse button to position the "blue" end of the
;                  profile line on one of the "back" faces of the
;                  wire-frame cube.   As the profile line is moved,
;                  The profile plot is dynamically updated.
;
;
;        "Probe" mode:
;                  In Probe mode, click and drag a mouse button over
;                  an object in the main draw window.   The actual
;                  X-Y-Z location within the data volume is displayed
;                  in the three text widgets.   Also, the data value
;                  at that 3D location is displayed in the status
;                  window, above the main draw window.   If the cursor
;                  is inside the wire-frame cube, but not on any object,
;                  then the status window displays "No data value", and
;                  the three text widgets are empty.  If the cursor is
;                  outside the wire-frame cube, then the status window
;                  and text widgets are empty.
;
;           "X":
;                  Use this text widget to enter the X coordinate for
;                  the probe.
;
;           "Y":
;                  Use this text widget to enter the Y coordinate for
;                  the probe.
;
;           "Z":
;                  Use this text widget to enter the Z coordinate for
;                  the probe.
;
;
;        "View" mode:
;                  In view mode, a small window shows the orientation
;                  of the data cube in the current view.   As view
;                  parameters are changed, this window is dynamically
;                  updated.   The main draw window is then updated
;                  when the user clicks on "Display", or exits View
;                  mode.
;
;        "Display":
;                  Clicking on this button will cause the objects in
;                  the main view window to be drawn in the new view.
;                  If any view parameters have been changed since
;                  the last time the main view was updated, the main
;                  view will be automatically redrawn when the user
;                  exits View mode.
;
;        1st Rotation:
;                  Use this slider to set the angle of the first view
;                  rotation (in degrees).   The droplist widget adjacent
;                  to the slider indicates which axis this rotation is
;                  about.
;
;        2nd Rotation:
;                  Use this slider to set the angle of the second view
;                  rotation (in degrees).   The droplist widget adjacent
;                  to the slider indicates which axis this rotation is
;                  about.
;
;        "Zoom %":
;                  Use this slider to set the zoom factor percent.
;                  Depending upon the view rotations, Slicer3 may
;                  override this setting to ensure that all eight
;                  corners of the data cube are within the window.
;
;        "Z %":
;                  Use this slider to set a scale factor for the Z
;                  axis (to compensate for the data's aspect ratio).
;
;
; EXAMPLE:
;       Example 1:
;       ----------
;       Create a data save file suitable for dynamic loading into
;       Slicer3.
;
;
;               ; Store some 3D data in a variable called "data_1".
;               data_1 = INDGEN(20,30,40)
;
;               ; Store some 3D data in a variable called "data_2".
;               data_2 = FINDGEN(20,30,40)
;
;               ; Define the names for the datasets (their names will
;               ; appear in the "Data:" pull-down menu in Slicer3.
;
;               data_1_name = 'Test Data 1'
;               data_2_name = 'Data 2'
;
;               ; Select a data file name.
;               dataFile = PICKFILE()
;
;               ; Write the file.
;
;               GET_LUN, lun
;               OPENW, lun, dataFile
;
;               WRITEU, lun, SIZE(data_1)
;               WRITEU, lun, STRLEN(data_1_name)
;               WRITEU, lun, BYTE(data_1_name)
;               WRITEU, lun, data_1
;
;               WRITEU, lun, SIZE(data_2)
;               WRITEU, lun, STRLEN(data_2_name)
;               WRITEU, lun, BYTE(data_2_name)
;               WRITEU, lun, data_2
;
;               CLOSE, lun
;               FREE_LUN, lun
;
;
;       Example 2:
;       ----------
;       Run Slicer3 with data passed to it at startup.
;
;               ; Create some 3D data.
;               data = INDGEN(20,30,40)
;
;               ; Create a pointer to the data, and use the "/NO_COPY"
;               ; keyword to save memory.
;               h_data = PTR_NEW(data, /NO_COPY)
;
;               ; Start up Slicer3.
;               SLICER3, h_data, /MODAL
;
;               ; If the user did not interactively load any data into
;               ; Slicer3 (via the "File/Load" menu), then the original
;               ; pointer to the data still exists (and the original data
;               ; will still reside in memory).   To free it, use:
;
;               if PTR_VALID(h_data) then PTR_FREE, h_data
;
;               ; If the pointer is no longer valid, then that indicates
;               ; that the user interactively loaded data into Slicer3.
;               ; Any data that is loaded interactively is automatically
;               ; deleted when the user exits Slicer3.
;
;               ; Note that the last contents of the main view window in
;               ; Slicer3 still resides in the Z-buffer.   To access this
;               ; image after exiting Slicer3, perform the following actions:
;
;               current_device = !D.Name
;               SET_PLOT, 'Z'
;               image_buffer = TVRD()
;               depth_buffer = TVRD(CHANNEL=1, /WORDS)
;               SET_PLOT, current_device
;               TV, image_buffer
;
;               ; Note that the image contained in "image_buffer" will look
;               ; "correct" only if the colors loaded by Slicer3 have not
;               ; been changed since the user exited Slicer3.
;
;
; MODIFICATION HISTORY:
;       Daniel Carr - RSI, Fri Nov 22 15:43:36 MST 1996
;       Daniel Carr - RSI, Fri Jan 10 12:08:01 MST 1997
;          Fixed bugs and added muti-dataset capability.
;       Alan Youngblood, Daniel Carr - RSI, Wed Feb 11 10:07:32 MST 1998
;          Modified routine to use pointers.
;
;
;-

;******************************************************************************
; Procedure to reset IDL system variables to their initial start-up states.
pro Viz3D_Reset

   !Order = 0

   T3d, /Reset
   !P.T3d = 0
   !P.Position = [0.0, 0.0, 0.0, 0.0]
   !P.Clip = [0L, 0L, (!D.X_Size-1L), (!D.Y_Size-1L), 0L, 0L]
   !P.Region = [0.0, 0.0, 0.0, 0.0]
   !P.Background = 0L
   !P.Charsize = 0.0
   !P.Charthick = 0.0
   !P.Color = 255
   !P.Font = (-1L)
   !P.Linestyle = 0L
   !P.Multi = [0L, 0L, 0L, 0L, 0L]
   !P.Noclip = 0L
   !P.Noerase = 0L
   !P.Nsum = 0L
   !P.Psym = 0L
   !P.Subtitle = ''
   !P.Symsize = 0.0
   !P.Thick = 0.0
   !P.Title = ''
   !P.Ticklen = 0.02
   !P.Channel = 0

   !X.S = [0.0, 0.0]
   !X.Style = 0L
   !X.Range = [0.0, 0.0]
   !X.Margin = [10.0, 3.0]
   !X.Type = 0L
   !X.Ticks = 0L
   !X.Ticklen = 0.0
   !X.Thick = 0.0
   !X.Crange = [0.0, 0.0]
   !X.Omargin = [0.0, 0.0]
   !X.Window = [0.0, 0.0]
   !X.Region = [0.0, 0.0]
   !X.Charsize = 0.0
   !X.Minor = 0L
   !X.Tickv = Replicate(0.0, 30)
   !X.Tickname = Replicate('', 30)
   !X.Gridstyle = 0L
   !X.Tickformat = ''
   !X.Title = ''

   !Y.S = [0.0, 0.0]
   !Y.Style = 0L
   !Y.Range = [0.0, 0.0]
   !Y.Margin = [4.0, 2.0]
   !Y.Type = 0L
   !Y.Ticks = 0L
   !Y.Ticklen = 0.0
   !Y.Thick = 0.0
   !Y.Crange = [0.0, 0.0]
   !Y.Omargin = [0.0, 0.0]
   !Y.Window = [0.0, 0.0]
   !Y.Region = [0.0, 0.0]
   !Y.Charsize = 0.0
   !Y.Minor = 0L
   !Y.Tickv = Replicate(0.0, 30)
   !Y.Tickname = Replicate('', 30)
   !Y.Gridstyle = 0L
   !Y.Tickformat = ''
   !Y.Title = ''

   !Z.S = [0.0, 0.0]
   !Z.Style = 0L
   !Z.Range = [0.0, 0.0]
   !Z.Margin = [0.0, 0.0]
   !Z.Type = 0L
   !Z.Ticks = 0L
   !Z.Ticklen = 0.0
   !Z.Thick = 0.0
   !Z.Crange = [0.0, 0.0]
   !Z.Omargin = [0.0, 0.0]
   !Z.Window = [0.0, 0.0]
   !Z.Region = [0.0, 0.0]
   !Z.Charsize = 0.0
   !Z.Minor = 0L
   !Z.Tickv = Replicate(0.0, 30)
   !Z.Tickname = Replicate('', 30)
   !Z.Gridstyle = 0L
   !Z.Tickformat = ''
   !Z.Title = ''

end
;******************************************************************************


;******************************************************************************
; Procedure to set the differential shading.
pro Viz3D_DiffColor, sViz3DColors
   ; IDL colors common block (used by "LOADCT" and "STRETCH").
   common colors, r, g, b, cur_red, cur_green, cur_blue

   fdiffShade = Float(sViz3DColors.diffShade) / 100.0

   cR1 = (sViz3DColors.rSlice * (1.0 - fdiffshade))
   cG1 = (sViz3DColors.gSlice * (1.0 - fdiffshade))
   cB1 = (sViz3DColors.bSlice * (1.0 - fdiffshade))
   TVLCT, cR1, cG1, cB1, 0

   cR1 = (sViz3DColors.rSlice * (1.0 - fdiffshade)) + $
         (127.0 * fdiffshade)
   cG1 = (sViz3DColors.gSlice * (1.0 - fdiffshade)) + $
         (127.0 * fdiffshade)
   cB1 = (sViz3DColors.bSlice * (1.0 - fdiffshade)) + $
         (127.0 * fdiffshade)
   TVLCT, cR1, cG1, cB1, sViz3DColors.nColor

   cR1 = (sViz3DColors.rSlice * (1.0 - fdiffshade)) + $
         (255.0 * fdiffshade)
   cG1 = (sViz3DColors.gSlice * (1.0 - fdiffshade)) + $
         (255.0 * fdiffshade)
   cB1 = (sViz3DColors.bSlice * (1.0 - fdiffshade)) + $
         (255.0 * fdiffshade)
   TVLCT, cR1, cG1, cB1, (2 * sViz3DColors.nColor)

   TVLCT, cur_red, cur_green, cur_blue, /GET
   sViz3DColors.cR = cur_red
   sViz3DColors.cG = cur_green
   sViz3DColors.cB = cur_blue

   if (sViz3DColors.displayBits eq 24) then LOADCT, 0

end
;******************************************************************************


;******************************************************************************
; Function to load a color table and set the annotation colors.
; Returns the color state.
function Viz3D_LoadColor, CTAB=cTab, DIFFSHADE=diffShade
   ; IDL colors common block (used by "LOADCT" and "STRETCH").
   common colors, r, g, b, cur_red, cur_green, cur_blue

   ;Annotation color numbers :
   ; 0 - Black.
   ; 1 - Blue.
   ; 2 - Green.
   ; 3 - Red.
   ; 4 - Purple.
   ; 5 - Cyan.
   ; 6 - Yellow.
   ; 7 - White.

   nColor = ((!D.N_Colors < 256) - 9) / 5
   white  = LONG([255, 255, 255])
   yellow = LONG([255, 255,   0])
   cyan   = LONG([  0, 255, 255])
   purple = LONG([191,   0, 191])
   red    = LONG([255,   0,   0])
   green  = LONG([  0, 255,   0])
   blue   = LONG([ 63,  63, 255])
   black  = LONG([  0,   0,   0])

   if (N_ELEMENTS(cTab) le 3L) then cTab = [3,1,8]

   LOADCT, cTab(0)
   STRETCH, 0, (nColor-1)
   rSlice = cur_red(0:nColor-1)
   gSlice = cur_green(0:nColor-1)
   bSlice = cur_blue(0:nColor-1)

   LOADCT, cTab(1)
   STRETCH, 0, (nColor-1)
   rSurfc = cur_red(0:nColor-1)
   gSurfc = cur_green(0:nColor-1)
   bSurfc = cur_blue(0:nColor-1)

   LOADCT, cTab(2)
   STRETCH, 0, (nColor-1)
   rPrjct = cur_red(0:nColor-1)
   gPrjct = cur_green(0:nColor-1)
   bPrjct = cur_blue(0:nColor-1)

   if (!D.N_Colors gt 256) then displayBits = 24 else displayBits = 8

   cWhite24 = (white(2) * 256L^2L) + (white(1) * 256L) + white(0)
   cYellow24 = (yellow(2) * 256L^2L) + (yellow(1) * 256L) + yellow(0)
   cCyan24 = (cyan(2) * 256L^2L) + (cyan(1) * 256L) + cyan(0)
   cPurple24 = (purple(2) * 256L^2L) + (purple(1) * 256L) + purple(0)
   cRed24 = (red(2) * 256L^2L) + (red(1) * 256L) + red(0)
   cGreen24 = (green(2) * 256L^2L) + (green(1) * 256L) + green(0)
   cBlue24 = (blue(2) * 256L^2L) + (blue(1) * 256L) + blue(0)
   cBlack24 = (black(2) * 256L^2L) + (black(1) * 256L) + black(0)

   nCol = !D.N_Colors < 256

   cWhite = nCol - 1
   cYellow = nCol - 2
   cCyan = nCol - 3
   cPurple = nCol - 4
   cRed = nCol - 5
   cGreen = nCol - 6
   cBlue = nCol - 7
   cBlack = nCol - 8

   colors8 = [cBlack,cBlue,cGreen,cRed,cPurple,cCyan,cYellow,cWhite]
   colors24 = [cBlack24,cBlue24,cGreen24,cRed24,$
               cPurple24,cCyan24,cYellow24,cWhite24]

   TVLCT, rSurfc, gSurfc, bSurfc, (3 * nColor)
   TVLCT, rPrjct, gPrjct, bPrjct, (4 * nColor)
   TVLCT, REFORM(white, 1, 3), cWhite
   TVLCT, REFORM(yellow, 1, 3), cYellow
   TVLCT, REFORM(cyan, 1, 3), cCyan
   TVLCT, REFORM(purple, 1, 3), cPurple
   TVLCT, REFORM(red, 1, 3), cRed
   TVLCT, REFORM(green, 1, 3), cGreen
   TVLCT, REFORM(blue, 1, 3), cBlue
   TVLCT, REFORM(black, 1, 3), cBlack

   sViz3DColors = {SViz3DColors, nColor:nColor, diffShade:diffShade, $
                   colors8:colors8, colors24:colors24, $
                   displayBits:displayBits, $
                   cR:cur_red, cG:cur_green, cB:cur_blue, $
                   rSlice:rSlice, gSlice:gSlice, bSlice:bSlice}

   Viz3D_DiffColor, sViz3DColors

   TVLCT, cur_red, cur_green, cur_blue, /Get

   if (displayBits eq 24) then LOADCT, 0

   return, sViz3DColors

end
;******************************************************************************


;******************************************************************************
; Function to translate an annotation color
; number to an actual color table index.
function Viz3D_TransColor, sMainState, cIndex

   if ((sMainState.sColorState.displayBits eq 24) and (!D.Name NE 'Z')) then $
      return, sMainState.sColorState.colors24(cIndex) $
   else $
      return, sMainState.sColorState.colors8(cIndex)

end
;******************************************************************************


;******************************************************************************
; Function to set up the 3D view.
; Returns the view state.
function Viz3D_View, viewWin, $
                     XMAX=xMax, YMAX=yMax, ZMAX=zMax, $
                     ANG1=ang1, ANG2=ang2, ANG3=ang3, $
                     DIR1=dir1, DIR2=dir2, DIR3=dir3, $
                     ZOOM=zoomFac, ZSCALE=zScale, PERSP=pDist

   if (N_ELEMENTS(xMax) le 0L) then xMax = 1.0
   if (N_ELEMENTS(yMax) le 0L) then yMax = 1.0
   if (N_ELEMENTS(zMax) le 0L) then zMax = 1.0

   if (N_ELEMENTS(ang1) le 0L) then ang1 = 0.0
   if (N_ELEMENTS(ang2) le 0L) then ang2 = 0.0
   if (N_ELEMENTS(ang3) le 0L) then ang3 = 0.0
   if (N_ELEMENTS(dir1) le 0L) then dir1 = 'Z'
   if (N_ELEMENTS(dir2) le 0L) then dir2 = 'Y'
   if (N_ELEMENTS(dir3) le 0L) then dir3 = 'X'
   if (N_ELEMENTS(zoomFac) le 0L) then zoomFac = 1.0
   if (N_ELEMENTS(zScale) le 0L) then zScale = 1.0
   if (N_ELEMENTS(pDist) le 0L) then pDist = 0.0

   s_ang1 = Sin(ang1*!Dtor)
   c_ang1 = Cos(ang1*!Dtor)
   s_ang2 = Sin(ang2*!Dtor)
   c_ang2 = Cos(ang2*!Dtor)
   s_ang3 = Sin(ang3*!Dtor)
   c_ang3 = Cos(ang3*!Dtor)

   ident4 = FLTARR(4, 4)
   ident4([0,5,10,15]) = 1.0

   vTrans = ident4

   ; First translation.
   m4X4 = ident4
   m4X4([3,7,11]) = (-0.5)
   vTrans = TEMPORARY(vTrans) # m4X4

   ; Scale for aspect ratio.
   xRange = FLOAT(xMax)
   yRange = FLOAT(yMax)
   zRange = FLOAT(zMax) * (zScale > 1.0)
   maxRange = xRange > yRange > zRange
   xyzFac = [xRange, yRange, zRange] / maxRange
   m4X4 = ident4
   m4X4(0,0) = xyzFac(0)
   m4X4(1,1) = xyzFac(1)
   m4X4(2,2) = xyzFac(2)
   vTrans = TEMPORARY(vTrans) # m4X4

   ; Z scale.
   if (zScale lt 1.0) then begin
      m4X4 = ident4
      m4X4(2,2) = zScale
      vTrans = TEMPORARY(vTrans) # m4X4
   endif
   if (zScale gt 1.0) then begin
      zFac = 1.0 / zScale
      m4X4 = ident4
      m4X4(0,0) = zFac
      m4X4(1,1) = zFac
      vTrans = TEMPORARY(vTrans) # m4X4
   endif

   ; First rotation.
   m4x4 = ident4
   case dir1 of
      'X': begin ; Rotate about x.
         m4x4(1,1) = c_ang1
         m4x4(1,2) = s_ang1
         m4x4(2,1) = (-s_ang1)
         m4x4(2,2) = c_ang1
      end
      'Y': begin ; Rotate about y.
         m4x4(0,0) = c_ang1
         m4x4(0,2) = (-s_ang1)
         m4x4(2,0) = s_ang1
           m4x4(2,2) = c_ang1
      end
      'Z': begin ; Rotate about z.
         m4x4(0,0) = c_ang1
         m4x4(0,1) = s_ang1
         m4x4(1,0) = (-s_ang1)
         m4x4(1,1) = c_ang1
      end
   endcase
   vTrans = TEMPORARY(vTrans) # m4X4

   ; Second rotation.
   m4x4 = ident4
   case dir2 of
      'X': begin ; Rotate about x.
         m4x4(1,1) = c_ang2
         m4x4(1,2) = s_ang2
         m4x4(2,1) = (-s_ang2)
         m4x4(2,2) = c_ang2
      end
      'Y': begin ; Rotate about y.
         m4x4(0,0) = c_ang2
         m4x4(0,2) = (-s_ang2)
         m4x4(2,0) = s_ang2
         m4x4(2,2) = c_ang2
      end
      'Z': begin ; Rotate about z.
         m4x4(0,0) = c_ang2
         m4x4(0,1) = s_ang2
         m4x4(1,0) = (-s_ang2)
         m4x4(1,1) = c_ang2
      end
   endcase
   vTrans = TEMPORARY(vTrans) # m4X4

   ; Third rotation.
   m4x4 = ident4
   case dir3 of
      'X': begin ; Rotate about x.
         m4x4(1,1) = c_ang3
         m4x4(1,2) = s_ang3
         m4x4(2,1) = (-s_ang3)
         m4x4(2,2) = c_ang3
      end
      'Y': begin ; Rotate about y.
         m4x4(0,0) = c_ang3
         m4x4(0,2) = (-s_ang3)
         m4x4(2,0) = s_ang3
         m4x4(2,2) = c_ang3
      end
      'Z': begin ; Rotate about z.
         m4x4(0,0) = c_ang3
         m4x4(0,1) = s_ang3
         m4x4(1,0) = (-s_ang3)
         m4x4(1,1) = c_ang3
      end
   endcase
   vTrans = TEMPORARY(vTrans) # m4X4

   ; Perspective.
   if (pDist gt 0.0) then begin
      pD2 = SQRT(pDist)
      m4X4 = ident4
      m4X4(2,3) = (-1.0 / pD2)
      vTrans = TEMPORARY(vTrans) # m4X4
   endif

   ; Zoom.
   m4X4 = ident4
   m4X4(0,0) = zoomFac
   m4X4(1,1) = zoomFac
   m4X4(2,2) = zoomFac
   vTrans = TEMPORARY(vTrans) # m4X4

   ; Zoom down so that all 8 corners of the
   ; unit cube are within the view.
   corners = Fltarr(3,8)
   corners(*,0) = [0,0,0]
   corners(*,1) = [1,0,0]
   corners(*,2) = [1,1,0]
   corners(*,3) = [0,1,0]
   corners(*,4) = [0,0,1]
   corners(*,5) = [1,0,1]
   corners(*,6) = [1,1,1]
   corners(*,7) = [0,1,1]
   corners = VERT_T3D(corners, Matrix=vTrans, /NO_COPY)
   sFac = (1.0 / (1.0 + (2.0 * (MAX(ABS(corners(0:1,*))) - 0.5)))) < 1.0
   if (sFac lt 1.0) then begin
      m4X4 = ident4
      m4X4(0,0) = sFac
      m4X4(1,1) = sFac
      m4X4(2,2) = sFac
      vTrans = TEMPORARY(vTrans) # m4X4
   endif

   m4X4 = ident4
   m4X4(2,2) = (0.1)
   vTrans = TEMPORARY(vTrans) # m4X4

   ; Last translation.
   m4X4 = ident4
   m4X4([3,7,11]) = (0.5)
   vTrans = TEMPORARY(vTrans) # m4X4

   ; System variables.
   !P.T = vTrans
   !X.S = [0.0, (1.0 / xRange)]
   !Y.S = [0.0, (1.0 / yRange)]
   !Z.S = [0.0, (1.0 / FLOAT(zMax))]

   sViz3DView = {SViz3DView, viewWin:viewWin, $
                 xMax:xMax, yMax:yMax, zMax:zMax, $
                 ang1:ang1, ang2:ang2, ang3:ang3, $
                 dir1:dir1, dir2:dir2, dir3:dir3, $
                 zoomFac:zoomFac, zScale:zScale, pDist:pDist, $
                 vTrans:vTrans, invTrans:INVERT(vTrans)}

   return, sViz3DView

end
;******************************************************************************


;******************************************************************************
; Procedure to save the depth buffer information for the cube's outer surface.
; This depth info is used to find the XYZ data coordinate on the surface of
; the cube, given a user selected screen XY coordinate.
pro Viz3D_FillDepth, sMainState

   SET_PLOT, 'Z'
   image = TVRD()
   depth = TVRD(CHANNEL=1, /WORDS)

   ; Front faces.
   ERASE, 0
   norm = [[0.0,0.0,0.0], [0.0,0.0,(-1.0)]]
   norm = VERT_T3D(norm, /NO_COPY)
   if ((norm(2,1)-norm(2,0)) gt 0.0) then $
      POLYFILL, [0,1,1,0], [0,0,1,1], [0,0,0,0], /NORMAL, /T3D, COLOR=1

   norm = [[0.0,0.0,0.0], [0.0,(-1.0),0.0]]
   norm = VERT_T3D(norm, /NO_COPY)
   if ((norm(2,1)-norm(2,0)) gt 0.0) then $
      POLYFILL, [0,1,1,0], [0,0,0,0], [0,0,1,1], /NORMAL, /T3D, COLOR=2

   norm = [[0.0,0.0,0.0], [(-1.0),0.0,0.0]]
   norm = VERT_T3D(norm, /NO_COPY)
   if ((norm(2,1)-norm(2,0)) gt 0.0) then $
      POLYFILL, [0,0,0,0], [0,1,1,0], [0,0,1,1], /NORMAL, /T3D, COLOR=3

   norm = [[1.0,1.0,1.0], [1.0,1.0,2.0]]
   norm = VERT_T3D(norm, /NO_COPY)
   if ((norm(2,1)-norm(2,0)) gt 0.0) then $
      POLYFILL, [0,1,1,0], [0,0,1,1], [1,1,1,1], /NORMAL, /T3D, COLOR=1

   norm = [[1.0,1.0,1.0], [1.0,2.0,1.0]]
   norm = VERT_T3D(norm, /NO_COPY)
   if ((norm(2,1)-norm(2,0)) gt 0.0) then $
      POLYFILL, [0,1,1,0], [1,1,1,1], [0,0,1,1], /NORMAL, /T3D, COLOR=2

   norm = [[1.0,1.0,1.0], [2.0,1.0,1.0]]
   norm = VERT_T3D(norm, /NO_COPY)
   if ((norm(2,1)-norm(2,0)) gt 0.0) then $
      POLYFILL, [1,1,1,1], [0,1,1,0], [0,0,1,1], /NORMAL, /T3D, COLOR=3

   *(sMainState.hFrontImage) = TVRD()
   *(sMainState.hFrontDepth) = TVRD(CHANNEL=1, /WORDS)

   ; Back faces.
   ERASE, 0

   norm = [[0.0,0.0,0.0], [0.0,0.0,(-1.0)]]
   norm = VERT_T3D(norm, /NO_COPY)
   if ((norm(2,1)-norm(2,0)) lt 0.0) then $
      POLYFILL, [0,1,1,0], [0,0,1,1], [0,0,0,0], /NORMAL, /T3D, COLOR=1

   norm = [[0.0,0.0,0.0], [0.0,(-1.0),0.0]]
   norm = VERT_T3D(norm, /NO_COPY)
   if ((norm(2,1)-norm(2,0)) lt 0.0) then $
      POLYFILL, [0,1,1,0], [0,0,0,0], [0,0,1,1], /NORMAL, /T3D, COLOR=2

   norm = [[0.0,0.0,0.0], [(-1.0),0.0,0.0]]
   norm = VERT_T3D(norm, /NO_COPY)
   if ((norm(2,1)-norm(2,0)) lt 0.0) then $
      POLYFILL, [0,0,0,0], [0,1,1,0], [0,0,1,1], /NORMAL, /T3D, COLOR=3

   norm = [[1.0,1.0,1.0], [1.0,1.0,2.0]]
   norm = VERT_T3D(norm, /NO_COPY)
   if ((norm(2,1)-norm(2,0)) lt 0.0) then $
      POLYFILL, [0,1,1,0], [0,0,1,1], [1,1,1,1], /NORMAL, /T3D, COLOR=1

   norm = [[1.0,1.0,1.0], [1.0,2.0,1.0]]
   norm = VERT_T3D(norm, /NO_COPY)
   if ((norm(2,1)-norm(2,0)) lt 0.0) then $
      POLYFILL, [0,1,1,0], [1,1,1,1], [0,0,1,1], /NORMAL, /T3D, COLOR=2

   norm = [[1.0,1.0,1.0], [2.0,1.0,1.0]]
   norm = VERT_T3D(norm, /NO_COPY)
   if ((norm(2,1)-norm(2,0)) lt 0.0) then $
      POLYFILL, [1,1,1,1], [0,1,1,0], [0,0,1,1], /NORMAL, /T3D, COLOR=3

   *(sMainState.hBackImage) = TVRD()
   *(sMainState.hBackDepth) = TVRD(CHANNEL=1, /WORDS)

   TV, TEMPORARY(image)
   TV, TEMPORARY(depth), CHANNEL=1, /WORDS
   SET_PLOT, sMainState.screenDevice

end
;******************************************************************************


;******************************************************************************
; Function to check out the current data.
function Viz3D_GetData, sMainState, curData

   if (N_ELEMENTS(curData) LE 0L) then curData = sMainState.curData

   hData3D = (*(sMainState.hDataList))(curData)

   return, TEMPORARY(*hData3D)

end
;******************************************************************************


;******************************************************************************
; Function to check in the current data.
pro Viz3D_PutData, data3D, sMainState, curData

   if (N_ELEMENTS(curData) LE 0L) then curData = sMainState.curData

   hData3D = (*(sMainState.hDataList))(curData)

   *hData3D = TEMPORARY(data3D)

end
;******************************************************************************


;******************************************************************************
; Function to scale the data into the specified thereshold range.
function Viz3D_ScaleData, data3D, sMainState

   lTh = (*(sMainState.sThreshState.hLowThresh))(sMainState.curData)
   hTh = (*(sMainState.sThreshState.hHighThresh))(sMainState.curData)
   return, BYTSCL(TEMPORARY(data3D), MIN=lTh, MAX=hTh, TOP=255B)

end
;******************************************************************************


;******************************************************************************
; Procedure to compute the intersection outline of the cube and the
; oblique slicing plane.   This procedure also draws the outline in
; the current window (either the big view window, or the small slice
; status window).
; The VERT_PLANE keyword returns the intersection coordinates in the
; image (oblique) plane.   The VERT_3D keyword returns the intersection
; coordinates as absolute 3D coordinates.   The plane and 3D coordinates
; will be necessary when drawing the data on the oblique plane.

pro Viz3D_DrawSliceOblique, obliqCenter, obliqNormal, fillColor, $
       edgeColor, rangeX, rangeY, rangeZ, zScale, NORMCOLOR=normColor, $
       SKIP_FILL=skipFill, VERT_PLANE=vertPlane, VERT_3D=vert3D

   ; Scale the oblique normal vector.
   zDim = rangeZ * zScale
   maxDim = rangeX > rangeY > zDim
   sliceNormU = obliqNormal
   sliceNormU(0) = obliqNormal(0) * rangeX / maxDim
   sliceNormU(1) = obliqNormal(1) * rangeY / maxDim
   sliceNormU(2) = obliqNormal(2) * zDim / maxDim

   ; Determine the rotation angles and build the transformation matrices.

   if ((sliceNormU(0) eq 0.0) and (sliceNormU(1) eq 0.0)) then angZ = 0.0 $
   else angZ = ATAN(sliceNormU(1), sliceNormU(0))

   pDistance = SQRT(sliceNormU(0)^2 + sliceNormU(1)^2)
   if ((pDistance eq 0.0) and (sliceNormU(2) eq 0.0)) then angY = 0.0 $
   else angY = ATAN(sliceNormU(2), pDistance)

   sZ = SIN(-angZ)
   cZ = COS(-angZ)
   sY = SIN(angY-(!PI/2.0))
   cY = COS(angY-(!PI/2.0))

   ident4 = FLTARR(4, 4)
   ident4([0,5,10,15]) = 1.0

   vTrans1 = ident4
   ; First translation.
   m4X4 = ident4
   m4X4([3,7,11]) = (-obliqCenter)
   vTrans1 = TEMPORARY(vTrans1) # m4X4
   ; Z rotation.
   m4X4 = ident4
   m4x4(0,0) = cZ
   m4x4(0,1) = sZ
   m4x4(1,0) = (-sZ)
   m4x4(1,1) = cZ
   vTrans1 = TEMPORARY(vTrans1) # m4X4
   ; Y rotation.
   m4X4 = ident4
   m4x4(0,0) = cY
   m4x4(0,2) = (-sY)
   m4x4(2,0) = sY
   m4x4(2,2) = cY
   vTrans1 = TEMPORARY(vTrans1) # m4X4

   sZ = SIN(angZ)
   cZ = COS(angZ)
   sY = SIN((!PI/2.0)-angY)
   cY = COS((!PI/2.0)-angY)

   ; Create vTrans2 as the inverse of vTrans1.
   vTrans2 = ident4
   ; Y rotation.
   m4X4 = ident4
   m4x4(0,0) = cY
   m4x4(0,2) = (-sY)
   m4x4(2,0) = sY
   m4x4(2,2) = cY
   vTrans2 = TEMPORARY(vTrans2) # m4X4
   ; Z rotation.
   m4X4 = ident4
   m4x4(0,0) = cZ
   m4x4(0,1) = sZ
   m4x4(1,0) = (-sZ)
   m4x4(1,1) = cZ
   vTrans2 = TEMPORARY(vTrans2) # m4X4
   ; First translation.
   m4X4 = ident4
   m4X4([3,7,11]) = (+obliqCenter)
   vTrans2 = TEMPORARY(vTrans2) # m4X4

   verts = FLTARR(3, 8, /NOZERO)
   verts(0, 0) = [0,0,0]
   verts(0, 1) = [1,0,0]
   verts(0, 2) = [1,1,0]
   verts(0, 3) = [0,1,0]
   verts(0, 4) = [0,0,1]
   verts(0, 5) = [1,0,1]
   verts(0, 6) = [1,1,1]
   verts(0, 7) = [0,1,1]
   verts = VERT_T3D(verts, /NO_COPY, MATRIX=vTrans1)

   ; Search for intersection points.

   newXYZ = [0.0, 0.0, 0.0]

   ; Verts 0, 1.
   if ((verts(2,0)*verts(2,1)) le 0.0) then begin
      lenFac = verts(2,0) / (verts(2,0) - verts(2,1))
      newXYZ = [Temporary(newXYZ), verts(*,0) + $
         ((verts(*,1)-verts(*,0)) * lenFac)]
   endif
   ; Verts 0, 3.
   if ((verts(2,0)*verts(2,3)) le 0.0) then begin
      lenFac = verts(2,0) / (verts(2,0) - verts(2,3))
      newXYZ = [Temporary(newXYZ), verts(*,0) + $
         ((verts(*,3)-verts(*,0)) * lenFac)]
   endif
   ; Verts 0, 4.
   if ((verts(2,0)*verts(2,4)) le 0.0) then begin
      lenFac = verts(2,0) / (verts(2,0) - verts(2,4))
      newXYZ = [Temporary(newXYZ), verts(*,0) + $
         ((verts(*,4)-verts(*,0)) * lenFac)]
   endif
   ; Verts 1, 2.
   if ((verts(2,1)*verts(2,2)) le 0.0) then begin
      lenFac = verts(2,1) / (verts(2,1) - verts(2,2))
      newXYZ = [Temporary(newXYZ), verts(*,1) + $
         ((verts(*,2)-verts(*,1)) * lenFac)]
   endif
   ; Verts 1, 5.
   if ((verts(2,1)*verts(2,5)) le 0.0) then begin
      lenFac = verts(2,1) / (verts(2,1) - verts(2,5))
      newXYZ = [Temporary(newXYZ), verts(*,1) + $
         ((verts(*,5)-verts(*,1)) * lenFac)]
   endif
   ; Verts 2, 3.
   if ((verts(2,2)*verts(2,3)) le 0.0) then begin
      lenFac = verts(2,2) / (verts(2,2) - verts(2,3))
      newXYZ = [Temporary(newXYZ), verts(*,2) + $
         ((verts(*,3)-verts(*,2)) * lenFac)]
   endif
   ; Verts 2, 6.
   if ((verts(2,2)*verts(2,6)) le 0.0) then begin
      lenFac = verts(2,2) / (verts(2,2) - verts(2,6))
      newXYZ = [Temporary(newXYZ), verts(*,2) + $
         ((verts(*,6)-verts(*,2)) * lenFac)]
   endif
   ; Verts 3, 7.
   if ((verts(2,3)*verts(2,7)) le 0.0) then begin
      lenFac = verts(2,3) / (verts(2,3) - verts(2,7))
      newXYZ = [Temporary(newXYZ), verts(*,3) + $
         ((verts(*,7)-verts(*,3)) * lenFac)]
   endif
   ; Verts 4, 5.
   if ((verts(2,4)*verts(2,5)) le 0.0) then begin
      lenFac = verts(2,4) / (verts(2,4) - verts(2,5))
      newXYZ = [Temporary(newXYZ), verts(*,4) + $
         ((verts(*,5)-verts(*,4)) * lenFac)]
   endif
   ; Verts 4, 7.
   if ((verts(2,4)*verts(2,7)) le 0.0) then begin
      lenFac = verts(2,4) / (verts(2,4) - verts(2,7))
      newXYZ = [Temporary(newXYZ), verts(*,4) + $
         ((verts(*,7)-verts(*,4)) * lenFac)]
   endif
   ; Verts 5, 6.
   if ((verts(2,5)*verts(2,6)) le 0.0) then begin
      lenFac = verts(2,5) / (verts(2,5) - verts(2,6))
      newXYZ = [Temporary(newXYZ), verts(*,5) + $
         ((verts(*,6)-verts(*,5)) * lenFac)]
   endif
   ; Verts 6, 7.
   if ((verts(2,6)*verts(2,7)) le 0.0) then begin
      lenFac = verts(2,6) / (verts(2,6) - verts(2,7))
      newXYZ = [Temporary(newXYZ), verts(*,6) + $
         ((verts(*,7)-verts(*,6)) * lenFac)]
   endif

   if (N_ELEMENTS(newXYZ) le 3L) then return

   newXYZ = newXYZ(3:*)
   nPoints = N_ELEMENTS(newXYZ) / 3
   newXYZ = REFORM(TEMPORARY(newXYZ), 3, nPoints)

   ; newXYZ is now a 3,n array of 3D intersection vertices.

   ; Sort the vertices into counter-clockwise order.

   newXMin = MIN(newXYZ(0,*), MAX=newXMax)
   newYMin = MIN(newXYZ(1,*), MAX=newYMax)
   cX = (newXMin + newXMax) / 2.0
   cY = (newYMin + newYMax) / 2.0

   angDiff = ATAN(newXYZ(1,*)-cY, newXYZ(0,*)-cX)

   negInd = WHERE(angDiff lt 0.0)
   if (negInd(0) ge 0L) then angDiff(negInd) = angDiff(negInd) + (2.0*!PI)
   sortInd = SORT(Temporary(angDiff))
   newXYZ = newXYZ(*, Temporary(sortInd))

   ; Save the intersection vertices on the image (oblique) plane.
   vertPlane = newXYZ
   ; Save the intersection vertices as 3D data coordinates.
   newXYZ = VERT_T3D(newXYZ, MATRIX=vTrans2, /NO_COPY)
   vert3D = newXYZ
   ; vertPlane and vert3D are returned to the calling procedure via
   ; the VERT_PLANE and VERT_3D keywords.

   ; Draw the intersection outline (optionally filled).
   if not(KEYWORD_SET(skipFill)) then $
      POLYFILL, newXYZ, /NORMAL, /T3D, COLOR=fillColor
   PLOTS, newXYZ, /NORMAL, /T3D, COLOR=edgeColor
   PLOTS, newXYZ(*,0), /NORMAL, /T3D, COLOR=edgeColor, /CONTINUE

   if (N_ELEMENTS(normColor) gt 0L) then begin
      ; Draw the surface normal vector.
      sliceNormU = obliqNormal
      sliceNormU(0) = obliqNormal(0) * maxDim / rangeX
      sliceNormU(1) = obliqNormal(1) * maxDim / rangeY
      sliceNormU(2) = obliqNormal(2) * maxDim / zDim

      normPoint = obliqCenter + sliceNormU

      if (normPoint(0) lt 0.0) then begin
         lenFac = obliqCenter(0) / (obliqCenter(0) - normPoint(0))
         normPoint = obliqCenter + (sliceNormU * lenFac)
      endif
      if (normPoint(0) gt 1.0) then begin
         lenFac = (1.0 - obliqCenter(0)) / (normPoint(0) - obliqCenter(0))
         normPoint = obliqCenter + (sliceNormU * lenFac)
      endif

      if (normPoint(1) lt 0.0) then begin
         lenFac = obliqCenter(1) / (obliqCenter(1) - normPoint(1))
         normPoint = obliqCenter + (sliceNormU * lenFac)
      endif
      if (normPoint(1) gt 1.0) then begin
         lenFac = (1.0 - obliqCenter(1)) / (normPoint(1) - obliqCenter(1))
         normPoint = obliqCenter + (sliceNormU * lenFac)
      endif

      if (normPoint(2) lt 0.0) then begin
         lenFac = obliqCenter(2) / (obliqCenter(2) - normPoint(2))
         normPoint = obliqCenter + (sliceNormU * lenFac)
      endif
      if (normPoint(2) gt 1.0) then begin
         lenFac = (1.0 - obliqCenter(2)) / (normPoint(2) - obliqCenter(2))
         normPoint = obliqCenter + (sliceNormU * lenFac)
      endif

      PLOTS, [obliqCenter(0), normPoint(0)], $
             [obliqCenter(1), normPoint(1)], $
             [obliqCenter(2), normPoint(2)], $
         /NORMAL, /T3D, COLOR=normColor
   endif

end
;******************************************************************************


;******************************************************************************
; Procedure to draw the cube and current slicing plane outline in the
; small slice status window.
pro Viz3D_SliceShow, sMainState

   fillColor = Viz3D_TransColor(sMainState, sMainState.sSliceState.fillColor)
   edgeColor = Viz3D_TransColor(sMainState, sMainState.sSliceState.edgeColor)
   normColor = Viz3D_TransColor(sMainState, sMainState.sSliceState.normColor)
   backColor = Viz3D_TransColor(sMainState, sMainState.backColor)

   WSET, sMainState.sSliceState.sliceWin
   ERASE

   case (sMainState.sSliceState.planeMode*sMainState.sSliceState.orthoDir) of
      0: begin ; Oblique
         Viz3D_DrawSliceOblique, sMainState.sSliceState.obliqCenter, $
            sMainState.sSliceState.obliqNormal, fillColor, edgeColor, $
            sMainState.sViewState.xMax, sMainState.sViewState.yMax, $
            sMainState.sViewState.zMax, sMainState.sViewState.zScale, $
            NORMCOLOR=normColor
      end
      1: begin ; X
         POLYFILL, sMainState.sSliceState.orthoPos, [0,1,1,0], [0,0,1,1], $
            /NORMAL, /T3D, COLOR=fillColor
         PLOTS, sMainState.sSliceState.orthoPos, [0,1,1,0,0], [0,0,1,1,0], $
            /NORMAL, /T3D, COLOR=edgeColor
      end
      2: begin ; Y
         POLYFILL, [0,1,1,0], sMainState.sSliceState.orthoPos, [0,0,1,1], $
            /NORMAL, /T3D, COLOR=fillColor
         PLOTS, [0,1,1,0,0], sMainState.sSliceState.orthoPos, [0,0,1,1,0], $
            /NORMAL, /T3D, COLOR=edgeColor
      end
      3: begin ; Z
         POLYFILL, [0,1,1,0], [0,0,1,1], sMainState.sSliceState.orthoPos, $
            /NORMAL, /T3D, COLOR=fillColor
         PLOTS, [0,1,1,0,0], [0,0,1,1,0], sMainState.sSliceState.orthoPos, $
            /NORMAL, /T3D, COLOR=edgeColor
      end
   endcase

   ; Draw the cube outline in the current window.
   Viz3D_DrawCube, sMainState, /SKIPBACK, /SKIPAXIS, /DIRECT

   WSET, sMainState.mainWin
end
;******************************************************************************


;******************************************************************************
; Procedure to draw the cube outline in the current window,
pro Viz3D_DrawCube, sMainState, SKIPCUBE=skipCube, SKIPAXIS=skipAxis, $
    SKIPFRONT=skipFront, SKIPBACK=skipBack, DIRECT=directDraw, AXIS=drawAxis

   axisOn = sMainState.axisOn
   if (KEYWORD_SET(drawAxis)) then axisOn = 1

   cubeColor = Viz3D_TransColor(sMainState, sMainState.cubeColor)
   axisColor = Viz3D_TransColor(sMainState, sMainState.axisColor)

   if ((axisOn) and not(KEYWORD_SET(skipAxis))) then $
      lineColor = axisColor else lineColor = cubeColor

   verts = FLTARR(3, 8, /NOZERO)
   verts(0, 0) = [0,0,0]
   verts(0, 1) = [1,0,0]
   verts(0, 2) = [1,1,0]
   verts(0, 3) = [0,1,0]
   verts(0, 4) = [0,0,1]
   verts(0, 5) = [1,0,1]
   verts(0, 6) = [1,1,1]
   verts(0, 7) = [0,1,1]
   verts = VERT_T3D(verts, /NO_COPY)

   cubeDrawn = 0
   if (((sMainState.cubeOn) and not(KEYWORD_SET(skipCube))) or $
      KEYWORD_SET(directDraw)) then begin
      ; Draw the cube.

      cubeDrawn = 1
      lCol = [lineColor, lineColor, cubeColor, cubeColor, lineColor]

      if not(KEYWORD_SET(skipBack)) then begin
         ; Draw the "hidden" lines as dotted.

         cp = CROSSP(verts(*,1)-verts(*,0), verts(*,4)-verts(*,0))
         hideFace0154 = (cp(2) lt 0.0)
         cp = CROSSP(verts(*,4)-verts(*,0), verts(*,3)-verts(*,0))
         hideFace0473 = (cp(2) lt 0.0)
         cp = CROSSP(verts(*,3)-verts(*,0), verts(*,1)-verts(*,0))
         hideFace0321 = (cp(2) lt 0.0)
         cp = CROSSP(verts(*,2)-verts(*,6), verts(*,7)-verts(*,6))
         hideFace6237 = (cp(2) lt 0.0)
         cp = CROSSP(verts(*,5)-verts(*,6), verts(*,2)-verts(*,6))
         hideFace6514 = (cp(2) lt 0.0)
         cp = CROSSP(verts(*,7)-verts(*,6), verts(*,5)-verts(*,6))
         hideFace6745 = (cp(2) lt 0.0)

         if (hideFace0154 and hideFace0321) then $
            PLOTS, verts(*, [0,1]), /NORMAL, T3D=0, COLOR=lineColor, $
               LINESTYLE=1
         if (hideFace0473 and hideFace0321) then $
            PLOTS, verts(*, [0,3]), /NORMAL, T3D=0, COLOR=lineColor, $
               LINESTYLE=1
         if (hideFace0154 and hideFace0473) then $
            PLOTS, verts(*, [0,4]), /NORMAL, T3D=0, COLOR=lineColor, $
               LINESTYLE=1

         if (hideFace0321 and hideFace6514) then $
            PLOTS, verts(*, [1,2]), /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace0321 and hideFace6237) then $
            PLOTS, verts(*, [2,3]), /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace0154 and hideFace6514) then $
            PLOTS, verts(*, [1,5]), /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace0473 and hideFace6237) then $
            PLOTS, verts(*, [3,7]), /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace6745 and hideFace0473) then $
            PLOTS, verts(*, [4,7]), /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace6745 and hideFace0154) then $
            PLOTS, verts(*, [4,5]), /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace6745 and hideFace6514) then $
            PLOTS, verts(*, [5,6]), /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace6745 and hideFace6237) then $
            PLOTS, verts(*, [6,7]), /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace6514 and hideFace6237) then $
            PLOTS, verts(*, [2,6]), /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
      endif

      if not(KEYWORD_SET(skipFront)) then begin
         ; Draw the remaining cube edges as solid (linestyle 0).

         cp = CROSSP(verts(*,1)-verts(*,0), verts(*,4)-verts(*,0))
         if (cp(2) gt 0.0) then $
            PLOTS, verts(*, [0,1,5,4,0]), /NORMAL, T3D=0, COLOR=lCol

         cp = CROSSP(verts(*,4)-verts(*,0), verts(*,3)-verts(*,0))
         if (cp(2) gt 0.0) then $
            PLOTS, verts(*, [0,4,7,3,0]), /NORMAL, T3D=0, COLOR=lCol

         cp = CROSSP(verts(*,3)-verts(*,0), verts(*,1)-verts(*,0))
         if (cp(2) gt 0.0) then $
            PLOTS, verts(*, [0,3,2,1,0]), /NORMAL, T3D=0, COLOR=lCol

         cp = CROSSP(verts(*,7)-verts(*,6), verts(*,5)-verts(*,6))
         if (cp(2) gt 0.0) then $
            PLOTS, verts(*, [6,7,4,5,6]), /NORMAL, T3D=0, COLOR=cubeColor

         cp = CROSSP(verts(*,5)-verts(*,6), verts(*,2)-verts(*,6))
         if (cp(2) gt 0.0) then $
            PLOTS, verts(*, [6,5,1,2,6]), /NORMAL, T3D=0, COLOR=cubeColor

         cp = CROSSP(verts(*,2)-verts(*,6), verts(*,7)-verts(*,6))
         if (cp(2) gt 0.0) then $
            PLOTS, verts(*, [6,2,3,7,6]), /NORMAL, T3D=0, COLOR=cubeColor
      endif
   endif

   if (((axisOn) and not(KEYWORD_SET(skipAxis))) and $
      (not(KEYWORD_SET(skipFront)) and not(KEYWORD_SET(directDraw)))) $
      then begin
      ; Draw the axis.

      if (cubeDrawn eq 0) then begin
         PLOTS, verts(*, [0,1]), /NORMAL, T3D=0, COLOR=axisColor
         PLOTS, verts(*, [0,3]), /NORMAL, T3D=0, COLOR=axisColor
         PLOTS, verts(*, [0,4]), /NORMAL, T3D=0, COLOR=axisColor
      endif

      cp = CROSSP(verts(*,5)-verts(*,1), verts(*,0)-verts(*,1))
      showPlane1 = (cp(2) gt 0.0)
      cp = CROSSP(verts(*,0)-verts(*,1), verts(*,2)-verts(*,1))
      showPlane2 = (cp(2) gt 0.0)
      cp = CROSSP(verts(*,2)-verts(*,1), verts(*,5)-verts(*,1))
      showPlane3 = (cp(2) gt 0.0)
      if ((showPlane1+showPlane2+showPlane3) ge 1) then begin
         PLOTS, [0.95,1.0,0.95], [0.01, 0.0, -0.01], [0,0,0], /NORMAL, $
            /T3D, COLOR=axisColor
         PLOTS, [0.95,1.0,0.95], [0,0,0], [0.01, 0.0, -0.01], /NORMAL, $
            /T3D, COLOR=axisColor
         XYOUTS, 1.05, 0.0, Z=0.0, 'X', /NORMAL, /T3D, $
            CHARSIZE=24.0/FLOAT(!D.Y_Ch_Size), COLOR=axisColor, $
            ALIGNMENT=0.5, TEXT_AXES=2
      endif

      cp = CROSSP(verts(*,0)-verts(*,3), verts(*,7)-verts(*,3))
      showPlane1 = (cp(2) gt 0.0)
      cp = CROSSP(verts(*,7)-verts(*,3), verts(*,2)-verts(*,3))
      showPlane2 = (cp(2) gt 0.0)
      cp = CROSSP(verts(*,2)-verts(*,3), verts(*,0)-verts(*,3))
      showPlane3 = (cp(2) gt 0.0)
      if ((showPlane1+showPlane2+showPlane3) ge 1) then begin
         PLOTS, [0.01, 0.0, -0.01], [0.95,1.0,0.95], [0,0,0], /NORMAL, $
            /T3D, COLOR=axisColor
         PLOTS, [0,0,0], [0.95,1.0,0.95], [0.01, 0.0, -0.01], /NORMAL, $
            /T3D, COLOR=axisColor
         XYOUTS, 0.0, 1.05, Z=0.0, 'Y', /NORMAL, /T3D, $
            CHARSIZE=24.0/FLOAT(!D.Y_Ch_Size), COLOR=axisColor, $
            ALIGNMENT=0.5, TEXT_AXES=1
      endif

      cp = CROSSP(verts(*,5)-verts(*,4), verts(*,7)-verts(*,4))
      showPlane1 = (cp(2) gt 0.0)
      cp = CROSSP(verts(*,7)-verts(*,4), verts(*,0)-verts(*,4))
      showPlane2 = (cp(2) gt 0.0)
      cp = CROSSP(verts(*,0)-verts(*,4), verts(*,5)-verts(*,4))
      showPlane3 = (cp(2) gt 0.0)
      if ((showPlane1+showPlane2+showPlane3) ge 1) then begin
         PLOTS, [0,0,0], [0.01, 0.0, -0.01], [0.95,1.0,0.95], /NORMAL, $
            /T3D, COLOR=axisColor
         PLOTS, [0.01, 0.0, -0.01], [0,0,0], [0.95,1.0,0.95], /NORMAL, $
            /T3D, COLOR=axisColor
         XYOUTS, 0.0, 0.0, Z=1.05, 'Z', /NORMAL, /T3D, $
            CHARSIZE=24.0/FLOAT(!D.Y_Ch_Size), COLOR=axisColor, $
            ALIGNMENT=0.5, TEXT_AXES=0
      endif
   endif

   EMPTY

end
;******************************************************************************


;******************************************************************************
pro Viz3D_DrawData, sMainState

   SET_PLOT, 'Z'
   Viz3D_DrawCube, sMainState, /SKIPFRONT
   img = TVRD()
   SET_PLOT, sMainState.screenDevice
   WSET, sMainState.pixWin
   if (sMainState.sColorState.displayBits eq 8) then $
      TV, TEMPORARY(img) $
   else begin
      img3 = BYTARR(3, sMainState.winX, sMainState.winY, /NOZERO)
      img3(0, *, *) = sMainState.sColorState.cR(img)
      img3(1, *, *) = sMainState.sColorState.cG(img)
      img3(2, *, *) = sMainState.sColorState.cB(Temporary(img))
      TV, Temporary(img3), TRUE=1
   endelse
   Viz3D_DrawCube, sMainState, /SKIPBACK
   EMPTY
   WSET, sMainState.mainWin
   DEVICE, COPY=[0, 0, sMainState.winX, sMainState.winY, 0, 0, $
                 sMainState.pixWin]
   sMainState.cleanupView = 0
   sMainState.cleanupBuffer = 0

end
;******************************************************************************


;******************************************************************************
; Procedure to draw the block outline.
pro Viz3D_BlockShow, sMainState, DIRECT=directDraw

   if not(KEYWORD_SET(directDraw)) then $
      WSET, sMainState.sBlockState.blockWin

   x1 = sMainState.sBlockState.c1(0)
   y1 = sMainState.sBlockState.c1(1)
   z1 = sMainState.sBlockState.c1(2)
   x2 = sMainState.sBlockState.c2(0)
   y2 = sMainState.sBlockState.c2(1)
   z2 = sMainState.sBlockState.c2(2)
   blockColor = $
      Viz3D_TransColor(sMainState, sMainState.sBlockState.blockColor)

   ; Draw the cube outline in the current window.
   if not(KEYWORD_SET(directDraw)) then begin
      ERASE
      Viz3D_DrawCube, sMainState, /DIRECT, /AXIS
   endif else begin

      c1Color = Viz3D_TransColor(sMainState, sMainState.sBlockState.c1Color)
      c2Color = Viz3D_TransColor(sMainState, sMainState.sBlockState.c2Color)
      xM = sMainState.sViewState.xMax
      yM = sMainState.sViewState.yMax
      zM = sMainState.sViewState.zMax

      norm = [[x1,y1,z1], [x1,y1,z1-1]]
      norm = VERT_T3D(norm, /NO_COPY)
      if ((norm(2,1)-norm(2,0)) gt 0.0) then begin
         PLOTS, [x1,x1], [y1,y1], [z1, 0], /DATA, /T3D, COLOR=c1Color
         PLOTS, [x1,x1], [y1,y1], [z1,zM], /DATA, /T3D, COLOR=c1Color, $
            LINESTYLE=1
      endif else begin
         PLOTS, [x1,x1], [y1,y1], [z1, 0], /DATA, /T3D, COLOR=c1Color, $
            LINESTYLE=1
         PLOTS, [x1,x1], [y1,y1], [z1,zM], /DATA, /T3D, COLOR=c1Color
      endelse
      norm = [[x1,y1,z1], [x1,y1-1,z1]]
      norm = VERT_T3D(norm, /NO_COPY)
      if ((norm(2,1)-norm(2,0)) gt 0.0) then begin
         PLOTS, [x1,x1], [y1, 0], [z1,z1], /DATA, /T3D, COLOR=c1Color
         PLOTS, [x1,x1], [y1,yM], [z1,z1], /DATA, /T3D, COLOR=c1Color, $
            LINESTYLE=1
      endif else begin
         PLOTS, [x1,x1], [y1, 0], [z1,z1], /DATA, /T3D, COLOR=c1Color, $
            LINESTYLE=1
         PLOTS, [x1,x1], [y1,yM], [z1,z1], /DATA, /T3D, COLOR=c1Color
      endelse
      norm = [[x1,y1,z1], [x1-1,y1,z1]]
      norm = VERT_T3D(norm, /NO_COPY)
      if ((norm(2,1)-norm(2,0)) gt 0.0) then begin
         PLOTS, [x1, 0], [y1,y1], [z1,z1], /DATA, /T3D, COLOR=c1Color
         PLOTS, [x1,xM], [y1,y1], [z1,z1], /DATA, /T3D, COLOR=c1Color, $
            LINESTYLE=1
      endif else begin
         PLOTS, [x1, 0], [y1,y1], [z1,z1], /DATA, /T3D, COLOR=c1Color, $
            LINESTYLE=1
         PLOTS, [x1,xM], [y1,y1], [z1,z1], /DATA, /T3D, COLOR=c1Color
      endelse

      norm = [[x1,y1,z1], [x1,y1,z1-1]]
      norm = VERT_T3D(norm, /NO_COPY)
      if ((norm(2,1)-norm(2,0)) gt 0.0) then begin
         PLOTS, [x2,x2], [y2,y2], [z2, 0], /DATA, /T3D, COLOR=c2Color
         PLOTS, [x2,x2], [y2,y2], [z2,zM], /DATA, /T3D, COLOR=c2Color, $
            LINESTYLE=1
      endif else begin
         PLOTS, [x2,x2], [y2,y2], [z2, 0], /DATA, /T3D, COLOR=c2Color, $
            LINESTYLE=1
         PLOTS, [x2,x2], [y2,y2], [z2,zM], /DATA, /T3D, COLOR=c2Color
      endelse
      norm = [[x1,y1,z1], [x1,y1-1,z1]]
      norm = VERT_T3D(norm, /NO_COPY)
      if ((norm(2,1)-norm(2,0)) gt 0.0) then begin
         PLOTS, [x2,x2], [y2, 0], [z2,z2], /DATA, /T3D, COLOR=c2Color
         PLOTS, [x2,x2], [y2,yM], [z2,z2], /DATA, /T3D, COLOR=c2Color, $
            LINESTYLE=1
      endif else begin
         PLOTS, [x2,x2], [y2, 0], [z2,z2], /DATA, /T3D, COLOR=c2Color, $
            LINESTYLE=1
         PLOTS, [x2,x2], [y2,yM], [z2,z2], /DATA, /T3D, COLOR=c2Color
      endelse
      norm = [[x1,y1,z1], [x1-1,y1,z1]]
      norm = VERT_T3D(norm, /NO_COPY)
      if ((norm(2,1)-norm(2,0)) gt 0.0) then begin
         PLOTS, [x2, 0], [y2,y2], [z2,z2], /DATA, /T3D, COLOR=c2Color
         PLOTS, [x2,xM], [y2,y2], [z2,z2], /DATA, /T3D, COLOR=c2Color, $
            LINESTYLE=1
      endif else begin
         PLOTS, [x2, 0], [y2,y2], [z2,z2], /DATA, /T3D, COLOR=c2Color, $
            LINESTYLE=1
         PLOTS, [x2,xM], [y2,y2], [z2,z2], /DATA, /T3D, COLOR=c2Color
      endelse

   endelse

   PLOTS, [x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1], [z1,z1,z1,z1,z1], $
      /DATA, /T3D, COLOR=blockColor
   PLOTS, [x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1], [z2,z2,z2,z2,z2], $
      /DATA, /T3D, COLOR=blockColor
   PLOTS, [x1,x1], [y1,y1], [z1,z2], /DATA, /T3D, COLOR=blockColor
   PLOTS, [x2,x2], [y1,y1], [z1,z2], /DATA, /T3D, COLOR=blockColor
   PLOTS, [x2,x2], [y2,y2], [z1,z2], /DATA, /T3D, COLOR=blockColor
   PLOTS, [x1,x1], [y2,y2], [z1,z2], /DATA, /T3D, COLOR=blockColor

   EMPTY

   if not(KEYWORD_SET(directDraw)) then $
      WSET, sMainState.mainWin

end
;******************************************************************************


;******************************************************************************
; Procedure to draw the cube in the small view window.
pro Viz3D_ViewShow, sMainState

   WSET, sMainState.sViewState.viewWin

   ; Draw the cube outline in the current window.
   ERASE
   Viz3D_DrawCube, sMainState, /DIRECT, /AXIS

   WSET, sMainState.mainWin

end
;******************************************************************************


;******************************************************************************
; Procedure to draw the histogram and threshold value in the
; small iso-surface window.
pro Viz3D_SurfShow, sMainState, CALC_HIST=calcHist

   WSET, sMainState.sSurfState.surfWin

   tColor = Viz3D_TransColor(sMainState, sMainState.sSurfState.tColor)
   histColor = Viz3D_TransColor(sMainState, sMainState.sSurfState.histColor)
   axisColor = Viz3D_TransColor(sMainState, sMainState.sSurfState.axisColor)

   saveP = !P
   saveX = !X
   saveY = !Y
   saveZ = !Z

   x1 = 0.03
   y1 = 0.05
   x2 = 0.97
   y2 = 0.95

   lTh = (*(sMainState.sThreshState.hLowThresh))(sMainState.curData)
   hTh = (*(sMainState.sThreshState.hHighThresh))(sMainState.curData)
   surfThresh = *(sMainState.sSurfState.hSurfThresh)
   surfThresh(sMainState.curData) = $
      (surfThresh(sMainState.curData) > lTh) < hTh
   sT = surfThresh(sMainState.curData)
   *(sMainState.sSurfState.hSurfThresh) = surfThresh

   tPos = (sT - lTh) / (hTh - lTh)
   tPos = (tPos * (x2-x1)) + x1

   lRangeHist = *(sMainState.sSurfState.hRangeHist)
   rangeHist = *(lRangeHist(sMainState.curData))

   if (KEYWORD_SET(calcHist)) then begin
      data3D = Viz3D_GetData(sMainState)
      bSize = (hTh - lTh) / 200.0
      n_samples = N_ELEMENTS(data3D)
      IF (n_samples GT 32767L) THEN BEGIN
         n_samples = 32767L + ROUND(SQRT(n_samples - 32767L))
         s_index = RANDOMU(s, 1)
         s_index = RANDOMU(s, n_samples)
         s_index = LONG(TEMPORARY(s_index) * FLOAT(n_samples))
         hist_data = FLOAT(data3D[TEMPORARY(s_index)])
         rangeHist = $
            HISTOGRAM(TEMPORARY(hist_data), MIN=lTh, MAX=hTh, $
               BINSIZE=bSize) > 1.0
      ENDIF ELSE BEGIN
         rangeHist = $
            HISTOGRAM(FLOAT(data3D), MIN=lTh, MAX=hTh, $
               BINSIZE=bSize) > 1.0
      ENDELSE
      rangeHist = ALOG(FLOAT(TEMPORARY(rangeHist))) > 0.0
      Viz3D_PutData, data3D, sMainState

      PLOT, rangeHist, T3D=0, /DATA, XSTYLE=5, YSTYLE=5, /NODATA, $
         POSITION=[x1,y1,x2,y2], MIN_VALUE=1.0
      PLOTS, [x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1], T3D=0, $
         /NORMAL, COLOR=axisColor
      OPLOT, rangeHist, T3D=0, COLOR=histColor, MIN_VALUE=1.0, PSYM=3
   
      WSET, sMainState.sSurfState.surfPix
      DEVICE, COPY=[0, 0, 192, 128, 0, 0, sMainState.sSurfState.surfWin]
      WSET, sMainState.sSurfState.surfWin

   endif else begin
      PLOT, rangeHist, T3D=0, /DATA, XSTYLE=5, YSTYLE=5, /NODATA, $
         POSITION=[x1,y1,x2,y2], MIN_VALUE=1.0, /NOERASE
      DEVICE, COPY=[0, 0, 192, 128, 0, 0, sMainState.sSurfState.surfPix]
   endelse

   PLOTS, [tPos, tPos], [y1, y2], T3D=0, /NORMAL, COLOR=tColor
   EMPTY

   *(lRangeHist(sMainState.curData)) = rangeHist
   *(sMainState.sSurfState.hRangeHist) = lRangeHist

   WIDGET_CONTROL, sMainState.sSurfState.wSurfText, $
      SET_VALUE=STRTRIM(STRING(sT),2)

   !P = TEMPORARY(saveP)
   !X = TEMPORARY(saveX)
   !Y = TEMPORARY(saveY)
   !Z = TEMPORARY(saveZ)

   WSET, sMainState.mainWin
end
;******************************************************************************


;******************************************************************************
; Procedure to draw the histogram and threshold values in the
; small threshold window.
pro Viz3D_ThreshShow, sMainState, DYNAMIC=dynUpdate

   WSET, sMainState.sThreshState.threshWin
   lowColor = Viz3D_TransColor(sMainState, sMainState.sThreshState.lowColor)
   highColor = Viz3D_TransColor(sMainState, sMainState.sThreshState.highColor)
   transColor = $
      Viz3D_TransColor(sMainState, sMainState.sThreshState.transColor)
   histColor = Viz3D_TransColor(sMainState, sMainState.sThreshState.histColor)
   axisColor = Viz3D_TransColor(sMainState, sMainState.sThreshState.axisColor)
   backColor = Viz3D_TransColor(sMainState, 0)

   saveP = !P
   saveX = !X
   saveY = !Y
   saveZ = !Z

   x1 = 0.03
   y1 = 0.05
   x2 = 0.97
   y2 = 0.95

   lTh = (*(sMainState.sThreshState.hLowThresh))(sMainState.curData)
   hTh = (*(sMainState.sThreshState.hHighThresh))(sMainState.curData)
   vTp = (*(sMainState.sThreshState.hTransVal))(sMainState.curData)

   minD = (*(sMainState.sThreshState.hMinData))(sMainState.curData)
   maxD = (*(sMainState.sThreshState.hMaxData))(sMainState.curData)

   lowPos = (lTh - minD) / (maxD - minD)
   lowPos = (lowPos * (x2-x1)) + x1
   highPos = (hTh - minD) / (maxD - minD)
   highPos = (highPos * (x2-x1)) + x1

   transPos = (vTp - minD) / (maxD - minD)
   transPos = (transPos * (x2-x1)) + x1

   lDataHist = *(sMainState.sThreshState.hDataHist)
   DataHist = *(lDataHist(sMainState.curData))

   if (KEYWORD_SET(dynUpdate)) then begin
      PLOT, dataHist, T3D=0, /DATA, XSTYLE=5, YSTYLE=5, /NODATA, $
         POSITION=[x1,y1,x2,y2], MIN_VALUE=1.0, /NOERASE
      DEVICE, COPY=[0, 0, 192, 128, 0, 0, sMainState.sThreshState.threshPix]
   endif else begin
      PLOT, dataHist, T3D=0, /DATA, XSTYLE=5, YSTYLE=5, /NODATA, $
         POSITION=[x1,y1,x2,y2], MIN_VALUE=1.0
      PLOTS, [x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1], T3D=0, /NORMAL, COLOR=axisColor
      OPLOT, dataHist, T3D=0, COLOR=histColor, MIN_VALUE=1.0, PSYM=3

      WSET, sMainState.sThreshState.threshPix
      DEVICE, COPY=[0, 0, 192, 128, 0, 0, sMainState.sThreshState.threshWin]
      WSET, sMainState.sThreshState.threshWin

   endelse

   PLOTS, [transPos, transPos], [y1, y2], T3D=0, /NORMAL, COLOR=transColor
   PLOTS, [transPos, transPos], [y1, y2], T3D=0, /NORMAL, COLOR=backColor, $
          LINESTYLE=1
   PLOTS, [lowPos, lowPos], [y1, y2], T3D=0, /NORMAL, COLOR=lowColor
   PLOTS, [highPos, highPos], [y1, y2], T3D=0, /NORMAL, COLOR=highColor


   WIDGET_CONTROL, sMainState.sThreshState.wThreshLowText, $
      SET_VALUE=STRTRIM(STRING(lTh),2)
   WIDGET_CONTROL, sMainState.sThreshState.wThreshHighText, $
      SET_VALUE=STRTRIM(STRING(hTh),2)
   WIDGET_CONTROL, sMainState.sThreshState.wThreshTransText, $
      SET_VALUE=STRTRIM(STRING(vTp),2)

   !P = TEMPORARY(saveP)
   !X = TEMPORARY(saveX)
   !Y = TEMPORARY(saveY)
   !Z = TEMPORARY(saveZ)

   WSET, sMainState.mainWin
end
;******************************************************************************


;******************************************************************************
; Procedure to plot profile data in the small profile window.
pro Viz3D_ProfShow, sMainState, DIRECT=directDraw

   lineColor = Viz3D_TransColor(sMainState, sMainState.sProfState.lineColor)
   axisColor = Viz3D_TransColor(sMainState, sMainState.sProfState.axisColor)
   markColor1 = Viz3D_TransColor(sMainState, sMainState.sProfState.markColor1)
   markColor2 = Viz3D_TransColor(sMainState, sMainState.sProfState.markColor2)

   data3D = Viz3D_GetData(sMainState)

   if (sMainState.sProfState.profType eq 0) then begin ; Orthogonal profile.
      case sMainState.sProfState.profDir of
         1: begin
               profData = data3D(sMainState.sProfState.x1Ortho, $
                               sMainState.sProfState.y1Ortho, *)
               yData = FINDGEN(sMainState.szData(3))
            end
         2: begin
               profData = data3D(sMainState.sProfState.x1Ortho, *, $
                               sMainState.sProfState.z1Ortho)
               yData = FINDGEN(sMainState.szData(2))
            end
         3: begin
               profData = data3D(*, sMainState.sProfState.y1Ortho, $
                               sMainState.sProfState.z1Ortho)
               yData = FINDGEN(sMainState.szData(1))
            end
      endcase
   endif else begin ; Oblique profile.
      samples = MAX(sMainState.szData(1:3)) * 2
      yData = FINDGEN(samples)
      xInter = yData / FLOAT(samples - 1)
      yInter = xInter
      zInter = xInter
      xInter = ((sMainState.sProfState.x2Obliq - $
                 sMainState.sProfState.x1Obliq) * $
                TEMPORARY(xInter)) + sMainState.sProfState.x1Obliq
      yInter = ((sMainState.sProfState.y2Obliq - $
                 sMainState.sProfState.y1Obliq) * $
                TEMPORARY(yInter)) + sMainState.sProfState.y1Obliq
      zInter = ((sMainState.sProfState.z2Obliq - $
                 sMainState.sProfState.z1Obliq) * $
                TEMPORARY(zInter)) + sMainState.sProfState.z1Obliq
      profData = INTERPOLATE(data3D, TEMPORARY(xInter), $
                    TEMPORARY(yInter), TEMPORARY(zInter))
   endelse

   Viz3D_PutData, data3D, sMainState

   saveP = !P
   saveX = !X
   saveY = !Y
   saveZ = !Z

   maxY = N_ELEMENTS(yData) - 1

   lTh = (*(sMainState.sThreshState.hLowThresh))(sMainState.curData)
   hTh = (*(sMainState.sThreshState.hHighThresh))(sMainState.curData)

   if not(KEYWORD_SET(directDraw)) then begin
      WSET, sMainState.sProfState.profPix
      PLOT, profData, yData, TICKLEN=(1), XSTYLE=1, YSTYLE=5, $
         XRANGE=[lTh, hTh], YRANGE=[0, maxY], XMARGIN=[1,1], YMARGIN=[2,1], $
         /NODATA, COLOR=axisColor, CHARSIZE=0.5
      PLOTS, [lTh,lTh], [0,maxY], /DATA, T3D=0, COLOR=axisColor
      PLOTS, [hTh,hTh], [0,maxY], /DATA, T3D=0, COLOR=axisColor
      PLOTS, [lTh,hTh], [0,0], /DATA, T3D=0, COLOR=markColor1, THICK=2
      PLOTS, [lTh,hTh], [maxY,maxY], /DATA, T3D=0, COLOR=markColor2, THICK=2
   endif
   WSET, sMainState.sProfState.profWin
   DEVICE, COPY=[0, 0, 192, sMainState.sProfState.profWinY, 0, 0, $
                 sMainState.sProfState.profPix]

   PLOT, TEMPORARY(profData), TEMPORARY(yData), TICKLEN=(1), XMARGIN=[1,1], $
      YMARGIN=[2,1], /NOERASE, COLOR=lineColor, XSTYLE=5, YSTYLE=5, $
      XRANGE=[lTh, hTh], YRANGE=[0, maxY], CHARSIZE=0.5, /DATA, T3D=0
   EMPTY

   !P = TEMPORARY(saveP)
   !X = TEMPORARY(saveX)
   !Y = TEMPORARY(saveY)
   !Z = TEMPORARY(saveZ)

   WSET, sMainState.mainWin
end
;******************************************************************************


;******************************************************************************
; Function to return the transparent value given a threshold percent.
; Compensate for the differential shading bands using planeDir.
function Viz3D_TranspValu, sMainState, planeDir
   lTh = (*(sMainState.sThreshState.hLowThresh))(sMainState.curData)
   hTh = (*(sMainState.sThreshState.hHighThresh))(sMainState.curData)
   vTp = (*(sMainState.sThreshState.hTransVal))(sMainState.curData)
   transVal = (vTp - lTh) / (hTh - lTh)
   transVal = (transVal * FLOAT(sMainState.sColorState.nColor-1)) > 0.0

   return, ROUND(transVal) + (planeDir * sMainState.sColorState.nColor)
end
;******************************************************************************


;******************************************************************************
; Add a new graphic to the display list.
function Viz3D_Add_Graphic, sMainState, graphic, NO_COPY=noCopy

   gType = TAG_NAMES(graphic, /STRUCTURE_NAME)
   case gType of
      'VIZ3D_ORTHO_PLANE': begin
         gName = 'Slice: Orthogonal, '
         case graphic.orthoDir of
            1: begin
               gName = TEMPORARY(gName) + 'X, ' + $
               STRTRIM(STRING(ROUND(graphic.orthoPos * $
                  sMainState.sViewState.xMax)), 2)
            end
            2: begin
               gName = TEMPORARY(gName) + 'Y, ' + $
               STRTRIM(STRING(ROUND(graphic.orthoPos * $
                  sMainState.sViewState.yMax)), 2)
            end
            3: begin
               gName = TEMPORARY(gName) + 'Z, ' + $
               STRTRIM(STRING(ROUND(graphic.orthoPos * $
                  sMainState.sViewState.zMax)), 2)
            end
         endcase
      end
      'VIZ3D_OBLIQ_PLANE': begin
         gName = 'Slice: Oblique'
      end
      'VIZ3D_BLOCK': begin
         gName = 'Block: '
         case graphic.blockMode of
            0: gName = TEMPORARY(gName) + 'Subract, ('
            1: gName = TEMPORARY(gName) + 'Add, ('
         endcase
         gName = TEMPORARY(gName) + STRTRIM(STRING(graphic.c1(0)), 2) + ','
         gName = TEMPORARY(gName) + STRTRIM(STRING(graphic.c1(1)), 2) + ','
         gName = TEMPORARY(gName) + STRTRIM(STRING(graphic.c1(2)), 2) + '), ('
         gName = TEMPORARY(gName) + STRTRIM(STRING(graphic.c2(0)), 2) + ','
         gName = TEMPORARY(gName) + STRTRIM(STRING(graphic.c2(1)), 2) + ','
         gName = TEMPORARY(gName) + STRTRIM(STRING(graphic.c2(2)), 2) + ')'
      end
      'VIZ3D_SURF': begin
         gName = 'Surface: '
         case graphic.surfSide of
            0: gName = TEMPORARY(gName) + 'Low, '
            1: gName = TEMPORARY(gName) + 'High, '
         endcase
         gName = TEMPORARY(gName) + STRTRIM(STRING(graphic.surfThresh), 2)
      end
      'VIZ3D_PROJ': begin
         gName = 'Projection: '
         case graphic.projType of
            0: gName = TEMPORARY(gName) + 'Max, '
            1: gName = TEMPORARY(gName) + 'Avg, '
         endcase
         case graphic.projReso of
            0: gName = TEMPORARY(gName) + 'Low'
            1: gName = TEMPORARY(gName) + 'Med'
            2: gName = TEMPORARY(gName) + 'High'
         endcase
      end
   endcase

   wDeleteBttn = WIDGET_BUTTON(sMainState.wDeleteMenu, VALUE=gName, $
                               UVALUE='wDeleteBttn')

   p_graphic = PTR_NEW(graphic, NO_COPY=KEYWORD_SET(noCopy))
   ; Add p_graphic to the list
   if (PTR_VALID(sMainState.hDisplayList)) then begin
      if (N_ELEMENTS(*(sMainState.hDisplayList)) GT 0) then begin
         *(sMainState.hDisplayList) = [*(sMainState.hDisplayList), p_graphic]
      endif else begin
         *(sMainState.hDisplayList) = [p_graphic]
      endelse
   endif
   return, p_graphic

end
;******************************************************************************


;******************************************************************************
; Procedure to draw the data on an orthogonal slicing plane.
pro Viz3D_OrthoPlaneDraw, sMainState, SKIP_ADD=skipAdd, SKIP_DRAW=skipDraw

   WIDGET_CONTROL, sMainState.wStatText, $
      SET_VALUE=('Drawing orthogonal plane ...')

   sliceDir = sMainState.sSliceState.orthoDir
   slicePos = sMainState.sSliceState.orthoPos

   ; Check out the data without making a copy.
   data3D = Viz3D_GetData(sMainState)

   case sliceDir of
      1: begin ; X
         pPos = ROUND(slicePos * sMainState.sViewState.xMax)
         dataPlane = REFORM(data3D(pPos, *, *))
         x = REPLICATE(slicePos, 4)
         y = [0.0,1.0,1.0,0.0]
         z = [0.0,0.0,1.0,1.0]
         pDim1 = sMainState.sViewState.yMax
         pDim2 = sMainState.sViewState.zMax
      end
      2: begin ; Y
         pPos = ROUND(slicePos * sMainState.sViewState.yMax)
         dataPlane = REFORM(data3D(*, pPos, *))
         x = [0.0,1.0,1.0,0.0]
         y = REPLICATE(slicePos, 4)
         z = [0.0,0.0,1.0,1.0]
         pDim1 = sMainState.sViewState.xMax
         pDim2 = sMainState.sViewState.zMax
      end
      3: begin ; Z
         pPos = ROUND(slicePos * sMainState.sViewState.zMax)
         dataPlane = REFORM(data3D(*, *, pPos))
         x = [0.0,1.0,1.0,0.0]
         y = [0.0,0.0,1.0,1.0]
         z = REPLICATE(slicePos, 4)
         pDim1 = sMainState.sViewState.xMax
         pDim2 = sMainState.sViewState.yMax
      end
   endcase

   ; Check the data back in again.
   Viz3D_PutData, data3D, sMainState

   ; Extract the data on the orthogonal plane.
   dataPlane = Viz3D_ScaleData(Temporary(dataPlane), sMainState)

   ; Scale the data to the current threshold.
   dataPlane = BYTSCL(TEMPORARY(dataPlane), MIN=0B, MAX=255B, $
                      TOP=sMainState.sColorState.nColor-1)

   ; Shift the data into the appropriate differential shading
   ; portion of the color table.
   dataPlane = TEMPORARY(dataPlane) + $
      BYTE((sliceDir-1) * sMainState.sColorState.nColor)

   ; Draw the slice in the Z buffer.

   SET_PLOT, 'Z'

   if (sMainState.sSliceState.drawMode) then begin ; Expose mode
      tempImage = TVRD()
      tempDepth = TVRD(/WORDS, CHANNEL=1)
      ERASE
      POLYFILL, x, y, z, PATTERN=dataPlane, /T3D, /NORMAL, $
         IMAGE_INTERP=sMainState.interpMode, $
         TRANSPARENT=Viz3D_TranspValu(sMainState, sliceDir-1), $
         IMAGE_COORD=[[0,0], [pDim1,0], [pDim1,pDim2], [0,pDim2]]
      planeDepth = TVRD(/WORDS, CHANNEL=1)
      planeIndex = WHERE(TEMPORARY(planeDepth) gt (-32765))
      tempDepth(TEMPORARY(planeIndex)) = (-32765)
      planeDepth = 0
      TV, TEMPORARY(tempDepth), /WORDS, CHANNEL=1
      TV, TEMPORARY(tempImage)
   endif

   POLYFILL, x, y, z, PATTERN=dataPlane, /T3D, /NORMAL, $
      IMAGE_INTERP=sMainState.interpMode, $
      TRANSPARENT=Viz3D_TranspValu(sMainState, sliceDir-1), $
      IMAGE_COORD=[[0,0], [pDim1,0], [pDim1,pDim2], [0,pDim2]]

   SET_PLOT, sMainState.screenDevice

   ; Update the big view window.
   if (NOT(KEYWORD_SET(skipDraw))) then $
      Viz3D_DrawData, sMainState

   ; Add orthogonal plane to sMainState.hDisplayList
   if (NOT(KEYWORD_SET(skipAdd))) then begin
      lTh = (*(sMainState.sThreshState.hLowThresh))(sMainState.curData)
      hTh = (*(sMainState.sThreshState.hHighThresh))(sMainState.curData)
      vTp = (*(sMainState.sThreshState.hTransVal))(sMainState.curData)
      graphic = {VIZ3D_ORTHO_PLANE, curData:sMainState.curData, $
                 orthoDir:sMainState.sSliceState.orthoDir, $
                 orthoPos:sMainState.sSliceState.orthoPos, $
                 drawMode:sMainState.sSliceState.drawMode, $
                 lowThresh:lTh, highThresh:hTh, transVal:vTp}
      h = Viz3D_Add_Graphic(sMainState, graphic, /NO_COPY)
   endif

   WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=('')

end
;******************************************************************************


;******************************************************************************
; Procedure to draw the data on an oblique slicing plane.
pro Viz3D_ObliqPlaneDraw, sMainState, vertPlane, vert3D, SKIP_ADD=skipAdd, $
       SKIP_DRAW=skipDraw

   WIDGET_CONTROL, sMainState.wStatText, $
      SET_VALUE=('Drawing oblique plane ...')

   sliceNormal = sMainState.sSliceState.obliqNormal
   sliceCenter = sMainState.sSliceState.obliqCenter

   ; Compute the maximum necessary slicing plane image dimension.
   maxLen = LONG(SQRT(sMainState.sViewState.xMax^2 + $
                      sMainState.sViewState.yMax^2 + $
                      sMainState.sViewState.zMax^2))

   ; Compute the amount necessary to shift the slicing plane to the center.
   maxXvert = MAX(vertPlane(0,*), MIN=minXvert)
   maxYvert = MAX(vertPlane(1,*), MIN=minYvert)
   shiftX = (maxXvert + minXvert) / 2.0
   shiftY = (maxYvert + minYvert) / 2.0

   maxP = SQRT(3.0) ; Slicing rectangle interpolation points
                    ; range between 0 and SQRT(3) in X and Y
                    ; (before shift and transformation).

   ; Compute the shifted interpolation points on the plane.
   ; (Normalized coordinates, with no rotation or translation).
   planePts = FLTARR(3, maxLen^2, /NOZERO)
   sF = FLOAT(maxLen-1L)
   planePts(0, *) = shiftX + $
      ((((FINDGEN(maxLen) # REPLICATE(1.0, maxLen)) / sF) - 0.5) * maxP)
   planePts(1, *) = shiftY + $
      ((((REPLICATE(1.0, maxLen) # FINDGEN(maxLen)) / sF) - 0.5) * maxP)
   planePts(2, *) = 0.0

   ; Scale the surface normal vector.
   zDim = sMainState.sViewState.zMax * sMainState.sViewState.zScale
   maxDim = sMainState.sViewState.xMax > sMainState.sViewState.yMax > zDim
   sliceNormU = sliceNormal
   sliceNormU(0) = sliceNormal(0) * sMainState.sViewState.xMax / maxDim
   sliceNormU(1) = sliceNormal(1) * sMainState.sViewState.yMax / maxDim
   sliceNormU(2) = sliceNormal(2) * zDim / maxDim

   ; Compute the rotations and transformation matrix.

   if ((sliceNormU(0) eq 0.0) and (sliceNormU(1) eq 0.0)) then angZ = 0.0 $
   else angZ = ATAN(sliceNormU(1), sliceNormU(0))

   pDistance = SQRT(sliceNormU(0)^2 + sliceNormU(1)^2)
   if ((pDistance eq 0.0) and (sliceNormU(2) eq 0.0)) then angY = 0.0 $
   else angY = ATAN(sliceNormU(2), pDistance)

   sZ = SIN(angZ)
   cZ = COS(angZ)
   sY = SIN((!PI/2.0)-angY)
   cY = COS((!PI/2.0)-angY)

   ident4 = FLTARR(4, 4)
   ident4([0,5,10,15]) = 1.0

   vTrans = ident4
   ; Y rotation.
   m4X4 = ident4
   m4x4(0,0) = cY
   m4x4(0,2) = (-sY)
   m4x4(2,0) = sY
   m4x4(2,2) = cY
   vTrans = TEMPORARY(vTrans) # m4X4
   ; Z rotation.
   m4X4 = ident4
   m4x4(0,0) = cZ
   m4x4(0,1) = sZ
   m4x4(1,0) = (-sZ)
   m4x4(1,1) = cZ
   vTrans = TEMPORARY(vTrans) # m4X4
   ; Translation to the slice center.
   m4X4 = ident4
   m4X4([3,7,11]) = (+sliceCenter)
   vTrans = TEMPORARY(vTrans) # m4X4

   ; Transform the plane points.
   planePts = VERT_T3D(planePts, MATRIX=vTrans, /NO_COPY)
   planePts(0, *) = planePts(0, *) * Float(sMainState.sViewState.xMax)
   planePts(1, *) = planePts(1, *) * Float(sMainState.sViewState.yMax)
   planePts(2, *) = planePts(2, *) * Float(sMainState.sViewState.zMax)

   ; Check out the data without making a copy.
   data3D = Viz3D_GetData(sMainState)

   ; Interpolate the data on the oblique plane using the plane points,
   dataPlane = $
      Interpolate(data3D, planePts(0,*), planePts(1,*), planePts(2,*))

   ; Check the data back in.
   Viz3D_PutData, data3D, sMainState

   ; Free the memory used by planePts.
   planePts = 0

   ; Scale the oblique plane data using the current thresholds.
   dataPlane = Viz3D_ScaleData(Temporary(dataPlane), sMainState)

   ; Scale the oblique plane data into the appropriate part of the color table
   ; (for differential shading).
   dataPlane = BYTSCL(TEMPORARY(dataPlane), MIN=0B, MAX=255B, $
                      TOP=sMainState.sColorState.nColor-1)
   maxDir = MAX(ABS(sliceNormU), maxInd)
   dataPlane = TEMPORARY(dataPlane) + $
      BYTE(maxInd * sMainState.sColorState.nColor)

   ; Reform dataPlane from a one dimensional array of (maxLen^2) elements
   ; to a 2D array with dimensions (maxLen, maxLen).
   dataPlane = REFORM(TEMPORARY(dataPlane), maxLen, maxLen)

   ; USE the Z buffer.

   SET_PLOT, 'Z'

   ; Compute the image coordinates corresponding to the planar cube/plane
   ; intersection points.

   s_vertPlane = vertPlane(0:1, *)
   s_vertPlane(0,*) = (s_vertPlane(0,*) / maxP) + 0.5 - (shiftX / maxP)
   s_vertPlane(1,*) = (s_vertPlane(1,*) / maxP) + 0.5 - (shiftY / maxP)
   s_vertPlane = TEMPORARY(s_vertPlane) * sF

   if (sMainState.sSliceState.drawMode) then begin ; Expose mode
      tempImage = TVRD()
      tempDepth = TVRD(/WORDS, CHANNEL=1)
      ERASE
      POLYFILL, vert3D, PATTERN=dataPlane, /T3D, /NORMAL, $
         IMAGE_INTERP=sMainState.interpMode, $
         TRANSPARENT=Viz3D_TranspValu(sMainState, maxInd), $
         IMAGE_COORD=s_vertPlane
      planeDepth = TVRD(/WORDS, CHANNEL=1)
      planeIndex = WHERE(planeDepth ge (-32765))
      if (planeIndex(0) ge 0L) then begin
         tempDepth(planeIndex) = planeDepth(planeIndex)
      endif
      planeDepth = 0
      TV, TEMPORARY(tempDepth), /WORDS, CHANNEL=1
      TV, TEMPORARY(tempImage)
   endif

   ; Draw the data on the oblique plane.
   POLYFILL, vert3D, PATTERN=dataPlane, /T3D, /NORMAL, $
      IMAGE_INTERP=sMainState.interpMode, $
      TRANSPARENT=Viz3D_TranspValu(sMainState, maxInd), $
      IMAGE_COORD=s_vertPlane

   SET_PLOT, sMainState.screenDevice

   ; Update the big view window.
   if (NOT(KEYWORD_SET(skipDraw))) then $
      Viz3D_DrawData, sMainState

   ; Add oblique plane to sMainState.hDisplayList
   if (NOT(KEYWORD_SET(skipAdd))) then begin
      lTh = (*(sMainState.sThreshState.hLowThresh))(sMainState.curData)
      hTh = (*(sMainState.sThreshState.hHighThresh))(sMainState.curData)
      vTp = (*(sMainState.sThreshState.hTransVal))(sMainState.curData)
      graphic = {VIZ3D_OBLIQ_PLANE, curData:sMainState.curData, $
                 hVertPlane:PTR_NEW(), hVert3D:PTR_NEW(), $
                 drawMode:sMainState.sSliceState.drawMode, $
                 obliqNormal:sMainState.sSliceState.obliqNormal, $
                 obliqCenter:sMainState.sSliceState.obliqCenter, $
                 lowThresh:lTh, highThresh:hTh, transVal:vTp}
      h = Viz3D_Add_Graphic(sMainState, graphic, /NO_COPY)
      graphic = *h
      graphic.hVertPlane = PTR_NEW(vertPlane, /NO_COPY)
      graphic.hVert3D = PTR_NEW(vert3D, /NO_COPY)
      *h = graphic
   endif

   WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=('')
end
;******************************************************************************


;******************************************************************************
; Procedure to draw the data on a block.
pro Viz3D_BlockDraw, sMainState, SKIP_ADD=skipAdd, SKIP_DRAW=skipDraw

   WIDGET_CONTROL, sMainState.wStatText, $
      SET_VALUE=('Drawing block ...')

   x1 = FLOAT(sMainState.sBlockState.c1(0) < sMainState.sBlockState.c2(0))
   x2 = FLOAT(sMainState.sBlockState.c1(0) > sMainState.sBlockState.c2(0))
   y1 = FLOAT(sMainState.sBlockState.c1(1) < sMainState.sBlockState.c2(1))
   y2 = FLOAT(sMainState.sBlockState.c1(1) > sMainState.sBlockState.c2(1))
   z1 = FLOAT(sMainState.sBlockState.c1(2) < sMainState.sBlockState.c2(2))
   z2 = FLOAT(sMainState.sBlockState.c1(2) > sMainState.sBlockState.c2(2))

   ; Return if any block dimension is zero.
   if (((x1 eq x2) or (y1 eq y2)) or (z1 eq z2)) then return

   ; Check out the data without making a copy.
   data3D = Viz3D_GetData(sMainState)

   SET_PLOT, 'Z'

   if (sMainState.sBlockState.blockMode) then begin ; Add mode.
      dim1 = (x2 - x1)
      dim2 = (y2 - y1)
      faceDir = 2
      subImg = Viz3D_ScaleData(REFORM(data3D(x1:x2, y1:y2, z1)), sMainState)
      subImg = BYTSCL(TEMPORARY(subImg), MIN=0B, MAX=255B, $
                      TOP=sMainState.sColorState.nColor-1)
      subImg = TEMPORARY(subImg) + $
                  BYTE(faceDir * sMainState.sColorState.nColor)
      POLYFILL, [x1,x2,x2,x1], [y1,y1,y2,y2], [z1,z1,z1,z1], PATTERN=subImg, $
         /T3D, /DATA, IMAGE_INTERP=sMainState.interpMode, $
         TRANSPARENT=Viz3D_TranspValu(sMainState, faceDir), $
         IMAGE_COORD=[[0,0],[dim1,0],[dim1,dim2],[0,dim2]]
      subImg = Viz3D_ScaleData(REFORM(data3D(x1:x2, y1:y2, z2)), sMainState)
      subImg = BYTSCL(TEMPORARY(subImg), MIN=0B, MAX=255B, $
                      TOP=sMainState.sColorState.nColor-1)
      subImg = TEMPORARY(subImg) + $
                  BYTE(faceDir * sMainState.sColorState.nColor)
      POLYFILL, [x1,x2,x2,x1], [y1,y1,y2,y2], [z2,z2,z2,z2], PATTERN=subImg, $
         /T3D, /DATA, IMAGE_INTERP=sMainState.interpMode, $
         TRANSPARENT=Viz3D_TranspValu(sMainState, faceDir), $
         IMAGE_COORD=[[0,0],[dim1,0],[dim1,dim2],[0,dim2]]

      dim1 = (x2 - x1)
      dim2 = (z2 - z1)
      faceDir = 1
      subImg = Viz3D_ScaleData(REFORM(data3D(x1:x2, y1, z1:z2)), sMainState)
      subImg = BYTSCL(TEMPORARY(subImg), MIN=0B, MAX=255B, $
                      TOP=sMainState.sColorState.nColor-1)
      subImg = TEMPORARY(subImg) + $
                  BYTE(faceDir * sMainState.sColorState.nColor)
      POLYFILL, [x1,x2,x2,x1], [y1,y1,y1,y1], [z1,z1,z2,z2], PATTERN=subImg, $
         /T3D, /DATA, IMAGE_INTERP=sMainState.interpMode, $
         TRANSPARENT=Viz3D_TranspValu(sMainState, faceDir), $
         IMAGE_COORD=[[0,0],[dim1,0],[dim1,dim2],[0,dim2]]
      subImg = Viz3D_ScaleData(REFORM(data3D(x1:x2, y2, z1:z2)), sMainState)
      subImg = BYTSCL(TEMPORARY(subImg), MIN=0B, MAX=255B, $
                      TOP=sMainState.sColorState.nColor-1)
      subImg = TEMPORARY(subImg) + $
                  BYTE(faceDir * sMainState.sColorState.nColor)
      POLYFILL, [x1,x2,x2,x1], [y2,y2,y2,y2], [z1,z1,z2,z2], PATTERN=subImg, $
         /T3D, /DATA, IMAGE_INTERP=sMainState.interpMode, $
         TRANSPARENT=Viz3D_TranspValu(sMainState, faceDir), $
         IMAGE_COORD=[[0,0],[dim1,0],[dim1,dim2],[0,dim2]]

      dim1 = (y2 - y1)
      dim2 = (z2 - z1)
      faceDir = 0
      subImg = Viz3D_ScaleData(REFORM(data3D(x1, y1:y2, z1:z2)), sMainState)
      subImg = BYTSCL(TEMPORARY(subImg), MIN=0B, MAX=255B, $
                      TOP=sMainState.sColorState.nColor-1)
      subImg = TEMPORARY(subImg) + $
                  BYTE(faceDir * sMainState.sColorState.nColor)
      POLYFILL, [x1,x1,x1,x1], [y1,y2,y2,y1], [z1,z1,z2,z2], PATTERN=subImg, $
         /T3D, /DATA, IMAGE_INTERP=sMainState.interpMode, $
         TRANSPARENT=Viz3D_TranspValu(sMainState, faceDir), $
         IMAGE_COORD=[[0,0],[dim1,0],[dim1,dim2],[0,dim2]]
      subImg = Viz3D_ScaleData(REFORM(data3D(x2, y1:y2, z1:z2)), sMainState)
      subImg = BYTSCL(TEMPORARY(subImg), MIN=0B, MAX=255B, $
                      TOP=sMainState.sColorState.nColor-1)
      subImg = TEMPORARY(subImg) + $
                  BYTE(faceDir * sMainState.sColorState.nColor)
      POLYFILL, [x2,x2,x2,x2], [y1,y2,y2,y1], [z1,z1,z2,z2], PATTERN=subImg, $
         /T3D, /DATA, IMAGE_INTERP=sMainState.interpMode, $
         TRANSPARENT=Viz3D_TranspValu(sMainState, faceDir), $
         IMAGE_COORD=[[0,0],[dim1,0],[dim1,dim2],[0,dim2]]


   endif else begin ; Subtract mode.
      dims = FLOAT([sMainState.sViewState.xMax, sMainState.sViewState.yMax, $
                    sMainState.sViewState.zMax])
      dims = [[dims], [dims]]

      imageBuffer1 = TVRD()
      depthBuffer1 = TVRD(CHANNEL=1, /WORDS)
      POLYFILL, [x1,x2,x2,x1], [y1,y1,y2,y2], [z1,z1,z1,z1], /T3D, /DATA
      POLYFILL, [x1,x2,x2,x1], [y1,y1,y2,y2], [z2,z2,z2,z2], /T3D, /DATA
      POLYFILL, [x1,x2,x2,x1], [y1,y1,y1,y1], [z1,z1,z2,z2], /T3D, /DATA
      POLYFILL, [x1,x2,x2,x1], [y2,y2,y2,y2], [z1,z1,z2,z2], /T3D, /DATA
      POLYFILL, [x1,x1,x1,x1], [y1,y2,y2,y1], [z1,z1,z2,z2], /T3D, /DATA
      POLYFILL, [x2,x2,x2,x2], [y1,y2,y2,y1], [z1,z1,z2,z2], /T3D, /DATA
      depthMask1 = (TVRD(CHANNEL=1, /WORDS) gt depthBuffer1)

      ERASE

      norm = [[x1,y1,z1], [x1,y1,z1-1.0]] / dims
      norm = VERT_T3D(norm, /NO_COPY)
      if ((norm(2,1)-norm(2,0)) lt 0.0) then begin
         dim1 = (x2 - x1)
         dim2 = (y2 - y1)
         faceDir = 2
         subImg = $
            Viz3D_ScaleData(REFORM(data3D(x1:x2, y1:y2, z1)), sMainState)
         subImg = BYTSCL(TEMPORARY(subImg), MIN=0B, MAX=255B, $
                         TOP=sMainState.sColorState.nColor-1)
         subImg = TEMPORARY(subImg) + $
                     BYTE(faceDir * sMainState.sColorState.nColor)
         POLYFILL, [x1,x2,x2,x1], [y1,y1,y2,y2], [z1,z1,z1,z1], $
            PATTERN=subImg, /T3D, /DATA, IMAGE_INTERP=sMainState.interpMode, $
            TRANSPARENT=Viz3D_TranspValu(sMainState, faceDir), $
            IMAGE_COORD=[[0,0],[dim1,0],[dim1,dim2],[0,dim2]]
      endif

      norm = [[x1,y1,z1], [x1,y1-1.0,z1]] / dims
      norm = VERT_T3D(norm, /NO_COPY)
      if ((norm(2,1)-norm(2,0)) lt 0.0) then begin
         dim1 = (x2 - x1)
         dim2 = (z2 - z1)
         faceDir = 1
         subImg = $
            Viz3D_ScaleData(REFORM(data3D(x1:x2, y1, z1:z2)), sMainState)
         subImg = BYTSCL(TEMPORARY(subImg), MIN=0B, MAX=255B, $
                         TOP=sMainState.sColorState.nColor-1)
         subImg = TEMPORARY(subImg) + $
                     BYTE(faceDir * sMainState.sColorState.nColor)
         POLYFILL, [x1,x2,x2,x1], [y1,y1,y1,y1], [z1,z1,z2,z2], $
            PATTERN=subImg, /T3D, /DATA, IMAGE_INTERP=sMainState.interpMode, $
            TRANSPARENT=Viz3D_TranspValu(sMainState, faceDir), $
            IMAGE_COORD=[[0,0],[dim1,0],[dim1,dim2],[0,dim2]]
      endif

      norm = [[x1,y1,z1], [x1-1.0,y1,z1]] / dims
      norm = VERT_T3D(norm, /NO_COPY)
      if ((norm(2,1)-norm(2,0)) lt 0.0) then begin
         dim1 = (y2 - y1)
         dim2 = (z2 - z1)
         faceDir = 0
         subImg = $
            Viz3D_ScaleData(REFORM(data3D(x1, y1:y2, z1:z2)), sMainState)
         subImg = BYTSCL(TEMPORARY(subImg), MIN=0B, MAX=255B, $
                         TOP=sMainState.sColorState.nColor-1)
         subImg = TEMPORARY(subImg) + $
                     BYTE(faceDir * sMainState.sColorState.nColor)
         POLYFILL, [x1,x1,x1,x1], [y1,y2,y2,y1], [z1,z1,z2,z2], $
            PATTERN=subImg, /T3D, /DATA, IMAGE_INTERP=sMainState.interpMode, $
            TRANSPARENT=Viz3D_TranspValu(sMainState, faceDir), $
            IMAGE_COORD=[[0,0],[dim1,0],[dim1,dim2],[0,dim2]]
      endif

      norm = [[x2,y2,z2], [x2,y2,z2+1.0]] / dims
      norm = VERT_T3D(norm, /NO_COPY)
      if ((norm(2,1)-norm(2,0)) lt 0.0) then begin
         dim1 = (x2 - x1)
         dim2 = (y2 - y1)
         faceDir = 2
         subImg = $
            Viz3D_ScaleData(REFORM(data3D(x1:x2, y1:y2, z2)), sMainState)
         subImg = BYTSCL(TEMPORARY(subImg), MIN=0B, MAX=255B, $
                         TOP=sMainState.sColorState.nColor-1)
         subImg = TEMPORARY(subImg) + $
                     BYTE(faceDir * sMainState.sColorState.nColor)
         POLYFILL, [x1,x2,x2,x1], [y1,y1,y2,y2], [z2,z2,z2,z2], $
            PATTERN=subImg, /T3D, /DATA, IMAGE_INTERP=sMainState.interpMode, $
            TRANSPARENT=Viz3D_TranspValu(sMainState, faceDir), $
            IMAGE_COORD=[[0,0],[dim1,0],[dim1,dim2],[0,dim2]]
      endif

      norm = [[x2,y2,z2], [x2,y2+1.0,z2]] / dims
      norm = VERT_T3D(norm, /NO_COPY)
      if ((norm(2,1)-norm(2,0)) lt 0.0) then begin
         dim1 = (x2 - x1)
         dim2 = (z2 - z1)
         faceDir = 1
         subImg = $
            Viz3D_ScaleData(REFORM(data3D(x1:x2, y2, z1:z2)), sMainState)
         subImg = BYTSCL(TEMPORARY(subImg), MIN=0B, MAX=255B, $
                         TOP=sMainState.sColorState.nColor-1)
         subImg = TEMPORARY(subImg) + $
                     BYTE(faceDir * sMainState.sColorState.nColor)
         POLYFILL, [x1,x2,x2,x1], [y2,y2,y2,y2], [z1,z1,z2,z2], $
            PATTERN=subImg, /T3D, /DATA, IMAGE_INTERP=sMainState.interpMode, $
            TRANSPARENT=Viz3D_TranspValu(sMainState, faceDir), $
            IMAGE_COORD=[[0,0],[dim1,0],[dim1,dim2],[0,dim2]]
      endif

      norm = [[x2,y2,z2], [x2+1.0,y2,z2]] / dims
      norm = VERT_T3D(norm, /NO_COPY)
      if ((norm(2,1)-norm(2,0)) lt 0.0) then begin
         dim1 = (y2 - y1)
         dim2 = (z2 - z1)
         faceDir = 0
         subImg = $
            Viz3D_ScaleData(REFORM(data3D(x2, y1:y2, z1:z2)), sMainState)
         subImg = BYTSCL(TEMPORARY(subImg), MIN=0B, MAX=255B, $
                         TOP=sMainState.sColorState.nColor-1)
         subImg = TEMPORARY(subImg) + $
                     BYTE(faceDir * sMainState.sColorState.nColor)
         POLYFILL, [x2,x2,x2,x2], [y1,y2,y2,y1], [z1,z1,z2,z2], $
            PATTERN=subImg, /T3D, /DATA, IMAGE_INTERP=sMainState.interpMode, $
            TRANSPARENT=Viz3D_TranspValu(sMainState, faceDir), $
            IMAGE_COORD=[[0,0],[dim1,0],[dim1,dim2],[0,dim2]]
      endif

      imageBuffer2 = TVRD()
      depthBuffer2 = TVRD(CHANNEL=1, /WORDS)
      depthMask2 = (depthBuffer2 gt (-32765)) and $
                   (depthBuffer2 lt depthBuffer1)
      depthMask = WHERE(TEMPORARY(depthMask1) * TEMPORARY(depthMask2))

      if (depthMask(0) ge 0L) then begin
         imageBuffer1(depthMask) = imageBuffer2(depthMask)
         depthBuffer1(depthMask) = depthBuffer2(depthMask)
      endif

      TV, imageBuffer1
      TV, depthBuffer1, CHANNEL=1, /WORDS

   endelse

   ; Put the data back.
   Viz3D_PutData, data3D, sMainState

   SET_PLOT, sMainState.screenDevice

   ; Update the big view window.
   if (NOT(KEYWORD_SET(skipDraw))) then $
      Viz3D_DrawData, sMainState

   ; Add block to sMainState.hDisplayList
   if (NOT(KEYWORD_SET(skipAdd))) then begin
      lTh = (*(sMainState.sThreshState.hLowThresh))(sMainState.curData)
      hTh = (*(sMainState.sThreshState.hHighThresh))(sMainState.curData)
      vTp = (*(sMainState.sThreshState.hTransVal))(sMainState.curData)
      graphic = {VIZ3D_BLOCK, curData:sMainState.curData, $
                 c1:sMainState.sBlockState.c1, $
                 c2:sMainState.sBlockState.c2, $
                 blockMode:sMainState.sBlockState.blockMode, $
                 lowThresh:lTh, highThresh:hTh, transVal:vTp}
      h = Viz3D_Add_Graphic(sMainState, graphic, /NO_COPY)
   endif

   WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=('')
end
;******************************************************************************


;******************************************************************************
; Procedure to draw an iso-surface.
pro Viz3D_SurfDraw, sMainState, hVertList, hPolyList, hColorList, $
       SKIP_ADD=skipAdd, SKIP_DRAW=skipDraw, CALC_SURF=calcSurf

   if (KEYWORD_SET(calcSurf)) then begin
      WIDGET_CONTROL, sMainState.wStatText, $
         SET_VALUE=('Calculating iso-surface ...')

      ; Check out the data without making a copy.
      data3D = Viz3D_GetData(sMainState)

      surfThresh = *(sMainState.sSurfState.hSurfThresh)
      sT = surfThresh(sMainState.curData)
      SHADE_VOLUME, data3D, sT, vertList, polyList, $
                    LOW=(sMainState.sSurfState.surfSide eq 0)
      Viz3D_PutData, data3D, sMainState

      if (N_ELEMENTS(polyList) le 0L) then begin ; No iso-surface to render.
         WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=''
         return
      endif

      if (polyList(0) lt 0L) then begin ; No iso-surface to render.
         WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=''
         return
      endif

      if (sMainState.curData eq sMainState.sSurfState.curShade) then begin
         ; Light-source shading.
         colorList = (-1L)
      endif else begin
         data3D = $
            Viz3D_GetData(sMainState, sMainState.sSurfState.curShade)
         shadeData = data3D
         Viz3D_PutData, data3D, sMainState, sMainState.sSurfState.curShade
         colorList = INTERPOLATE(TEMPORARY(shadeData), $
            vertList(0,*), vertList(1,*), vertList(2,*), MISSING=0)
         lowThresh = *(sMainState.sThreshState.hLowThresh)
         highThresh = *(sMainState.sThreshState.hHighThresh)
         lTh = lowThresh(sMainState.sSurfState.curShade)
         hTh = highThresh(sMainState.sSurfState.curShade)
         colorList = BYTSCL(TEMPORARY(colorList), MIN=lTh, MAX=hTh, $
            TOP=sMainState.sColorState.nColor-1)
         colorList = TEMPORARY(colorList) + $
            BYTE(3 * sMainState.sColorState.nColor)
      endelse

      WIDGET_CONTROL, sMainState.wStatText, $
         SET_VALUE=(STRTRIM((SIZE(vertList))(2),2) + ' Verts, ' + $
         STRTRIM((SIZE(polyList))(1)/4,2) + ' Polys.')
      hVertList = PTR_NEW(vertList, /NO_COPY)
      hPolyList = PTR_NEW(polyList, /NO_COPY)
      hColorList = PTR_NEW(colorList, /NO_COPY)
   endif else begin
      WIDGET_CONTROL, sMainState.wStatText, $
         SET_VALUE=('Drawing iso-surface ...')
   endelse

   SET_PLOT, 'Z'

   SET_SHADING, /GOURAUD, LIGHT=[-0.5, 0.5, 0.5], REJECT=0, $
      VALUES=[(sMainState.sColorState.nColor*3), $
              (sMainState.sColorState.nColor*5)-1]

   if ((*hColorList)(0) lt 0B) then begin ; Light-source shading.
      nothing = POLYSHADE(*hVertList, *hPolyList, /DATA, /T3D, $
         TOP=((sMainState.sColorState.nColor*4)-1))
   endif else begin ; Data shading.
      nothing = POLYSHADE(*hVertList, *hPolyList, /DATA, /T3D, $
         SHADES=(*hColorList), TOP=((sMainState.sColorState.nColor*4)-1))
   endelse


   SET_PLOT, sMainState.screenDevice

   ; Update the big view window.
   if (NOT(KEYWORD_SET(skipDraw))) then $
      Viz3D_DrawData, sMainState

   ; Add surface to sMainState.hDisplayList
   if (NOT(KEYWORD_SET(skipAdd))) then begin
      sT = (*(sMainState.sSurfState.hSurfThresh))(sMainState.curData)
      graphic = {VIZ3D_SURF, curData:sMainState.curData, $
                 hVertList:PTR_NEW(/ALLOCATE_HEAP), $
                 hPolyList:PTR_NEW(/ALLOCATE_HEAP), $
                 hColorList:PTR_NEW(/ALLOCATE_HEAP), $
                 surfThresh:sT, surfSide:sMainState.sSurfState.surfSide}
      h = Viz3D_Add_Graphic(sMainState, graphic, /NO_COPY)
      if (KEYWORD_SET(calcSurf)) then begin
         graphic = *h
         *((*h).hVertList) = TEMPORARY(*hVertList)
         *((*h).hPolyList) = TEMPORARY(*hPolyList)
         *((*h).hColorList) = TEMPORARY(*hColorList)
         PTR_FREE, hVertList
         PTR_FREE, hPolyList
         PTR_FREE, hColorList
      endif
   endif

   WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=('')
end
;******************************************************************************


;******************************************************************************
; Procedure to draw a projection.
pro Viz3D_ProjDraw, sMainState, SKIP_ADD=skipAdd, SKIP_DRAW=skipDraw

   WIDGET_CONTROL, sMainState.wStatText, $
      SET_VALUE=('Calculating projection ...')

   SET_PLOT, 'Z'
   imageBuffer = TVRD()
   depthBuffer = TVRD(CHANNEL=1, /WORDS)

   annoIndex = WHERE(imageBuffer gt (5*sMainState.sColorState.nColor))
   if (annoIndex(0) ge 0L) then depthBuffer(TEMPORARY(annoIndex)) = (-32765)

   minDim = sMainState.sViewState.xMax < sMainState.sViewState.yMax < $
            sMainState.sViewState.zMax
   maxDim = sMainState.sViewState.xMax > sMainState.sViewState.yMax > $
            sMainState.sViewState.zMax
   avgDim = Float(maxDim + minDim) / 2.0
   case sMainState.sProjState.projReso of
      0: begin ; Low.
         xS = ROUND(FLOAT(avgDim) / 2.0) > 2
         yS = ROUND(FLOAT(avgDim) / 2.0) > 2
         zS = ROUND(FLOAT(avgDim) / 2.0) > 2
      end
      1: begin ; Medium.
         xS = ROUND(FLOAT(avgDim) / 1.5) > 2
         yS = ROUND(FLOAT(avgDim) / 1.5) > 2
         zS = ROUND(FLOAT(avgDim) / 1.5) > 2
      end
      2: begin ; high.
         xS = ROUND(FLOAT(maxDim) / 1.0) > 2
         yS = ROUND(FLOAT(maxDim) / 1.0) > 2
         zS = ROUND(FLOAT(maxDim) / 1.0) > 2
      end
   endcase

   ; Calculate the projection.

   data3D_raw = Viz3D_GetData(sMainState)
   data3D = data3D_raw
   Viz3D_PutData, data3D_raw, sMainState
   data3D = Viz3D_ScaleData(data3D, sMainState)

   xyS = xS * yS

   cubePts = FLTARR(8,4)
   cubePts(0,*) = [0,0,0,1]
   cubePts(1,*) = [1,0,0,1]
   cubePts(2,*) = [1,1,0,1]
   cubePts(3,*) = [0,1,0,1]
   cubePts(4,*) = [0,0,1,1]
   cubePts(5,*) = [1,0,1,1]
   cubePts(6,*) = [1,1,1,1]
   cubePts(7,*) = [0,1,1,1]
   cubePts = cubepts # !P.T
   xMin = MIN(cubePts(*,0), MAX=xMax)
   yMin = MIN(cubePts(*,1), MAX=yMax)
   zMin = MIN(cubePts(*,2), MAX=zMax)
   zRange = zMax - zMin

   xInd = REFORM(((Findgen(xS) / Float(xS - 1L)) # Replicate(1.0, yS)), xyS)
   yInd = REFORM((Replicate(1.0, xS) # (Findgen(yS) / Float(yS - 1L))), xyS)

   index = FLTARR(xyS, 4, /NOZERO)
   index(0, 0) = TEMPORARY(xInd)
   index(0, 1) = TEMPORARY(yInd)
   index(*, 3) = 1.0

   kMax = zS - 1
   fkMax = FLOAT(kMax)

   pDepth = CONGRID(depthBuffer, Xs, Ys, /INTERP, /MINUS_ONE)
   if (sMainState.sProjState.projType) then $
      sumImage = FLTARR(xS, yS) $
   else sumImage = BYTARR(xS, yS)

   depthQ = (100.0 - sMainState.sProjState.depthQ) / 100.0
   for k=0, kMax do begin ; Project one plane at a time.

      fk = FLOAT(k)
      zVal = ((fk / fkMax) * zRange) + zMin
      index(*, 2) = zVal
      tIndex = index # sMainState.sViewState.invTrans

      indX = ((tIndex(*, 0) / tIndex(*, 3)) - !X.S(0)) / !X.S(1)
      indY = ((tIndex(*, 1) / tIndex(*, 3)) - !Y.S(0)) / !Y.S(1)
      indZ = ((tIndex(*, 2) / tIndex(*, 3)) - !Z.S(0)) / !Z.S(1)

      zPlane = ROUND(2.0 * ((((zVal - 0.5) > (-0.5)) < 0.5) * 32765.0))
      zInd = WHERE(pDepth lt zPlane)

      if (zInd(0) ge 0L) then begin ; Something new in front
         depthFac = BYTE((1.0 - ((zVal-zMin)/(zMax-zMin))) * depthQ * 255.0)
         tempImg = REFORM(((INTERPOLATE(data3D, indX, indY, indZ, $
                            MISSING=0B) > depthFac) - depthFac), xS, yS)

         lTh = (*(sMainState.sThreshState.hLowThresh))(sMainState.curData)
         hTh = (*(sMainState.sThreshState.hHighThresh))(sMainState.curData)
         vTp = (*(sMainState.sThreshState.hTransVal))(sMainState.curData)

         transVal = (FLOAT(vTp - lTh) / FLOAT(hTh - lTh)) * 255.0
         transVal = BYTE((transVal > 0.0) < 255.0)
         transInd = WHERE(tempImg lt transVal)
         if (transInd(0) ge 0L) then tempImg(TEMPORARY(transInd)) = 0B
         if (sMainState.sProjState.projType) then begin ; Avg projection.
            newVal = FLOAT(tempImg(zInd))
            setInd = WHERE(sumImage(zInd) lt newVal)
            sumImage(zInd) = sumImage(zInd) + (TEMPORARY(newVal) / fkMax)
            if (setInd(0) ge 0L) then begin
               setInd = zInd(TEMPORARY(setInd))
               pDepth(TEMPORARY(setInd)) = zPlane
            endif
         endif else begin ; Max projection.
            setInd = WHERE(sumImage(zInd) lt tempImg(zInd))
            if (setInd(0) ge 0L) then begin
               setInd = zInd(TEMPORARY(setInd))
               pDepth(TEMPORARY(setInd)) = zPlane
            endif
            sumImage = TEMPORARY(sumImage) > tempImg
         endelse
      endif
   endfor
   if (sMainState.sProjState.projType) then begin ; Scale the avg projection.
      minSum = MIN(sumImage, MAX=maxSum)
      if ((maxSum - minSum) gt 0.0) then $
         sumImage = (TEMPORARY(sumImage) - minSum) / (maxSum - minSum)
      minD = (*(sMainState.sThreshState.hMinData))(sMainState.curData)
      maxD = (*(sMainState.sThreshState.hMaxData))(sMainState.curData)
      sumImage = (TEMPORARY(sumImage) * (maxD - minD)) + minD
      sumImage = BYTSCL(TEMPORARY(sumImage), MIN=minD, MAX=maxD, $
         TOP=sMainState.sColorState.nColor-1)
   endif else begin
      sumImage = BYTSCL(TEMPORARY(sumImage), $
         MIN=0, MAX=255, TOP=sMainState.sColorState.nColor-1)
   endelse

   tempImg = 0
   index = 0
   tIndex = 0
   setInd = 0
   indX = 0
   indY = 0
   indZ = 0
   zInd = 0

   sBand = (FLOAT(!D.X_Size) / FLOAT(xS)) > (FLOAT(!D.Y_Size) / FLOAT(yS))
   sBand = ROUND(sBand / 2.0) > 3
   sumImage = CONGRID(TEMPORARY(sumImage), !D.X_Size, !D.Y_Size, $
                      /INTERP, /MINUS_ONE)
   sumImage = SMOOTH(TEMPORARY(sumImage), sBand, /EDGE)
   pDepth = CONGRID(TEMPORARY(pDepth), !D.X_Size, !D.Y_Size, $
                    /INTERP, /MINUS_ONE)
   pDepth = SMOOTH(TEMPORARY(pDepth), sBand, /EDGE)
   newInd = WHERE(pDepth ge depthBuffer)

   if (newInd(0) ge 0L) then begin
      sumImage = TEMPORARY(sumImage) + BYTE(sMainState.sColorState.nColor*4)
      sumImage = sumImage(newInd)
      pDepth = pDepth(newInd)
      imageBuffer(newInd) = TEMPORARY(sumImage)
      depthBuffer(TEMPORARY(newInd)) = TEMPORARY(pDepth)

      ERASE
      TV, TEMPORARY(imageBuffer)
      TV, TEMPORARY(depthBuffer), CHANNEL=1, /WORDS
   endif

   SET_PLOT, sMainState.screenDevice

   ; Update the big view window.
   if (NOT(KEYWORD_SET(skipDraw))) then $
      Viz3D_DrawData, sMainState

   ; Add projection to sMainState.hDisplayList
   if (NOT(KEYWORD_SET(skipAdd))) then begin
      lTh = (*(sMainState.sThreshState.hLowThresh))(sMainState.curData)
      hTh = (*(sMainState.sThreshState.hHighThresh))(sMainState.curData)
      vTp = (*(sMainState.sThreshState.hTransVal))(sMainState.curData)
      graphic = {VIZ3D_PROJ, curData:sMainState.curData, $
                 projType:sMainState.sProjState.projType, $
                 projReso:sMainState.sProjState.projReso, $
                 depthQ:sMainState.sProjState.depthQ, $
                 lowThresh:lTh, highThresh:hTh, transVal:vTp}
      h = Viz3D_Add_Graphic(sMainState, graphic, /NO_COPY)
   endif

   WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=('')
end
;******************************************************************************


;******************************************************************************
; Procedure to draw the probe in the big window.
pro Viz3D_ProbeDraw, sMainState

   lineColor = Viz3D_TransColor(sMainState, sMainState.sProbeState.lineColor)

   WSET, sMainState.mainWin
   DEVICE, COPY=[0, 0, sMainState.winX, sMainState.winY, $
                 0, 0, sMainState.pixWin]

   PLOTS, [sMainState.sProbeState.x, sMainState.sProbeState.x], $
          [sMainState.sProbeState.y, sMainState.sProbeState.y], $
          [0, sMainState.szData(3)-1], $
          /DATA, /T3D, COLOR=lineColor

   PLOTS, [sMainState.sProbeState.x, sMainState.sProbeState.x], $
          [0, sMainState.szData(2)-1], $
          [sMainState.sProbeState.z, sMainState.sProbeState.z], $
          /DATA, /T3D, COLOR=lineColor

   PLOTS, [0, sMainState.szData(1)-1], $
          [sMainState.sProbeState.y, sMainState.sProbeState.y], $
          [sMainState.sProbeState.z, sMainState.sProbeState.z], $
          /DATA, /T3D, COLOR=lineColor
   EMPTY

   data3D = Viz3D_GetData(sMainState)
   dataPt = INTERPOLATE(data3D, sMainState.sProbeState.x, $
                                sMainState.sProbeState.y, $
                                sMainState.sProbeState.z)
   Viz3D_PutData, data3D, sMainState

   WIDGET_CONTROL, sMainState.wStatText, $
      SET_VALUE=('Data: ' + STRTRIM(STRING(FLOAT(dataPt)), 2))
   WIDGET_CONTROL, sMainState.sProbeState.wProbeXText, $
      SET_VALUE=STRTRIM(STRING(sMainState.sProbeState.x), 2)
   WIDGET_CONTROL, sMainState.sProbeState.wProbeYText, $
      SET_VALUE=STRTRIM(STRING(sMainState.sProbeState.y), 2)
   WIDGET_CONTROL, sMainState.sProbeState.wProbeZText, $
      SET_VALUE=STRTRIM(STRING(sMainState.sProbeState.z), 2)

end
;******************************************************************************


;******************************************************************************
; Procedure to draw the profile line in the big window.
pro Viz3D_ProfDraw, sMainState

   lineColor = Viz3D_TransColor(sMainState, sMainState.sProfState.lineColor)
   markColor1 = Viz3D_TransColor(sMainState, sMainState.sProfState.markColor1)
   markColor2 = Viz3D_TransColor(sMainState, sMainState.sProfState.markColor2)

   WSET, sMainState.mainWin
   DEVICE, COPY=[0, 0, sMainState.winX, sMainState.winY, $
                 0, 0, sMainState.pixWin]
   if (sMainState.sProfState.profType) then begin
      if ((sMainState.sProfState.z2Obliq ne 0) and $
          (sMainState.sProfState.z2Obliq ne sMainState.szData(3)-1)) then $
         PLOTS, [sMainState.sProfState.x2Obliq, $
                 sMainState.sProfState.x2Obliq], $
                [sMainState.sProfState.y2Obliq, $
                 sMainState.sProfState.y2Obliq], $
                [0, sMainState.szData(3)-1], $
                /DATA, /T3D, COLOR=markColor2
      if ((sMainState.sProfState.y2Obliq ne 0) and $
          (sMainState.sProfState.y2Obliq ne sMainState.szData(2)-1)) then $
         PLOTS, [sMainState.sProfState.x2Obliq, $
                 sMainState.sProfState.x2Obliq], $
                [0, sMainState.szData(2)-1], $
                [sMainState.sProfState.z2Obliq, $
                 sMainState.sProfState.z2Obliq], $
                /DATA, /T3D, COLOR=markColor2
      if ((sMainState.sProfState.x2Obliq ne 0) and $
          (sMainState.sProfState.x2Obliq ne sMainState.szData(1)-1)) then $
         PLOTS, [0, sMainState.szData(1)-1], $
                [sMainState.sProfState.y2Obliq, $
                 sMainState.sProfState.y2Obliq], $
                [sMainState.sProfState.z2Obliq, $
                 sMainState.sProfState.z2Obliq], $
                /DATA, /T3D, COLOR=markColor2

      if ((sMainState.sProfState.z1Obliq ne 0) and $
          (sMainState.sProfState.z1Obliq ne sMainState.szData(3)-1)) then $
         PLOTS, [sMainState.sProfState.x1Obliq, $
                 sMainState.sProfState.x1Obliq], $
                [sMainState.sProfState.y1Obliq, $
                 sMainState.sProfState.y1Obliq], $
                [0, sMainState.szData(3)-1], $
                /DATA, /T3D, COLOR=markColor1
      if ((sMainState.sProfState.y1Obliq ne 0) and $
          (sMainState.sProfState.y1Obliq ne sMainState.szData(2)-1)) then $
         PLOTS, [sMainState.sProfState.x1Obliq, $
                 sMainState.sProfState.x1Obliq], $
                [0, sMainState.szData(2)-1], $
                [sMainState.sProfState.z1Obliq, $
                 sMainState.sProfState.z1Obliq], $
                /DATA, /T3D, COLOR=markColor1
      if ((sMainState.sProfState.x1Obliq ne 0) and $
          (sMainState.sProfState.x1Obliq ne sMainState.szData(1)-1)) then $
         PLOTS, [0, sMainState.szData(1)-1], $
                [sMainState.sProfState.y1Obliq, $
                 sMainState.sProfState.y1Obliq], $
                [sMainState.sProfState.z1Obliq, $
                 sMainState.sProfState.z1Obliq], $
                /DATA, /T3D, COLOR=markColor1

      PLOTS, [sMainState.sProfState.x1Obliq, sMainState.sProfState.x2Obliq], $
             [sMainState.sProfState.y1Obliq, sMainState.sProfState.y2Obliq], $
             [sMainState.sProfState.z1Obliq, sMainState.sProfState.z2Obliq], $
             /DATA, /T3D, COLOR=lineColor

      WIDGET_CONTROL, sMainState.wStatText, $
         SET_VALUE=('(' + $
            STRTRIM(STRING(FIX(sMainState.sProfState.x1Obliq)), 2) + ', ' + $
            STRTRIM(STRING(FIX(sMainState.sProfState.y1Obliq)), 2) + ', ' + $
            STRTRIM(STRING(FIX(sMainState.sProfState.z1Obliq)), 2) + $
               '), (' + $
            STRTRIM(STRING(FIX(sMainState.sProfState.x2Obliq)), 2) + ', ' + $
            STRTRIM(STRING(FIX(sMainState.sProfState.y2Obliq)), 2) + ', ' + $
            STRTRIM(STRING(FIX(sMainState.sProfState.z2Obliq)), 2) + ')')
   endif else begin
      PLOTS, [sMainState.sProfState.x1Ortho, sMainState.sProfState.x2Ortho], $
             [sMainState.sProfState.y1Ortho, sMainState.sProfState.y2Ortho], $
             [sMainState.sProfState.z1Ortho, sMainState.sProfState.z2Ortho], $
             /DATA, /T3D, COLOR=lineColor
      PLOTS, [sMainState.sProfState.x1Ortho], $
             [sMainState.sProfState.y1Ortho], $
             [sMainState.sProfState.z1Ortho], $
             /DATA, /T3D, COLOR=markColor1, PSYM=4, THICK=2
      PLOTS, [sMainState.sProfState.x2Ortho], $
             [sMainState.sProfState.y2Ortho], $
             [sMainState.sProfState.z2Ortho], $
             /DATA, /T3D, COLOR=markColor2, PSYM=4, THICK=2
      WIDGET_CONTROL, sMainState.wStatText, $
         SET_VALUE=('(' + $
            STRTRIM(STRING(FIX(sMainState.sProfState.x1Ortho)), 2) + ', ' + $
            STRTRIM(STRING(FIX(sMainState.sProfState.y1Ortho)), 2) + ', ' + $
            STRTRIM(STRING(FIX(sMainState.sProfState.z1Ortho)), 2) + $
               '), (' + $
            STRTRIM(STRING(FIX(sMainState.sProfState.x2Ortho)), 2) + ', ' + $
            STRTRIM(STRING(FIX(sMainState.sProfState.y2Ortho)), 2) + ', ' + $
            STRTRIM(STRING(FIX(sMainState.sProfState.z2Ortho)), 2) + ')')
   endelse

   EMPTY

end
;******************************************************************************


;******************************************************************************
; Redraw everything in the display list.
pro Viz3D_Redraw, sMainState

   SET_PLOT, 'Z'
   ERASE

   displayListIndex = 0
   if (PTR_VALID(sMainState.hDisplayList)) then begin
      if (N_ELEMENTS(*(sMainState.hDisplayList)) GT 0) then begin
         hID = (*(sMainState.hDisplayList))[0]  ; get first in list
      endif  else hID = PTR_NEW() ;null pointer, not valid
   endif else hID = PTR_NEW() ;null pointer, not valid
   while (PTR_VALID(hID)) do begin
      graphic = *hID
      gType = TAG_NAMES(graphic, /STRUCTURE_NAME)
      case gType of
         'VIZ3D_ORTHO_PLANE': begin
            altState = sMainState
            altState.curData = graphic.curData
            altState.sSliceState.orthoDir = graphic.orthoDir
            altState.sSliceState.orthoPos = graphic.orthoPos
            altState.sSliceState.drawMode = graphic.drawMode

            lowThresh = *(altState.sThreshState.hLowThresh)
            lowThresh(sMainState.curData) = graphic.lowThresh
            hLowThresh = PTR_NEW(lowThresh, /NO_COPY)
            altState.sThreshState.hLowThresh = hLowThresh

            highThresh = *(altState.sThreshState.hHighThresh)
            highThresh(sMainState.curData) = graphic.highThresh
            hHighThresh = PTR_NEW(highThresh, /NO_COPY)
            altState.sThreshState.hHighThresh = hHighThresh

            transVal = *(altState.sThreshState.hTransVal)
            transVal(sMainState.curData) = graphic.transVal
            hTransVal = PTR_NEW(transVal, /NO_COPY)
            altState.sThreshState.hTransVal = hTransVal

            Viz3D_OrthoPlaneDraw, TEMPORARY(altState), /SKIP_ADD, $
               /SKIP_DRAW

            PTR_FREE, hLowThresh
            PTR_FREE, hHighThresh
            PTR_FREE, hTransVal
         end
         'VIZ3D_OBLIQ_PLANE': begin
            altState = sMainState
            altState.curData = graphic.curData
            altState.sSliceState.drawMode = graphic.drawMode
            altState.sSliceState.obliqNormal = graphic.obliqNormal
            altState.sSliceState.obliqCenter = graphic.obliqCenter

            lowThresh = *(altState.sThreshState.hLowThresh)
            lowThresh(sMainState.curData) = graphic.lowThresh
            hLowThresh = PTR_NEW(lowThresh, /NO_COPY)
            altState.sThreshState.hLowThresh = hLowThresh

            highThresh = *(altState.sThreshState.hHighThresh)
            highThresh(sMainState.curData) = graphic.highThresh
            hHighThresh = PTR_NEW(highThresh, /NO_COPY)
            altState.sThreshState.hHighThresh = hHighThresh

            transVal = *(altState.sThreshState.hTransVal)
            transVal(sMainState.curData) = graphic.transVal
            hTransVal = PTR_NEW(transVal, /NO_COPY)
            altState.sThreshState.hTransVal = hTransVal

            vertPlane = *(graphic.hVertPlane)
            vert3D = *(graphic.hVert3D)
            Viz3D_ObliqPlaneDraw, TEMPORARY(altState), vertPlane, vert3D, $
               /SKIP_ADD, /SKIP_DRAW

            PTR_FREE, hLowThresh
            PTR_FREE, hHighThresh
            PTR_FREE, hTransVal
         end
         'VIZ3D_BLOCK': begin
            altState = sMainState
            altState.curData = graphic.curData
            altState.sBlockState.c1 = graphic.c1
            altState.sBlockState.c2 = graphic.c2
            altState.sBlockState.blockMode = graphic.blockMode

            lowThresh = *(altState.sThreshState.hLowThresh)
            lowThresh(sMainState.curData) = graphic.lowThresh
            hLowThresh = PTR_NEW(lowThresh, /NO_COPY)
            altState.sThreshState.hLowThresh = hLowThresh

            highThresh = *(altState.sThreshState.hHighThresh)
            highThresh(sMainState.curData) = graphic.highThresh
            hHighThresh = PTR_NEW(highThresh, /NO_COPY)
            altState.sThreshState.hHighThresh = hHighThresh

            transVal = *(altState.sThreshState.hTransVal)
            transVal(sMainState.curData) = graphic.transVal
            hTransVal = PTR_NEW(transVal, /NO_COPY)
            altState.sThreshState.hTransVal = hTransVal

            Viz3D_BlockDraw, TEMPORARY(altState), /SKIP_ADD, $
               /SKIP_DRAW

            PTR_FREE, hLowThresh
            PTR_FREE, hHighThresh
            PTR_FREE, hTransVal
         end
         'VIZ3D_SURF': begin
            altState = sMainState
            altState.curData = graphic.curData
            altState.sSurfState.surfSide = graphic.surfSide

            surfThresh = *(altState.sSurfState.hSurfThresh)
            surfThresh(sMainState.curData) = graphic.surfThresh
            hSurfThresh = PTR_NEW(surfThresh, /NO_COPY)
            altState.sSurfState.hSurfThresh = hSurfThresh

            Viz3D_SurfDraw, TEMPORARY(altState), graphic.hVertList, $
               graphic.hPolyList, graphic.hColorList, /SKIP_ADD, /SKIP_DRAW

            PTR_FREE, hSurfThresh
         end
         'VIZ3D_PROJ': begin
            altState = sMainState
            altState.curData = graphic.curData
            altState.sProjState.projReso = graphic.projReso
            altState.sProjState.projType = graphic.projType
            altState.sProjState.depthQ = graphic.depthQ

            lowThresh = *(altState.sThreshState.hLowThresh)
            lowThresh(sMainState.curData) = graphic.lowThresh
            hLowThresh = PTR_NEW(lowThresh, /NO_COPY)
            altState.sThreshState.hLowThresh = hLowThresh

            highThresh = *(altState.sThreshState.hHighThresh)
            highThresh(sMainState.curData) = graphic.highThresh
            hHighThresh = PTR_NEW(highThresh, /NO_COPY)
            altState.sThreshState.hHighThresh = hHighThresh

            transVal = *(altState.sThreshState.hTransVal)
            transVal(sMainState.curData) = graphic.transVal
            hTransVal = PTR_NEW(transVal, /NO_COPY)
            altState.sThreshState.hTransVal = hTransVal

            Viz3D_ProjDraw, TEMPORARY(altState), /SKIP_ADD, /SKIP_DRAW

            PTR_FREE, hLowThresh
            PTR_FREE, hHighThresh
            PTR_FREE, hTransVal
         end
      endcase
      *hID = graphic
      displayListIndex = displayListIndex + 1
      if (displayListIndex LT N_ELEMENTS(*(sMainState.hDisplayList))) then $
         hID = (*(sMainState.hDisplayList))[displayListIndex] $
      else hID = PTR_NEW()      ; null since end of list reached
   endwhile

   Viz3D_DrawData, sMainState

end
;******************************************************************************


;******************************************************************************
; Return an xyz data coordinate on the face of the cube.
; This procedure uses the buffers saved in sMainState.frontDepth,
; sMainState.backDepth, sMainState.frontImage, and sMainState.backImage
; (the Viz3D_FillDepth procedure computes them).
; The ON_FACE keyword returns the face that (dX, dY) maps to :
;   0- not on any face,  1- XY,   2- XZ,   3- YZ.
; (dX, dY) is the 2D screen coordinate (where the user clicked).
; invTrans is the inverse of the current view transformation
; (from sMainState.sViewState.invTrans).
function Viz3D_XYZPos, dX, dY, hDepth, hImage, invTrans, ON_FACE=onFace


   szDepth = SIZE(*hDepth)
   dX = (dX > 0) < (szDepth(1)-1)
   dY = (dY > 0) < (szDepth(2)-1)
   dZ = (*hDepth)(dX, dY)
   onFace = (*hImage)(dX, dY)


   zMax = (2L^15L) - 3L
   zMin = (-zmax)

   if (dZ le (zMin)) then return, [0.0, 0.0, 0.0] ; Not on any face.

   nX = Float(dX) / Float(!D.X_Size - 1)
   nY = Float(dY) / Float(!D.Y_Size - 1)
   nZ = (Float(dZ) - Float(zMin)) / Float(zMax - zMin)

   xyzNorm = [nX, nY, nZ, 1.0] # invTrans
   xyzNorm = xyzNorm / xyzNorm(3)

   xPos = (xyzNorm(0) - !X.S(0)) / !X.S(1)
   yPos = (xyzNorm(1) - !Y.S(0)) / !Y.S(1)
   zPos = (xyzNorm(2) - !Z.S(0)) / !Z.S(1)

   ; Return the 3D data coordinate on a face of the cube.
   return, [xPos, yPos, zPos]

end
;******************************************************************************


;******************************************************************************
; Event handler for slice events.
pro Viz3D_SliceEvent, event

   WIDGET_CONTROL, event.id, GET_UVALUE=uVal

   if (uVal eq 'wMainDraw') then begin
      wDrawBase = WIDGET_INFO(event.id, /PARENT)
      WIDGET_CONTROL, wDrawBase, GET_UVALUE=wMainBase
   endif else begin
      wMainBase = event.top
      WIDGET_CONTROL, /HOURGLASS
   endelse

   WIDGET_CONTROL, wMainBase, GET_UVALUE=sMainState, /NO_COPY

   case uVal of
      'wMainDraw': begin ; Big view window.
         if (event.release ge 1) then begin
            WIDGET_CONTROL, sMainState.wMainDraw, DRAW_MOTION_EVENTS=0, $
               /DRAW_BUTTON_EVENTS
            WIDGET_CONTROL, /HOURGLASS
         endif

         if (event.release gt 1) then begin ; Toggle plane direction.
            if (sMainState.sSliceState.planeMode eq 1) then begin
               case sMainState.sSliceState.orthoDir of
                  1: begin ; X, toggle to Y.
                     sMainState.sSliceState.orthoDir = 2
                     WIDGET_CONTROL,sMainState.sSliceState.wSliceYBttn, $
                        /SET_BUTTON
                  end
                  2: begin ; Y, toggle to Z.
                     sMainState.sSliceState.orthoDir = 3
                     WIDGET_CONTROL,sMainState.sSliceState.wSliceZBttn, $
                        /SET_BUTTON
                  end
                  3: begin ; Z, toggle to X.
                     sMainState.sSliceState.orthoDir = 1
                     WIDGET_CONTROL,sMainState.sSliceState.wSliceXBttn, $
                        /SET_BUTTON
                  end
               endcase
               Viz3D_SliceShow, sMainState
            endif
         endif

         if ((event.press eq 1) or $
            ((event.press eq 0) and (event.release eq 0))) then begin

            if (event.press eq 1) then begin
               ; Initiate dynamic plane positioning.
               WIDGET_CONTROL, sMainState.wMainDraw, /DRAW_MOTION_EVENTS, $
                  /DRAW_BUTTON_EVENTS
               sMainState.cleanupView = 1
            endif

            ; Compute a 3D cube surface location from the current
            ; cursor position.
            xyzPos = Viz3D_XYZPos(event.x, event.y, sMainState.hFrontDepth, $
               sMainState.hFrontImage, sMainState.sViewState.invTrans, $
               ON_FACE=onFace)
            ; Normalize the 3D coordinate.
            xyzPos = Temporary(xyzPos) / [sMainState.sViewState.xMax, $
               sMainState.sViewState.yMax, sMainState.sViewState.zMax]

            DEVICE, COPY=[0, 0, sMainState.winX, sMainState.winY, 0, 0, $
                          sMainState.pixWin]

            if (onFace gt 0) then begin
               ; Cursor position is currently somewhere on the cube.

               if (sMainState.sSliceState.planeMode eq 1) then begin
                  ; Drag an orthogonal plane.

                  case sMainState.sSliceState.orthoDir of
                     1: begin ; X
                        Plots, xyzPos(0), [0,1,1,0,0], [0,0,1,1,0], /NORMAL, $
                           /T3D, COLOR=Viz3D_TransColor(sMainState, $
                           sMainState.sSliceState.edgeColor)

                        WIDGET_CONTROL, sMainState.wStatText, $
                           SET_VALUE=('X: ' + $
                           STRING(xyzPos(0)*sMainState.sViewState.xMax))
                     end
                     2: begin ; Y
                        Plots, [0,1,1,0,0], xyzPos(1), [0,0,1,1,0], /NORMAL, $
                           /T3D, COLOR=Viz3D_TransColor(sMainState, $
                           sMainState.sSliceState.edgeColor)

                        WIDGET_CONTROL, sMainState.wStatText, $
                           SET_VALUE=('Y: ' + $
                           STRING(xyzPos(1)*sMainState.sViewState.yMax))
                     end
                     3: begin ; Z
                        Plots, [0,1,1,0,0], [0,0,1,1,0], xyzPos(2), /NORMAL, $
                           /T3D, COLOR=Viz3D_TransColor(sMainState, $
                           sMainState.sSliceState.edgeColor)

                        WIDGET_CONTROL, sMainState.wStatText, $
                           SET_VALUE=('Z: ' + $
                           STRING(xyzPos(2)*sMainState.sViewState.zMax))
                     end
                  endcase
                  EMPTY

               endif else begin ; Drag an oblique plane.

                  edgeColor = Viz3D_TransColor(sMainState, $
                                 sMainState.sSliceState.edgeColor)
                  normColor = Viz3D_TransColor(sMainState, $
                                 sMainState.sSliceState.normColor)

                  ; Move plane normal.

                  if (sMainState.sSliceState.obliqMove eq 0) then begin
                     sMainState.sSliceState.obliqNormal = $
                        xyzPos - sMainState.sSliceState.obliqCenter
                     maxDim = sMainState.sViewState.xMax > $
                              sMainState.sViewState.yMax > $
                              sMainState.sViewState.zMax
                     sMainState.sSliceState.obliqNormal(0) = $
                        sMainState.sSliceState.obliqNormal(0) * $
                        sMainState.sViewState.xMax / maxDim
                     sMainState.sSliceState.obliqNormal(1) = $
                        sMainState.sSliceState.obliqNormal(1) * $
                        sMainState.sViewState.yMax / maxDim
                     sMainState.sSliceState.obliqNormal(2) = $
                        sMainState.sSliceState.obliqNormal(2) * $
                        sMainState.sViewState.zMax / maxDim
                     Viz3D_DrawSliceOblique, $
                        sMainState.sSliceState.obliqCenter, $
                        sMainState.sSliceState.obliqNormal, 0, edgeColor, $
                        sMainState.sViewState.xMax, $
                        sMainState.sViewState.yMax, $
                        sMainState.sViewState.zMax, $
                        sMainState.sViewState.zScale, $
                        NORMCOLOR=normColor, /SKIP_FILL
                     sMainState.sSliceState.obliqNormal = $
                        sMainState.sSliceState.obliqNormal / $
                        SQRT(TOTAL(sMainState.sSliceState.obliqNormal^2))
                     WIDGET_CONTROL, sMainState.wStatText, $
                        SET_VALUE=('Normal: ' + $
                        STRING(sMainState.sSliceState.obliqNormal, $
                           FORMAT='(F6.3,1X,F6.3,1X,F6.3)'))

                  endif else begin ; Move plane center.

                     case onFace of
                        1: begin
                           sMainState.sSliceState.obliqCenter([0,1]) = $
                              xyzPos([0,1])
                        end
                        2: begin
                           sMainState.sSliceState.obliqCenter([0,2]) = $
                              xyzPos([0,2])
                        end
                        3: begin
                           sMainState.sSliceState.obliqCenter([1,2]) = $
                              xyzPos([1,2])
                        end
                     endcase
                     Viz3D_DrawSliceOblique, $
                        sMainState.sSliceState.obliqCenter, $
                        sMainState.sSliceState.obliqNormal, 0, edgeColor, $
                        sMainState.sViewState.xMax, $
                        sMainState.sViewState.yMax, $
                        sMainState.sViewState.zMax, $
                        sMainState.sViewState.zScale, $
                        NORMCOLOR=normColor, /SKIP_FILL
                     xyzPos = sMainState.sSliceState.obliqCenter * $
                              [sMainState.sViewState.xMax, $
                               sMainState.sViewState.yMax, $
                               sMainState.sViewState.zMax]
                     WIDGET_CONTROL, sMainState.wStatText, $
                        SET_VALUE=('Center: ' + $
                        STRING(xyzPos, FORMAT='(F6.3,1X,F6.3,1X,F6.3)'))
                  endelse
               endelse
            endif
         endif

         if (event.release eq 1) then begin ; Cleanup and draw slice.
            xyzPos = Viz3D_XYZPos(event.x, event.y, sMainState.hFrontDepth, $
               sMainState.hFrontImage, sMainState.sViewState.invTrans, $
               ON_FACE=onFace)
            xyzPos = Temporary(xyzPos) / [sMainState.sViewState.xMax, $
               sMainState.sViewState.yMax, sMainState.sViewState.zMax]

            if (TOTAL(ABS(xyzPos)) gt 0.0) then begin
               ; Plane is now positioned.

               if (sMainState.sSliceState.planeMode eq 1) then begin
                  ; Orthogonal plane.
                  case sMainState.sSliceState.orthoDir of
                     1: sMainState.sSliceState.orthoPos = xyzPos(0)
                     2: sMainState.sSliceState.orthoPos = xyzPos(1)
                     3: sMainState.sSliceState.orthoPos = xyzPos(2)
                  endcase
                  ; Update the small slice status window.
                  Viz3D_SliceShow, sMainState

                  ; Draw orthogonal slice with data in the big window.
                  Viz3D_OrthoPlaneDraw, sMainState
               endif else begin ; Oblique plane.
                  ; Normal moved.
                  if (sMainState.sSliceState.obliqMove eq 0) then begin
                     sMainState.sSliceState.obliqNormal = $
                        xyzPos - sMainState.sSliceState.obliqCenter
                     maxDim = sMainState.sViewState.xMax > $
                              sMainState.sViewState.yMax > $
                              sMainState.sViewState.zMax
                     sMainState.sSliceState.obliqNormal(0) = $
                        sMainState.sSliceState.obliqNormal(0) * $
                        sMainState.sViewState.xMax / maxDim
                     sMainState.sSliceState.obliqNormal(1) = $
                        sMainState.sSliceState.obliqNormal(1) * $
                        sMainState.sViewState.yMax / maxDim
                     sMainState.sSliceState.obliqNormal(2) = $
                        sMainState.sSliceState.obliqNormal(2) * $
                        sMainState.sViewState.zMax / maxDim
                  endif else begin ; Center moved.
                     case onFace of
                        1: begin
                           sMainState.sSliceState.obliqCenter([0,1]) = $
                              xyzPos([0,1])
                        end
                        2: begin
                           sMainState.sSliceState.obliqCenter([0,2]) = $
                              xyzPos([0,2])
                        end
                        3: begin
                           sMainState.sSliceState.obliqCenter([1,2]) = $
                              xyzPos([1,2])
                        end
                     endcase
                  endelse
                  fillColor = Viz3D_TransColor(sMainState, $
                                 sMainState.sSliceState.fillColor)
                  edgeColor = Viz3D_TransColor(sMainState, $
                                 sMainState.sSliceState.edgeColor)
                  normColor = Viz3D_TransColor(sMainState, $
                                 sMainState.sSliceState.normColor)
                  backColor = Viz3D_TransColor(sMainState, $
                                 sMainState.backColor)

                  ; Update the small slice status window.
                  Viz3D_SliceShow, sMainState
               endelse
            endif
         endif
      end
      'wSliceDraw': begin ; Event in small slice plane window.
         if (event.release gt 0) then begin
            ; Toggle the orthogonal plane direction.

            if (sMainState.sSliceState.planeMode eq 1) then begin
               case sMainState.sSliceState.orthoDir of
                  1: begin
                     sMainState.sSliceState.orthoDir = 2
                     WIDGET_CONTROL,sMainState.sSliceState.wSliceYBttn, $
                        /SET_BUTTON
                  end
                  2: begin
                     sMainState.sSliceState.orthoDir = 3
                     WIDGET_CONTROL,sMainState.sSliceState.wSliceZBttn, $
                        /SET_BUTTON
                  end
                  3: begin
                     sMainState.sSliceState.orthoDir = 1
                     WIDGET_CONTROL,sMainState.sSliceState.wSliceXBttn, $
                        /SET_BUTTON
                  end
               endcase
            endif else begin ; Reset the oblique slicing plane position.
               sMainState.sSliceState.obliqCenter = [0.5, 0.5, 0.5]
               sMainState.sSliceState.obliqNormal = [0.0, 0.0, 1.0]
            endelse

            ; Update the small slice status window.
            Viz3D_SliceShow, sMainState
            DEVICE, COPY=[0, 0, sMainState.winX, sMainState.winY, $
                          0, 0, sMainState.pixWin]

            WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=''
         endif
      end
      'wSliceOrthoBttn': begin ; Set orthogonal slice mode.
         sMainState.sSliceState.planeMode = 1
         WIDGET_CONTROL, sMainState.sSliceState.wSliceObliqBase, MAP=0
         WIDGET_CONTROL, sMainState.sSliceState.wSliceOrthoBase, MAP=1
         Viz3D_SliceShow, sMainState
         DEVICE, COPY=[0, 0, sMainState.winX, sMainState.winY, $
                       0, 0, sMainState.pixWin]
         WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=''
      end
      'wSliceObliqBttn': begin ; Set oblique slice mode.
         sMainState.sSliceState.planeMode = 0
         WIDGET_CONTROL, sMainState.sSliceState.wSliceOrthoBase, MAP=0
         WIDGET_CONTROL, sMainState.sSliceState.wSliceObliqBase, MAP=1
         fillColor = Viz3D_TransColor(sMainState, $
                        sMainState.sSliceState.fillColor)
         edgeColor = Viz3D_TransColor(sMainState, $
                        sMainState.sSliceState.edgeColor)
         normColor = Viz3D_TransColor(sMainState, $
                        sMainState.sSliceState.normColor)
         backColor = Viz3D_TransColor(sMainState, $
                        sMainState.backColor)
         Viz3D_SliceShow, sMainState
      end
      'wSliceDrawBttn': begin ; Set "Draw" mode for slices.
         sMainState.sSliceState.drawMode = 0
      end
      'wSliceExpoBttn': begin ; Set "Expose" mode for slices.
         sMainState.sSliceState.drawMode = 1
      end
      'wSliceXBttn': begin ; Set orthogonal slice direction X.
         sMainState.sSliceState.orthoDir = 1
         Viz3D_SliceShow, sMainState
      end
      'wSliceYBttn': begin ; Set orthogonal slice direction Y.
         sMainState.sSliceState.orthoDir = 2
         Viz3D_SliceShow, sMainState
      end
      'wSliceZBttn': begin ; Set orthogonal slice direction Z.
         sMainState.sSliceState.orthoDir = 3
         Viz3D_SliceShow, sMainState
      end
      'wSliceNormalBttn': begin ; Set oblique slice mode for normal movement.
         sMainState.sSliceState.obliqMove = 0
      end
      'wSliceCenterBttn': begin ; Set oblique slice mode for center movement.
         sMainState.sSliceState.obliqMove = 1
      end
      'wSliceObliqGoBttn': begin ; Draw the slice.
         fillColor = Viz3D_TransColor(sMainState, $
                        sMainState.sSliceState.fillColor)
         edgeColor = Viz3D_TransColor(sMainState, $
                        sMainState.sSliceState.edgeColor)
         Viz3D_DrawSliceOblique, sMainState.sSliceState.obliqCenter, $
            sMainState.sSliceState.obliqNormal, $
            fillColor, edgeColor, $
            sMainState.sViewState.xMax, $
            sMainState.sViewState.yMax, $
            sMainState.sViewState.zMax, $
            sMainState.sViewState.zScale, /SKIP_FILL, $
            VERT_PLANE=vertPlane, VERT_3D=vert3D
         ; Draw oblique slice with data in the big window.
         Viz3D_ObliqPlaneDraw, sMainState, vertPlane, vert3D
      end
      else:
   endcase

   WIDGET_CONTROL, wMainBase, SET_UVALUE=sMainState, /NO_COPY

end
;******************************************************************************


;******************************************************************************
; Event handler for block events.
pro Viz3D_BlockEvent, event

   WIDGET_CONTROL, event.id, GET_UVALUE=uVal
   if (uVal eq 'wMainDraw') then begin
      wDrawBase = WIDGET_INFO(event.id, /PARENT)
      WIDGET_CONTROL, wDrawBase, GET_UVALUE=wMainBase
   endif else begin
      wMainBase = event.top
      WIDGET_CONTROL, /HOURGLASS
   endelse

   WIDGET_CONTROL, wMainBase, GET_UVALUE=sMainState, /NO_COPY

   case uVal of
      'wMainDraw': begin ; Big draw window, move block.
         if (event.release ge 1) then begin
            WIDGET_CONTROL, sMainState.wMainDraw, DRAW_MOTION_EVENTS=0, $
               /DRAW_BUTTON_EVENTS
            WIDGET_CONTROL, /HOURGLASS
         endif

         if ((event.press ge 1) or $
            ((event.press eq 0) and (event.release eq 0))) then begin

            if (event.press ge 1) then begin
               ; Initiate dynamic block positioning.
               WIDGET_CONTROL, sMainState.wMainDraw, /DRAW_MOTION_EVENTS, $
                  /DRAW_BUTTON_EVENTS
               sMainState.cleanupView = 1
            endif
            if (event.press eq 1) then $
               sMainState.sBlockState.frontBack = 1
            if (event.press gt 1) then $
               sMainState.sBlockState.frontBack = 0

            ; Compute a 3D cube surface location from the current
            ; cursor position.
            xyzPos = Viz3D_XYZPos(event.x, event.y, sMainState.hFrontDepth, $
               sMainState.hFrontImage, sMainState.sViewState.invTrans, $
               ON_FACE=onFace)

            xyzPos = ROUND(xyzPos)

            if (onFace gt 0) then begin ; On cube, draw block.
               case (onFace) of
                  1: begin
                     if (sMainState.sBlockState.frontBack) then begin
                        sMainState.sBlockState.c1(0) = xyzPos(0)
                        sMainState.sBlockState.c1(1) = xyzPos(1)
                     endif else begin
                        sMainState.sBlockState.c2(0) = xyzPos(0)
                        sMainState.sBlockState.c2(1) = xyzPos(1)
                     endelse
                  end
                  2: begin
                     if (sMainState.sBlockState.frontBack) then begin
                        sMainState.sBlockState.c1(0) = xyzPos(0)
                        sMainState.sBlockState.c1(2) = xyzPos(2)
                     endif else begin
                        sMainState.sBlockState.c2(0) = xyzPos(0)
                        sMainState.sBlockState.c2(2) = xyzPos(2)
                     endelse
                  end
                  3: begin
                     if (sMainState.sBlockState.frontBack) then begin
                        sMainState.sBlockState.c1(1) = xyzPos(1)
                        sMainState.sBlockState.c1(2) = xyzPos(2)
                     endif else begin
                        sMainState.sBlockState.c2(1) = xyzPos(1)
                        sMainState.sBlockState.c2(2) = xyzPos(2)
                     endelse
                  end
               endcase
            endif
         endif

         DEVICE, COPY=[0, 0, sMainState.winX, sMainState.winY, 0, 0, $
                       sMainState.pixWin]
         Viz3D_BlockShow, sMainState, /DIRECT

         if (event.release ge 1) then begin ; Cleanup.
            Viz3D_BlockShow, sMainState
         endif

         WIDGET_CONTROL, sMainState.wStatText, $
            SET_VALUE=('(' + $
            STRTRIM(STRING(sMainState.sBlockState.c1(0)),2) + ', ' + $
            STRTRIM(STRING(sMainState.sBlockState.c1(1)),2) + ', ' + $
            STRTRIM(STRING(sMainState.sBlockState.c1(2)),2) + '), (' + $
            STRTRIM(STRING(sMainState.sBlockState.c2(0)),2) + ', ' + $
            STRTRIM(STRING(sMainState.sBlockState.c2(1)),2) + ', ' + $
            STRTRIM(STRING(sMainState.sBlockState.c2(2)),2) + ')')

      end
      'wBlockDraw': ; Small block window, do nothing.
      'wBlockAddBttn': sMainState.sBlockState.blockMode = 1
      'wBlockSubBttn': sMainState.sBlockState.blockMode = 0
      'wBlockGoBttn': begin ; Draw it.
         Viz3D_BlockDraw, sMainState
      end
   endcase

   WIDGET_CONTROL, wMainBase, SET_UVALUE=sMainState, /NO_COPY
end
;******************************************************************************


;******************************************************************************
; Event handler for surface events.
pro Viz3D_SurfEvent, event

   WIDGET_CONTROL, event.id, GET_UVALUE=uVal
   if (uVal ne 'wSurfDraw') then WIDGET_CONTROL, /HOURGLASS

   if (uVal eq 'wMainDraw') then begin
      wDrawBase = WIDGET_INFO(event.id, /PARENT)
      WIDGET_CONTROL, wDrawBase, GET_UVALUE=wMainBase
   endif else wMainBase = event.top

   WIDGET_CONTROL, wMainBase, GET_UVALUE=sMainState, /NO_COPY

   case uVal of
      'wMainDraw': ; Big draw window, do nothing.
      'wSurfDraw': begin
         if (event.release le 0) then begin
            ; Determine which to set, high or low.
            WSET, sMainState.sThreshState.threshWin
            norm = FLOAT(event.x) / FLOAT(!D.X_SIZE)
            norm = (norm - 0.03) / (0.97 - 0.03)
            lTh = (*(sMainState.sThreshState.hLowThresh))(sMainState.curData)
            hTh = (*(sMainState.sThreshState.hHighThresh))(sMainState.curData)
            newThresh = (norm * (hTh - lTh)) + lTh
            surfThresh = *(sMainState.sSurfState.hSurfThresh)
            surfThresh(sMainState.curData) = newThresh
            *(sMainState.sSurfState.hSurfThresh) = surfThresh
            WIDGET_CONTROL, event.id, /DRAW_MOTION_EVENTS, $
               /DRAW_BUTTON_EVENTS
            WSET, sMainState.mainWin
         endif else begin
            WIDGET_CONTROL, event.id, DRAW_MOTION_EVENTS=0, $
               /DRAW_BUTTON_EVENTS
            WIDGET_CONTROL, /HOURGLASS
         endelse

         Viz3D_SurfShow, sMainState

      end
      'wSurfText': begin
         WIDGET_CONTROL, event.id, GET_VALUE=newThresh
         lTh = (*(sMainState.sThreshState.hLowThresh))(sMainState.curData)
         hTh = (*(sMainState.sThreshState.hHighThresh))(sMainState.curData)
         newThresh = FLOAT(newThresh(0))
         newThresh = (newThresh > lTh) < hTh
         surfThresh = *(sMainState.sSurfState.hSurfThresh)
         surfThresh(sMainState.curData) = newThresh
         *(sMainState.sSurfState.hSurfThresh) = surfThresh
         WIDGET_CONTROL, event.id, SET_VALUE=STRTRIM(STRING(newThresh),2)
         Viz3D_SurfShow, sMainState
      end
      'wSurfLowBttn': sMainState.sSurfState.surfSide = 0
      'wSurfHighBttn': sMainState.sSurfState.surfSide = 1
      'wSurfShadeDrop': sMainState.sSurfState.curShade = event.index
      'wSurfGoBttn': Viz3D_SurfDraw, sMainState, /CALC_SURF
   endcase

   WIDGET_CONTROL, wMainBase, SET_UVALUE=sMainState, /NO_COPY
end
;******************************************************************************


;******************************************************************************
; Event handler for projection events.
pro Viz3D_ProjectEvent, event
   WIDGET_CONTROL, /HOURGLASS

   WIDGET_CONTROL, event.id, GET_UVALUE=uVal
   if (uVal eq 'wMainDraw') then begin
      wDrawBase = WIDGET_INFO(event.id, /PARENT)
      WIDGET_CONTROL, wDrawBase, GET_UVALUE=wMainBase
   endif else wMainBase = event.top

   WIDGET_CONTROL, wMainBase, GET_UVALUE=sMainState, /NO_COPY

   case uVal of
      'wMainDraw': ; Big draw window, do nothing.
      'wProjMaxBttn': sMainState.sProjState.projType = 0
      'wProjAvgBttn': sMainState.sProjState.projType = 1
      'wProjLowBttn': sMainState.sProjState.projReso = 0
      'wProjMedBttn': sMainState.sProjState.projReso = 1
      'wProjHighBttn': sMainState.sProjState.projReso = 2
      'wProjDQSlid': sMainState.sProjState.depthQ = event.value
      'wProjGoBttn': Viz3D_ProjDraw, sMainState
   endcase

   WIDGET_CONTROL, wMainBase, SET_UVALUE=sMainState, /NO_COPY
end
;******************************************************************************


;******************************************************************************
; Event handler for threshold events.
pro Viz3D_ThreshEvent, event

   WIDGET_CONTROL, event.id, GET_UVALUE=uVal
   if (uVal ne 'wThreshDraw') then WIDGET_CONTROL, /HOURGLASS

   if (uVal eq 'wMainDraw') then begin
      wDrawBase = WIDGET_INFO(event.id, /PARENT)
      WIDGET_CONTROL, wDrawBase, GET_UVALUE=wMainBase
   endif else wMainBase = event.top

   WIDGET_CONTROL, wMainBase, GET_UVALUE=sMainState, /NO_COPY

   lTh = (*(sMainState.sThreshState.hLowThresh))(sMainState.curData)
   hTh = (*(sMainState.sThreshState.hHighThresh))(sMainState.curData)
   vTp = (*(sMainState.sThreshState.hTransVal))(sMainState.curData)

   minD = (*(sMainState.sThreshState.hMinData))(sMainState.curData)
   maxD = (*(sMainState.sThreshState.hMaxData))(sMainState.curData)

   case uVal of
      'wMainDraw': ; Big draw window, do nothing.
      'wThreshDraw': begin ; Set new threshold value via draw widget.
         if (event.release le 0) then begin
            if (event.press eq 1) then sMainState.sThreshState.moveType = 0
            if (event.press ge 2) then sMainState.sThreshState.moveType = 1
            if (sMainState.sThreshState.moveType eq 0) then begin
               ; Move high/low threshold
               ; Determine which to set, high or low.
               WSET, sMainState.sThreshState.threshWin
               norm = FLOAT(event.x) / FLOAT(!D.X_SIZE)
               norm = (norm - 0.03) / (0.97 - 0.03)
               newThresh = (norm * (maxD - minD)) + minD
               newThresh = (newThresh > minD) < maxD

               ; Closer to min threshold.
               if (ABS(newThresh - lTh) lt ABS(newThresh - hTh)) then begin
                  newThresh = newThresh < hTh
                  lowThresh = *(sMainState.sThreshState.hLowThresh)
                  lowThresh(sMainState.curData) = newThresh
                  *(sMainState.sThreshState.hLowThresh) = lowThresh
               endif 
               ; Closer to max threshold.
               if (ABS(newThresh - lTh) gt ABS(newThresh - hTh)) then begin
                  newThresh = newThresh > lTh
                  highThresh = *(sMainState.sThreshState.hHighThresh)
                  highThresh(sMainState.curData) = newThresh
                  *(sMainState.sThreshState.hHighThresh) = highThresh
               endif
               if (ABS(newThresh - lTh) eq ABS(newThresh - hTh)) then begin
                  newThresh = (newThresh > minD) < maxD
                  if (newThresh gt hTh) then begin
                     highThresh = *(sMainState.sThreshState.hHighThresh)
                     highThresh(sMainState.curData) = newThresh
                     *(sMainState.sThreshState.hHighThresh) = highThresh
                  endif
                  if (newThresh lt lTh) then begin
                     lowThresh = *(sMainState.sThreshState.hLowThresh)
                     lowThresh(sMainState.curData) = newThresh
                     *(sMainState.sThreshState.hLowThresh) = lowThresh
                  endif
               endif

               ; Make sure the min and max threshold values are not equal.
               lowThresh = *(sMainState.sThreshState.hLowThresh)
               highThresh = *(sMainState.sThreshState.hHighThresh)
               if (lowThresh(sMainState.curData) eq $
                   highThresh(sMainState.curData)) then begin
                  lowThresh(sMainState.curData) = minD
                  highThresh(sMainState.curData) = maxD
               endif
               *(sMainState.sThreshState.hLowThresh) = lowThresh
               *(sMainState.sThreshState.hHighThresh) = highThresh

            endif else begin ; Move transparency.
               WSET, sMainState.sThreshState.threshWin
               norm = FLOAT(event.x) / FLOAT(!D.X_SIZE)
               norm = (norm - 0.03) / (0.97 - 0.03)
               newThresh = (norm * (maxD - minD)) + minD
               newThresh = (newThresh > minD) < maxD

               transVal = *(sMainState.sThreshState.hTransVal)
               transVal(sMainState.curData) = newThresh
               *(sMainState.sThreshState.hTransVal) = transVal
            endelse

            WIDGET_CONTROL, event.id, /DRAW_MOTION_EVENTS, $
               /DRAW_BUTTON_EVENTS
            WSET, sMainState.mainWin
         endif else begin
            WIDGET_CONTROL, event.id, DRAW_MOTION_EVENTS=0, $
               /DRAW_BUTTON_EVENTS
            WIDGET_CONTROL, /HOURGLASS
            if (event.release eq 1) then begin
               WIDGET_CONTROL, /HOURGLASS
               Viz3D_SurfShow, sMainState, /CALC_HIST
            endif
         endelse

         Viz3D_ThreshShow, sMainState, /DYNAMIC
      end
      'wThreshLowText': begin
         WIDGET_CONTROL, event.id, GET_VALUE=newThresh
         newThresh = FLOAT(newThresh(0))
         newThresh = (newThresh > minD) < hTh

         lowThresh = *(sMainState.sThreshState.hLowThresh)
         highThresh = *(sMainState.sThreshState.hHighThresh)
         lowThresh(sMainState.curData) = newThresh
         if (lowThresh(sMainState.curData) eq $
             highThresh(sMainState.curData)) then begin
            lowThresh(sMainState.curData) = minD
            highThresh(sMainState.curData) = maxD
         endif
         *(sMainState.sThreshState.hLowThresh) = lowThresh
         *(sMainState.sThreshState.hHighThresh) = highThresh

         WIDGET_CONTROL, event.id, SET_VALUE=STRTRIM(STRING(newThresh),2)
         Viz3D_ThreshShow, sMainState, /DYNAMIC
         Viz3D_SurfShow, sMainState, /CALC_HIST
      end
      'wThreshHighText': begin
         WIDGET_CONTROL, event.id, GET_VALUE=newThresh
         newThresh = FLOAT(newThresh(0))
         newThresh = (newThresh > lTh) < maxD

         lowThresh = *(sMainState.sThreshState.hLowThresh)
         highThresh = *(sMainState.sThreshState.hHighThresh)
         highThresh(sMainState.curData) = newThresh
         if (lowThresh(sMainState.curData) eq $
             highThresh(sMainState.curData)) then begin
            lowThresh(sMainState.curData) = minD
            highThresh(sMainState.curData) = maxD
         endif
         *(sMainState.sThreshState.hHighThresh) = highThresh
         *(sMainState.sThreshState.hLowThresh) = lowThresh

         WIDGET_CONTROL, event.id, SET_VALUE=STRTRIM(STRING(newThresh),2)
         Viz3D_ThreshShow, sMainState, /DYNAMIC
         Viz3D_SurfShow, sMainState, /CALC_HIST
      end
      'wThreshTransText': begin
         WIDGET_CONTROL, event.id, GET_VALUE=newThresh
         newThresh = FLOAT(newThresh(0))
         newThresh = (newThresh > minD) < maxD

         transVal = *(sMainState.sThreshState.hTransVal)
         transVal(sMainState.curData) = newThresh
         *(sMainState.sThreshState.hTransVal) = transVal

         WIDGET_CONTROL, event.id, SET_VALUE=STRTRIM(STRING(newThresh),2)
         Viz3D_ThreshShow, sMainState, /DYNAMIC
      end
   endcase

   WIDGET_CONTROL, wMainBase, SET_UVALUE=sMainState, /NO_COPY
end
;******************************************************************************


;******************************************************************************
; Event handler for profile events.
pro Viz3D_ProfileEvent, event

   WIDGET_CONTROL, event.id, GET_UVALUE=uVal
   if (uVal eq 'wMainDraw') then begin
      wDrawBase = WIDGET_INFO(event.id, /PARENT)
      WIDGET_CONTROL, wDrawBase, GET_UVALUE=wMainBase
   endif else begin
      wMainBase = event.top
      WIDGET_CONTROL, /HOURGLASS
   endelse

   WIDGET_CONTROL, wMainBase, GET_UVALUE=sMainState, /NO_COPY

   case uVal of
      'wProfDraw': ; Profile window, do nothing.
      'wMainDraw': begin ; Big draw window.
         if (event.release ge 1) then begin
            WIDGET_CONTROL, sMainState.wMainDraw, DRAW_MOTION_EVENTS=0, $
               /DRAW_BUTTON_EVENTS
            WIDGET_CONTROL, /HOURGLASS
         endif

         if ((event.press ge 1) or $
            ((event.press eq 0) and (event.release eq 0))) then begin

            if (event.press ge 1) then begin
               ; Initiate dynamic profile positioning.
               WIDGET_CONTROL, sMainState.wMainDraw, /DRAW_MOTION_EVENTS, $
                  /DRAW_BUTTON_EVENTS
               sMainState.cleanupView = 1
            endif
            if (event.press eq 1) then $
               sMainState.sProfState.frontBack = 1
            if (event.press gt 1) then $
               sMainState.sProfState.frontBack = 0

            ; Compute a 3D cube surface location from the current
            ; cursor position.
            if (sMainState.sProfState.frontBack) then $
               xyzPos = Viz3D_XYZPos(event.x, event.y, $
                  sMainState.hFrontDepth, sMainState.hFrontImage, $
                  sMainState.sViewState.invTrans, ON_FACE=onFace) $
            else $
               xyzPos = Viz3D_XYZPos(event.x, event.y, $
                  sMainState.hBackDepth, sMainState.hBackImage, $
                  sMainState.sViewState.invTrans, ON_FACE=onFace)

            xyzPos = ROUND(xyzPOS)

            if (onFace gt 0) then begin ; On cube, draw profile.
               if (sMainState.sProfState.profType eq 0) then begin
                  ; Orthogonal profile.
                  sMainState.sProfState.profDir = onFace
                  case onFace of
                     1: begin
                        sMainState.sProfState.x1Ortho = xyzPos(0)
                        sMainState.sProfState.x2Ortho = xyzPos(0)
                        sMainState.sProfState.y1Ortho = xyzPos(1)
                        sMainState.sProfState.y2Ortho = xyzPos(1)
                        sMainState.sProfState.z1Ortho = 0
                        sMainState.sProfState.z2Ortho = $
                           sMainState.szData(3) - 1
                     end
                     2: begin
                        sMainState.sProfState.x1Ortho = xyzPos(0)
                        sMainState.sProfState.x2Ortho = xyzPos(0)
                        sMainState.sProfState.y1Ortho = 0
                        sMainState.sProfState.y2Ortho = $
                           sMainState.szData(2) - 1
                        sMainState.sProfState.z1Ortho = xyzPos(2)
                        sMainState.sProfState.z2Ortho = xyzPos(2)
                     end
                     3: begin
                        sMainState.sProfState.x1Ortho = 0
                        sMainState.sProfState.x2Ortho = $
                           sMainState.szData(1) - 1
                        sMainState.sProfState.y1Ortho = xyzPos(1)
                        sMainState.sProfState.y2Ortho = xyzPos(1)
                        sMainState.sProfState.z1Ortho = xyzPos(2)
                        sMainState.sProfState.z2Ortho = xyzPos(2)
                     end
                  endcase
               endif else begin ; Oblique profile.
                  if (sMainState.sProfState.frontBack) then begin
                     sMainState.sProfState.x1Obliq = xyzPos(0)
                     sMainState.sProfState.y1Obliq = xyzPos(1)
                     sMainState.sProfState.z1Obliq = xyzPos(2)
                  endif else begin
                     sMainState.sProfState.x2Obliq = xyzPos(0)
                     sMainState.sProfState.y2Obliq = xyzPos(1)
                     sMainState.sProfState.z2Obliq = xyzPos(2)
                  endelse
               endelse

               if (onFace ne sMainState.sProfState.onFace) then begin
                  sMainState.sProfState.onFace = onFace
                  Viz3D_ProfShow, sMainState
               endif else Viz3D_ProfShow, sMainState, /DIRECT

               Viz3D_ProfDraw, sMainState

            endif
         endif
      end
      'wProfOrthoBttn': begin
         sMainState.sProfState.profType = 0
         Viz3D_ProfShow, sMainState
         Viz3D_ProfDraw, sMainState
      end
      'wProfObliqBttn': begin
         sMainState.sProfState.profType = 1
         Viz3D_ProfShow, sMainState
         Viz3D_ProfDraw, sMainState
      end
   endcase

   WIDGET_CONTROL, wMainBase, SET_UVALUE=sMainState, /NO_COPY
end
;******************************************************************************


;******************************************************************************
; Event handler for probe events.
pro Viz3D_ProbeEvent, event

   WIDGET_CONTROL, event.id, GET_UVALUE=uVal
   if (uVal eq 'wMainDraw') then begin
      wDrawBase = WIDGET_INFO(event.id, /PARENT)
      WIDGET_CONTROL, wDrawBase, GET_UVALUE=wMainBase
   endif else begin
      wMainBase = event.top
      WIDGET_CONTROL, /HOURGLASS
   endelse

   WIDGET_CONTROL, wMainBase, GET_UVALUE=sMainState, /NO_COPY

   case uVal of
      'wMainDraw': begin ; Big draw window.
         if (event.release ge 1) then begin
            WIDGET_CONTROL, sMainState.wMainDraw, DRAW_MOTION_EVENTS=0, $
               /DRAW_BUTTON_EVENTS
            WIDGET_CONTROL, /HOURGLASS
         endif

         if (event.press ge 1) then $
            WIDGET_CONTROL, sMainState.wMainDraw, /DRAW_MOTION_EVENTS, $
               /DRAW_BUTTON_EVENTS

         if ((event.press ge 1) or $
            ((event.press eq 0) and (event.release eq 0))) then begin

            eX = (event.x > 0) < (!D.X_Size - 1)
            eY = (event.y > 0) < (!D.Y_Size - 1)
            nX = Float(eX) / Float(!D.X_Size - 1)
            nY = Float(eY) / Float(!D.Y_Size - 1)
            SET_PLOT, 'Z'
            nZ = (TVRD(eX, eY, 1, 1, CHANNEL=1, /WORDS))(0)
            SET_PLOT, sMainState.screenDevice

            zBack = ((*(sMainState.hBackDepth))(eX, eY))(0)

            zFront = ((*(sMainState.hfrontDepth))(eX, eY))(0)

            if (zFront gt (-32765)) then begin

               if (nZ gt zBack) then begin
                  nZ = (FLOAT(nZ) + 32765.0) / (2.0 * 32765.0)

                  xyzNorm = [nX, nY, nZ, 1.0] # sMainState.sViewState.invTrans
                  xyzNorm = xyzNorm / xyzNorm(3)

                  xPos = (((xyzNorm(0) - !X.S(0)) / !X.S(1)) > 0.0) < $
                     FLOAT(sMainState.szData(1)-1)
                  yPos = (((xyzNorm(1) - !Y.S(0)) / !Y.S(1)) > 0.0) < $
                     FLOAT(sMainState.szData(2)-1)
                  zPos = (((xyzNorm(2) - !Z.S(0)) / !Z.S(1)) > 0.0) < $
                     FLOAT(sMainState.szData(3)-1)

                  sMainState.sProbeState.x = xPos
                  sMainState.sProbeState.y = yPos
                  sMainState.sProbeState.z = zPos

                  Viz3D_ProbeDraw, sMainState

               endif else begin
                  WIDGET_CONTROL, sMainState.wStatText, $
                     SET_VALUE=('No data value')
                  WSET, sMainState.mainWin
                  DEVICE, COPY=[0, 0, sMainState.winX, sMainState.winY, $
                                0, 0, sMainState.pixWin]
                  WIDGET_CONTROL, sMainState.sProbeState.wProbeXText, $
                     SET_VALUE=''
                  WIDGET_CONTROL, sMainState.sProbeState.wProbeYText, $
                     SET_VALUE=''
                  WIDGET_CONTROL, sMainState.sProbeState.wProbeZText, $
                     SET_VALUE=''
               endelse
            endif else begin
               WIDGET_CONTROL, sMainState.wStatText, $
                  SET_VALUE=('')
               WIDGET_CONTROL, sMainState.sProbeState.wProbeXText, $
                  SET_VALUE=''
               WIDGET_CONTROL, sMainState.sProbeState.wProbeYText, $
                  SET_VALUE=''
               WIDGET_CONTROL, sMainState.sProbeState.wProbeZText, $
                  SET_VALUE=''
               WSET, sMainState.mainWin
               DEVICE, COPY=[0, 0, sMainState.winX, sMainState.winY, $
                             0, 0, sMainState.pixWin]
            endelse
         endif
      end
      'wProbeText': begin
         WIDGET_CONTROL, sMainState.sProbeState.wProbeXText, GET_VALUE=xT
         WIDGET_CONTROL, sMainState.sProbeState.wProbeYText, GET_VALUE=yT
         WIDGET_CONTROL, sMainState.sProbeState.wProbeZText, GET_VALUE=zT
         sMainState.sProbeState.x = (FLOAT(xt(0)) > 0.0) < $
            FLOAT(sMainState.szData(1)-1)
         sMainState.sProbeState.y = (FLOAT(yt(0)) > 0.0) < $
            FLOAT(sMainState.szData(2)-1)
         sMainState.sProbeState.z = (FLOAT(zt(0)) > 0.0) < $
            FLOAT(sMainState.szData(3)-1)
         xT = STRTRIM(STRING(sMainState.sProbeState.x), 2)
         yT = STRTRIM(STRING(sMainState.sProbeState.y), 2)
         zT = STRTRIM(STRING(sMainState.sProbeState.z), 2)
         WIDGET_CONTROL, sMainState.sProbeState.wProbeXText, SET_VALUE=xT
         WIDGET_CONTROL, sMainState.sProbeState.wProbeYText, SET_VALUE=yT
         WIDGET_CONTROL, sMainState.sProbeState.wProbeZText, SET_VALUE=zT
         Viz3D_ProbeDraw, sMainState
      end
   endcase

   WIDGET_CONTROL, wMainBase, SET_UVALUE=sMainState, /NO_COPY
end
;******************************************************************************


;******************************************************************************
; Event handler for view events.
pro Viz3D_ViewEvent, event
   WIDGET_CONTROL, /HOURGLASS

   WIDGET_CONTROL, event.id, GET_UVALUE=uVal
   if (uVal eq 'wMainDraw') then begin
      wDrawBase = WIDGET_INFO(event.id, /PARENT)
      WIDGET_CONTROL, wDrawBase, GET_UVALUE=wMainBase
   endif else wMainBase = event.top

   WIDGET_CONTROL, wMainBase, GET_UVALUE=sMainState, /NO_COPY

   Viz3D_CleanView, sMainState, /SKIP_BUFFER

   redraw_flag = 0
   case uVal of
      'wMainDraw': ; Big draw window, do nothing.
      'wViewDraw': ; Small view window, do nothing.
      'wRot1Slid': begin
         if (sMainState.sViewState.ang1 ne event.value) then redraw_flag = 1
         sMainState.sViewState.ang1 = event.value
      end
      'wRot2Slid': begin
         if (sMainState.sViewState.ang2 ne event.value) then redraw_flag = 1
         sMainState.sViewState.ang2 = event.value
      end
      'wRot3Slid': begin
         if (sMainState.sViewState.ang3 ne event.value) then redraw_flag = 1
         sMainState.sViewState.ang3 = event.value
      end
      'wRot1Drop': begin
         case event.index of
            0: rot_axis = 'X'
            1: rot_axis = 'Y'
            2: rot_axis = 'Z'
         endcase
         if (sMainState.sViewState.dir1 ne rot_axis) then redraw_flag = 1
         sMainState.sViewState.dir1 = rot_axis
      end
      'wRot2Drop': begin
         case event.index of
            0: rot_axis = 'X'
            1: rot_axis = 'Y'
            2: rot_axis = 'Z'
         endcase
         if (sMainState.sViewState.dir2 ne rot_axis) then redraw_flag = 1
         sMainState.sViewState.dir2 = rot_axis
      end
      'wRot3Drop': begin
         case event.index of
            0: rot_axis = 'X'
            1: rot_axis = 'Y'
            2: rot_axis = 'Z'
         endcase
         if (sMainState.sViewState.dir3 ne rot_axis) then redraw_flag = 1
         sMainState.sViewState.dir3 = rot_axis
      end
      'wZoomSlid': begin
         zoomVal = FLOAT(event.value) / 100.0
         if (sMainState.sViewState.zoomFac ne zoomVal) then redraw_flag = 1
         sMainState.sViewState.zoomFac = zoomVal
      end
      'wZFacSlid': begin
         zFac = FLOAT(event.value) / 100.0
         if (sMainState.sViewState.zScale ne zFac) then redraw_flag = 1
         sMainState.sViewState.zScale = zFac
      end
      'wPerspSlid': begin
         perspDist = FLOAT(event.value)
         if (sMainState.sViewState.pDist ne perspDist) then redraw_flag = 1
         sMainState.sViewState.pDist = perspDist
      end
      'wViewDispBttn': begin
         sMainState.cleanupBuffer = 1
         Viz3D_CleanBuffer, sMainState
      end
   endcase

   if (redraw_flag) then begin
      sMainState.sViewState = $
         Viz3D_View(sMainState.sViewState.viewWin, $
            XMAX=sMainState.sViewState.xMax, $
            YMAX=sMainState.sViewState.yMax, $
            ZMAX=sMainState.sViewState.zMax, $
            ANG1=sMainState.sViewState.ang1, $
            ANG2=sMainState.sViewState.ang2, $
            ANG3=sMainState.sViewState.ang3, $
            DIR1=sMainState.sViewState.dir1, $
            DIR2=sMainState.sViewState.dir2, $
            DIR3=sMainState.sViewState.dir3, $
            ZOOM=sMainState.sViewState.zoomFac , $
            ZSCALE=sMainState.sViewState.zScale , $
            PERSP=sMainState.sViewState.pDist)
      Viz3D_ViewShow, sMainState
      sMainState.cleanupBuffer = 1
   endif

   WIDGET_CONTROL, wMainBase, SET_UVALUE=sMainState, /NO_COPY
end
;******************************************************************************


;******************************************************************************
; Differential shading event handler.
pro Viz3D_DiffShadeEvent, event
   WIDGET_CONTROL, /HOURGLASS

   WIDGET_CONTROL, event.id, GET_UVALUE=uVal

   case uVal of
      'wDiffSlid': begin
         WIDGET_CONTROL, event.top, GET_UVALUE=leader
         WIDGET_CONTROL, leader, GET_UVALUE=sMainState, /NO_COPY
         sMainState.sColorState.diffShade = event.value
         sViz3DColors = sMainState.sColorState
         Viz3D_DiffColor, sViz3DColors
         sMainState.sColorState = TEMPORARY(sViz3DColors)

         if (event.drag eq 0) then begin
            WSET, sMainState.mainWin
            DEVICE, COPY=[0, 0, sMainState.winX, sMainState.winY, 0, 0, $
                          sMainState.pixWin]
         endif

         WIDGET_CONTROL, leader, SET_UVALUE=sMainState, /NO_COPY
      end
      'wDiffCloseBttn': WIDGET_CONTROL, event.top, /DESTROY
   endcase
end
;******************************************************************************


;******************************************************************************
; Options event handler.
pro Viz3D_OptEvent, event
   WIDGET_CONTROL, /HOURGLASS

   WIDGET_CONTROL, event.id, GET_UVALUE=uVal
   WIDGET_CONTROL, event.top, GET_UVALUE=hState
   sOptState = *hState

   WIDGET_CONTROL, sOptState.wOptSizeText, GET_VALUE=szText
   DEVICE, GET_SCREEN_SIZE=scrSize
   wSize = (FIX(szText(0)) > 256) < (scrSize(1) - 80)
   WIDGET_CONTROL, sOptState.wOptSizeText, SET_VALUE=STRING(wSize)
   sOptState.winXY = wSize

   case uVal of
      'wOptAxisOnBttn': sOptState.AxisOn = 1
      'wOptAxisOffBttn': sOptState.AxisOn = 0
      'wOptCubeOnBttn': sOptState.CubeOn = 1
      'wOptCubeOffBttn': sOptState.CubeOn = 0
      'wOptSizeText':
      'wOptOkBttn': begin
         sOptState.Ok = 1
         WIDGET_CONTROL, event.top, /DESTROY
      end
      'wOptCanBttn': begin
         sOptState.Ok = 0
         WIDGET_CONTROL, event.top, /DESTROY
      end
   endcase

   *hState = sOptState
end
;******************************************************************************


;******************************************************************************
; Main event handler.
pro Viz3D_Event, event

ON_IOERROR, IO_BAD

   WIDGET_CONTROL, /HOURGLASS

   WIDGET_CONTROL, event.id, GET_UVALUE=uVal
   if (uVal eq 'wMainDraw') then begin
      wDrawBase = WIDGET_INFO(event.id, /PARENT)
      WIDGET_CONTROL, wDrawBase, GET_UVALUE=mainID
   endif else mainID = event.top

   WIDGET_CONTROL, mainID, GET_UVALUE=sMainState, /NO_COPY

   if ((uVal ne 'wQuitBttn') AND (uVal ne 'wDelAllBttn')) then $
      Viz3D_CleanView, sMainState

   case uVal of
      'wQuitBttn': begin ; Exit.
         WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
         WIDGET_CONTROL, mainID, /Destroy
         return
      end

      'wLoadBttn': begin ; Load new data from a file.
         ; Get file name.
         dataFile = PICKFILE(GROUP=event.top, /READ, FILTER='*.dat', $
                             /MUST_EXIST)
         WIDGET_CONTROL, /HOURGLASS
         if (dataFile(0) eq '') then begin ; No file, user cancelled.
            WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
            return
         endif

         ; Open file.
         GET_LUN, lun
         OPENR, lun, dataFile

         ; Read in and check the data.

         goRead = 1
         count = 0

         ; Keep reading as long as there is valid data.
         while (goRead) do begin
            count = count + 1

            fileStat = FSTAT(lun)

            if (fileStat.cur_ptr eq fileStat.size) then begin ; At EOF.
               WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
               CLOSE, lun
               FREE_LUN, lun
               return
            endif

            if ((fileStat.size - fileStat.cur_ptr) lt 29L) then begin
               ; Not enough bytes in file for data "header".
               ans = DIALOG_MESSAGE( $
                  'File does not contain valid 3-D data.', /ERROR)
               WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
               CLOSE, lun
               FREE_LUN, lun
               return
            endif

            ; Read and check the number of dimensions.
            nDims = 0L
            READU, lun, nDims
            if ((ndims NE 3L) and (count eq 1)) then begin
               nDims = SWAP_ENDIAN(TEMPORARY(nDims))
               if (ndims NE 3L) then begin
                  ans = DIALOG_MESSAGE( $
                     'File does not contain valid 3-D data.', /ERROR)
                  WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
                  CLOSE, lun
                  FREE_LUN, lun
                  return
               endif
               swapData = 1
            endif else swapData = 0
            if ((ndims NE 3L) and (count gt 1)) then begin
               ans = DIALOG_MESSAGE( $
                  'File does not contain valid 3-D data.', /ERROR)
               WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
               CLOSE, lun
               FREE_LUN, lun
               return
            endif

            ; Read the rest of the data "size" info.
            sizeDims = LONARR(3, /NOZERO)
            READU, lun, sizeDims
            dataType = 0L
            READU, lun, dataType
            nElems = 0L
            READU, lun, nElems

            ; Read the number of characters in data name.
            nNameChars = 0L
            READU, lun, nNameChars

            ; Byte swap if necessary.
            if (swapData) then begin
               sizeDims = SWAP_ENDIAN(TEMPORARY(sizeDims))
               dataType = SWAP_ENDIAN(TEMPORARY(dataType))
               nElems = SWAP_ENDIAN(TEMPORARY(nElems))
               nNameChars = SWAP_ENDIAN(TEMPORARY(nNameChars))
            endif

            ; Check the size of the data dimensions.
            if (MIN(sizeDims) LE 1L) then begin
               ans = DIALOG_MESSAGE( $
                  'File does not contain valid 3-D data.', /ERROR)
               WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
               CLOSE, lun
               FREE_LUN, lun
               return
            endif

            ; Check the data type.
            if ((dataType LT 1L) OR (dataType GT 5)) then begin
               ans = DIALOG_MESSAGE( $
                  'File does not contain valid 3-D data.', /ERROR)
               WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
               CLOSE, lun
               FREE_LUN, lun
               return
            endif

            ; Check the number of elements.
            if ((sizeDims(0)*sizeDims(1)*sizeDims(2)) NE nElems) then begin
               ans = DIALOG_MESSAGE( $
                  'File does not contain valid 3-D data.', /ERROR)
               WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
               CLOSE, lun
               FREE_LUN, lun
               return
            endif

            ; Check the number of characters in data name.
            if (nNameChars LT 1L) then begin
               ans = DIALOG_MESSAGE( $
                  'File contains invalid 3-D data name.', /ERROR)
               WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
               CLOSE, lun
               FREE_LUN, lun
               return
            endif

            ; See if there is enough bytes left in the file for data name.
            fileStat = FSTAT(lun)
            if ((fileStat.size - fileStat.cur_ptr) LT nNameChars) then begin
               ans = DIALOG_MESSAGE( $
                  'File does not contain valid 3-D data.', /ERROR)
               WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
               CLOSE, lun
               FREE_LUN, lun
               return
            endif

            ; Read the data name.
            dataName = BYTARR(nNameChars, /NOZERO)
            READU, lun, dataName
            dataName = STRING(dataName)

            case dataType of
               1: nBytes = nElems
               2: nBytes = nElems * 2L
               3: nBytes = nElems * 4L
               4: nBytes = nElems * 4L
               5: nBytes = nElems * 8L
            endcase

            fileStat = FSTAT(lun)

            ; If there are enough bytes left in the file, then read the data.
            if ((fileStat.size - fileStat.cur_ptr) ge nBytes) then begin

               if (count eq 1) then begin
                  ; Make the first data active.
                  sMainState.curData = 0
                  sMainState.sSurfState.curShade = 0

                  ; Something valid to read, so clean out the existing data.
                  lData3D = *(sMainState.hDataList)
                  for i=0, N_ELEMENTS(lData3D)-1L do begin
                     if (PTR_VALID(lData3D(i))) then $
                        PTR_FREE, lData3D(i)
                  endfor
                  lDataHist = *(sMainState.sThreshState.hDataHist)
                  for i=0, N_ELEMENTS(lDataHist)-1L do begin
                     if (PTR_VALID(lDataHist(i))) then $
                        PTR_FREE, lDataHist(i)
                  endfor
                  lRangeHist = *(sMainState.sSurfState.hRangeHist)
                  for i=0, N_ELEMENTS(lRangeHist)-1L do begin
                     if (PTR_VALID(lRangeHist(i))) then $
                        PTR_FREE, lRangeHist(i)
                  endfor

                  sizeDims1 = sizeDims
               endif else begin
                  ; Check if Nth data is the same size as the first data.
                  index = WHERE(sizeDims1 ne sizeDims)
                  if (index(0) ge 0L) then begin
                     ans = DIALOG_MESSAGE( $
                        'Data dimensions are not consistent.', /ERROR)
                     WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
                     CLOSE, lun
                     FREE_LUN, lun
                     return
                  endif
               endelse

               ; Read the data.
               case dataType of
                  1: data3D = BYTARR(sizeDims(0), sizeDims(1), sizeDims(2), $
                                 /NOZERO)
                  2: data3D = INTARR(sizeDims(0), sizeDims(1), sizeDims(2), $
                                 /NOZERO)
                  3: data3D = LONARR(sizeDims(0), sizeDims(1), sizeDims(2), $
                                 /NOZERO)
                  4: data3D = FLTARR(sizeDims(0), sizeDims(1), sizeDims(2), $
                                 /NOZERO)
                  5: data3D = DBLARR(sizeDims(0), sizeDims(1), sizeDims(2), $
                                 /NOZERO)
               endcase
               READU, lun, data3D
               if (swapData) then data3D = SWAP_ENDIAN(TEMPORARY(data3D))

               if (count eq 1) then begin

                  ; First data, reset everything.

                  minData = MIN(data3D, MAX=maxData)
                  minData = DOUBLE(minData)
                  maxData = DOUBLE(maxData)
                  if (minData eq maxData) then maxData(0) = maxData(0) + 1.0D
                  bSize = (maxData - minData) / 200.0
                  n_samples = N_ELEMENTS(data3D)
                  IF (n_samples GT 32767L) THEN BEGIN
                     n_samples = 32767L + ROUND(SQRT(n_samples - 32767L))
                     s_index = RANDOMU(s, 1)
                     s_index = RANDOMU(s, n_samples)
                     s_index = LONG(TEMPORARY(s_index) * FLOAT(n_samples))
                     hist_data = FLOAT(data3D[TEMPORARY(s_index)])
                     dataHist = $
                        HISTOGRAM(TEMPORARY(hist_data), MIN=minData, MAX=maxData, $
                           BINSIZE=bSize) > 1.0
                  ENDIF ELSE BEGIN
                     dataHist = $
                        HISTOGRAM(FLOAT(data3D), MIN=minData, MAX=maxData, $
                           BINSIZE=bSize) > 1.0
                  ENDELSE
                  dataHist = ALOG(FLOAT(TEMPORARY(dataHist))) > 0.0

                  lDataHist = PTR_NEW(dataHist)
                  lRangeHist = PTR_NEW(dataHist, /NO_COPY)
                  *(sMainState.sThreshState.hDataHist) = lDataHist
                  *(sMainState.sSurfState.hRangeHist) = lRangeHist

                  surfThresh = (minData + maxData) / 2.0D
                  *(sMainState.sSurfState.hSurfThresh) = surfThresh

                  szData = SIZE(data3D)
                  lData3D = PTR_NEW(data3D, /NO_COPY)
                  *(sMainState.hDataList) = lData3D

                  WIDGET_CONTROL, sMainState.wDataModeBase, SENSITIVE=0
                  WIDGET_CONTROL, sMainState.sSurfState.wSurfShadeBase, $
                     SENSITIVE=0
                  WIDGET_CONTROL, sMainState.wDataDrop, $
                     SET_DROPLIST_SELECT=0, SET_VALUE=dataName(0)
                  WIDGET_CONTROL, sMainState.sSurfState.wSurfShadeDrop, $
                     SET_DROPLIST_SELECT=0, SET_VALUE='Light-source'

                  *(sMainState.hDataName) = dataName
                  *(sMainState.sThreshState.hMinData) = minData
                  *(sMainState.sThreshState.hMaxData) = maxData
                  *(sMainState.sThreshState.hTransVal) = minData
                  *(sMainState.sThreshState.hLowThresh) = minData
                  *(sMainState.sThreshState.hHighThresh) = maxData

                  sMainState.szData = szData
                  sMainState.sViewState = $
                     Viz3D_View(sMainState.sViewState.viewWin, $
                        XMAX=(szData(1)-1), $
                        YMAX=(szData(2)-1), $
                        ZMAX=(szData(3)-1), $
                        ANG1=sMainState.sViewState.ang1, $
                        ANG2=sMainState.sViewState.ang2, $
                        ANG3=sMainState.sViewState.ang3, $
                        DIR1=sMainState.sViewState.dir1, $
                        DIR2=sMainState.sViewState.dir2, $
                        DIR3=sMainState.sViewState.dir3, $
                        ZOOM=sMainState.sViewState.zoomFac , $
                        ZSCALE=sMainState.sViewState.zScale , $
                        PERSP=sMainState.sViewState.pDist)

                  sMainState.sBlockState.c1 = $
                     [szData(1)-1, szData(2)-1, szData(3)-1] / 3
                  sMainState.sBlockState.c2 = $
                     (2 * [szData(1)-1, szData(2)-1, szData(3)-1]) / 3

                  sMainState.sProfState.x1Ortho = FLOAT(szData(1) - 1L) / 2.0
                  sMainState.sProfState.y1Ortho = FLOAT(szData(2) - 1L) / 2.0
                  sMainState.sProfState.z1Ortho = 0.0
                  sMainState.sProfState.x2Ortho = FLOAT(szData(1) - 1L) / 2.0
                  sMainState.sProfState.y2Ortho = FLOAT(szData(2) - 1L) / 2.0
                  sMainState.sProfState.z2Ortho = FLOAT(szData(3) - 1L)
                  sMainState.sProfState.x1Obliq = FLOAT(szData(1) - 1L) / 2.0
                  sMainState.sProfState.y1Obliq = FLOAT(szData(2) - 1L) / 2.0
                  sMainState.sProfState.z1Obliq = 0.0
                  sMainState.sProfState.x2Obliq = FLOAT(szData(1) - 1L) / 2.0
                  sMainState.sProfState.y2Obliq = FLOAT(szData(2) - 1L) / 2.0
                  sMainState.sProfState.z2Obliq = FLOAT(szData(3) - 1L)

                  sMainState.sProbeState.x = FLOAT(szData(1) - 1L) / 2.0
                  sMainState.sProbeState.y = FLOAT(szData(2) - 1L) / 2.0
                  sMainState.sProbeState.z = FLOAT(szData(3) - 1L) / 2.0

                  ; Clean out the display list.
                  for i=0, N_ELEMENTS(*(sMainState.hDisplayList)) -1 do begin
                     hGraphic = (*(sMainState.hDisplayList))[i]
                     IF (TAG_NAMES( *hGraphic, /STRUCTURE_NAME ) EQ 'VIZ3D_SURF') $
                     THEN BEGIN
                        PTR_FREE, (*hGraphic).hVertList
                        PTR_FREE, (*hGraphic).hPolyList
                        PTR_FREE, (*hGraphic).hColorList
                     ENDIF
                     IF (TAG_NAMES( *hGraphic, /STRUCTURE_NAME ) EQ $
                        'VIZ3D_OBLIQ_PLANE') THEN BEGIN
                        PTR_FREE, (*hGraphic).hVertPlane
                        PTR_FREE, (*hGraphic).hVert3D
                     ENDIF
                     PTR_FREE, hGraphic
;                    PTR_FREE, (*(sMainState.hDisplayList))[i]
                  endfor
                  PTR_FREE, sMainState.hDisplayList
                  sMainState.hDisplayList = PTR_NEW(/ALLOCATE_HEAP)
                  ; Repair "Delete" menu.
                  delBttn = WIDGET_INFO(sMainState.wDeleteMenu, /CHILD)
                  WHILE (WIDGET_INFO(delBttn, /VALID_ID)) DO BEGIN
                     WIDGET_CONTROL, delBttn, /DESTROY
                     delBttn = WIDGET_INFO(sMainState.wDeleteMenu, /CHILD)
                  ENDWHILE
                  ; Erase.
                  SET_PLOT, 'Z'
                  ERASE
                  Viz3D_DrawData, sMainState
                  sMainState.cleanupBuffer = 1
                  Viz3D_ThreshShow, sMainState
                  Viz3D_SurfShow, sMainState, /CALC_HIST
                  Viz3D_CleanBuffer, sMainState
                  WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=''
                  WIDGET_CONTROL, sMainState.sProbeState.wProbeXText, $
                     SET_VALUE=''
                  WIDGET_CONTROL, sMainState.sProbeState.wProbeYText, $
                     SET_VALUE=''
                  WIDGET_CONTROL, sMainState.sProbeState.wProbeZText, $
                     SET_VALUE=''

               endif else begin

                  ; Additional data.

                  minData = MIN(data3D, MAX=maxData)
                  minData = DOUBLE(minData)
                  maxData = DOUBLE(maxData)
                  if (minData eq maxData) then maxData(0) = maxData(0) + 1.0D

                  bSize = (maxData - minData) / 200.0
                  n_samples = N_ELEMENTS(data3D)
                  IF (n_samples GT 32767L) THEN BEGIN
                     n_samples = 32767L + ROUND(SQRT(n_samples - 32767L))
                     s_index = RANDOMU(s, 1)
                     s_index = RANDOMU(s, n_samples)
                     s_index = LONG(TEMPORARY(s_index) * FLOAT(n_samples))
                     hist_data = FLOAT(data3D[TEMPORARY(s_index)])
                     dataHist = $
                        HISTOGRAM(TEMPORARY(hist_data), MIN=minData, MAX=maxData, $
                           BINSIZE=bSize) > 1.0
                  ENDIF ELSE BEGIN
                     dataHist = $
                        HISTOGRAM(FLOAT(data3D), MIN=minData, MAX=maxData, $
                           BINSIZE=bSize) > 1.0
                  ENDELSE

                  dataHist = ALOG(FLOAT(TEMPORARY(dataHist))) > 0.0
                  lDataHist = *(sMainState.sThreshState.hDataHist)
                  lDataHist = [TEMPORARY(lDataHist), $
                     PTR_NEW(dataHist)]
                  *(sMainState.sThreshState.hDataHist) = lDataHist
                  lRangeHist = *(sMainState.sSurfState.hRangeHist)
                  lRangeHist = [TEMPORARY(lRangeHist), $
                     PTR_NEW(dataHist, /NO_COPY)]
                  *(sMainState.sSurfState.hRangeHist) = lRangeHist

                  lData3D = *(sMainState.hDataList)
                  lData3D = [TEMPORARY(lData3D), $
                     PTR_NEW(data3D, /NO_COPY)]
                  *(sMainState.hDataList) = lData3D
                  dataNames = *(sMainState.hDataName)
                  dataNames = [TEMPORARY(dataNames), dataName]

                  WIDGET_CONTROL, sMainState.wDataModeBase, SENSITIVE=1
                  WIDGET_CONTROL, sMainState.sSurfState.wSurfShadeBase, $
                     SENSITIVE=1
                  WIDGET_CONTROL, sMainState.wDataDrop, $
                     SET_DROPLIST_SELECT=0, SET_VALUE=dataNames
                  sName = dataNames
                  sName(0) = 'Light-source'
                  WIDGET_CONTROL, sMainState.sSurfState.wSurfShadeDrop, $
                     SET_DROPLIST_SELECT=0, SET_VALUE=sName

                  *(sMainState.hDataName) = dataNames

                  lMinData = *(sMainState.sThreshState.hMinData)
                  lMaxData = *(sMainState.sThreshState.hMaxData)
                  lMinData = [TEMPORARY(lMinData), minData]
                  lMaxData = [TEMPORARY(lMaxData), maxData]

                  surfThresh = (lMinData + lMaxData) / 2.0D
                  *(sMainState.sSurfState.hSurfThresh) = surfThresh

                  *(sMainState.sThreshState.hTransVal) = lMinData
                  *(sMainState.sThreshState.hLowThresh) = lMinData
                  *(sMainState.sThreshState.hHighThresh) = lMaxData

                  *(sMainState.sThreshState.hMinData) = lMinData
                  *(sMainState.sThreshState.hMaxData) = lMaxData

               endelse
            endif else begin
               goRead = 0
            endelse

         endwhile

         CLOSE, lun
         FREE_LUN, lun

      end

      'wSaveTiffBttn': begin ; Save tiff image.
         tiffFile = PICKFILE(GROUP=event.top, /WRITE, FILTER='*.tif', $
                             FILE='slicer3.tif')
         WIDGET_CONTROL, /HOURGLASS
         if (tiffFile(0) eq '') then begin
            WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
            return
         endif
         tiffFiles = FINDFILE(tiffFile)
         if (tiffFiles(0) ne '') then begin
            ans = DIALOG_MESSAGE('File exists.   Overwrite ?', /QUESTION)
            WIDGET_CONTROL, /HOURGLASS
            if (ans ne 'Yes') then begin
               WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
               return
            endif
         endif
         WSET, sMainState.pixWin
         if (sMainState.sColorState.displayBits eq 24) then begin
            iR = TVRD(CHANNEL=1)
            iG = TVRD(CHANNEL=2)
            iB = TVRD(CHANNEL=3)

            TIFF_WRITE, tiffFile(0), PLANARCONFIG=2, $
               RED=TEMPORARY(iR), GREEN=TEMPORARY(iG), BLUE=TEMPORARY(iB)
         endif else begin
            image = TVRD()
            TVLCT, cR, cG, cB, /GET
            TIFF_WRITE, tiffFile(0), TEMPORARY(image), $
                        RED=cR, GREEN=cG, BLUE=cB
         endelse
         WSET, sMAinState.mainWin
      end
      'wSavePSBttn': begin ; Save Postscript image.
         psFile = PICKFILE(GROUP=event.top, /WRITE, FILTER='*.ps', $
                             FILE='slicer3.ps')
         WIDGET_CONTROL, /HOURGLASS
         if (psFile(0) eq '') then begin
            WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
            return
         endif
         psFiles = FINDFILE(psFile)
         if (psFiles(0) ne '') then begin
            ans = DIALOG_MESSAGE('File exists.   Overwrite ?', /QUESTION)
            WIDGET_CONTROL, /HOURGLASS
            if (ans ne 'Yes') then begin
               WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
               return
            endif
         endif
         WSET, sMainState.pixWin
         image = TVRD()
         set_plot,'ps'
         device,filename=psFile(0),/color,bits=8, $
                xsize=20,ysize=20,/encapsulated
;                xsize=16,ysize=16,xoffset=2,yoffset=2
         image(where(image eq 0B))=255B
         TV,image
         device,/close
         set_plot,'x'
         WSET, sMAinState.mainWin
      end

      'wSaveSubsetBttn': begin ; Save subset.
         dataFile = PICKFILE(GROUP=event.top, /WRITE, FILTER='*.dat', $
                             FILE='slicer3.dat')
         WIDGET_CONTROL, /HOURGLASS
         if (dataFile(0) eq '') then begin
            WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
            return
         endif
         dataFiles = FINDFILE(dataFile)
         if (dataFiles(0) ne '') then begin
            ans = DIALOG_MESSAGE('File exists.   Overwrite ?', /QUESTION)
            WIDGET_CONTROL, /HOURGLASS
            if (ans ne 'Yes') then begin
               WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
               return
            endif
         endif

         x1 = sMainState.sBlockState.c1(0) < sMainState.sBlockState.c2(0)
         y1 = sMainState.sBlockState.c1(1) < sMainState.sBlockState.c2(1)
         z1 = sMainState.sBlockState.c1(2) < sMainState.sBlockState.c2(2)
         x2 = sMainState.sBlockState.c1(0) > sMainState.sBlockState.c2(0)
         y2 = sMainState.sBlockState.c1(1) > sMainState.sBlockState.c2(1)
         z2 = sMainState.sBlockState.c1(2) > sMainState.sBlockState.c2(2)

         GET_LUN, lun
         OPENW, lun, dataFile

         lData3D = *(sMainState.hDataList)
         dataNames = *(sMainState.hDataName)

         for i=0, N_ELEMENTS(lData3D)-1L do begin
            WRITEU, lun, SIZE((*(lData3D(i)))(x1:x2, y1:y2, z1:z2))
            WRITEU, lun, STRLEN(dataNames(i))
            WRITEU, lun, BYTE(dataNames(i))
            WRITEU, lun, (*(lData3D(i)))(x1:x2, y1:y2, z1:z2)
         endfor

         *(sMainState.hDataName) = dataNames
         *(sMainState.hDataList) = lData3D

         CLOSE, lun
         FREE_LUN, lun
      end

      'wOptionsBttn': begin ; Options.
         wOptBase = WIDGET_BASE(TITLE='Slicer3', GROUP_LEADER=event.top, $
                       /COLUMN, XPAD=1, YPAD=1, SPACE=1, /MODAL)

         wOptAxisBase = WIDGET_BASE(wOptBase, /ROW, /FRAME, $
                                    XPAD=1, YPAD=1, SPACE=1)
         wOptAxisLab = WIDGET_LABEL(wOptAxisBase, VALUE='Axis')
         wOptAxisExBase = WIDGET_BASE(wOptAxisBase, /ROW, /EXCLUSIVE, $
                                      XPAD=1, YPAD=1, SPACE=1)
         wOptAxisOnBttn = WIDGET_BUTTON(wOptAxisExBase, VALUE='On', $
                             UVALUE='wOptAxisOnBttn', /NO_RELEASE)
         wOptAxisOffBttn = WIDGET_BUTTON(wOptAxisExBase, VALUE='Off', $
                              UVALUE='wOptAxisOffBttn', /NO_RELEASE)
         if (sMainState.axisOn) then $
            WIDGET_CONTROL, wOptAxisOnBttn, /SET_BUTTON $
         else WIDGET_CONTROL, wOptAxisOffBttn, /SET_BUTTON

         wOptCubeBase = WIDGET_BASE(wOptBase, /ROW, /FRAME, $
                                    XPAD=1, YPAD=1, SPACE=1)
         wOptCubeLab = WIDGET_LABEL(wOptCubeBase, VALUE='Cube')
         wOptCubeExBase = WIDGET_BASE(wOptCubeBase, /ROW, /EXCLUSIVE, $
                                      XPAD=1, YPAD=1, SPACE=1)
         wOptCubeOnBttn = WIDGET_BUTTON(wOptCubeExBase, VALUE='On', $
                             UVALUE='wOptCubeOnBttn', /NO_RELEASE)
         wOptCubeOffBttn = WIDGET_BUTTON(wOptCubeExBase, VALUE='Off', $
                              UVALUE='wOptCubeOffBttn', /NO_RELEASE)
         if (sMainState.cubeOn) then $
            WIDGET_CONTROL, wOptcubeOnBttn, /SET_BUTTON $
         else WIDGET_CONTROL, wOptcubeOffBttn, /SET_BUTTON

         wOptSizeBase = WIDGET_BASE(wOptBase, /ROW, XPAD=1, YPAD=1, $
                                    SPACE=1)
         wOptSizeLab = WIDGET_LABEL(wOptSizeBase, VALUE='Window Size')
         wOptSizeText = WIDGET_TEXT(wOptSizeBase, UVALUE='wOptSizeText', $
                           VALUE=STRING(sMAinState.winX), /EDITABLE)

         wOptBttnBase = WIDGET_BASE(wOptBase, /ROW, XPAD=1, YPAD=1, $
                                    SPACE=1)
         wOptOkBttn = WIDGET_BUTTON(wOptBttnBase, VALUE='Ok', $
                         UVALUE='wOptOkBttn', /ALIGN_LEFT)
         wOptCanBttn = WIDGET_BUTTON(wOptBttnBase, VALUE='Cancel', $
                          UVALUE='wOptCanBttn', /ALIGN_RIGHT)

         sOptState = {SOptState, Ok:0, axisOn:sMainState.axisOn, $
                      cubeOn:sMainState.cubeOn, winXY:sMainState.winX, $
                      wOptSizeText:wOptSizeText}
         hState = PTR_NEW(sOptState, /NO_COPY)
         WIDGET_CONTROL, wOptBase, SET_UVALUE=hState

         WIDGET_CONTROL, wOptBase, /REALIZE
         XMANAGER, 'Viz3D_Event', wOptBase, $
            EVENT_HANDLER='Viz3D_OptEvent', /MODAL
         WIDGET_CONTROL, /HOURGLASS

         sOptState = *hState
         PTR_FREE, hState
         if (sOptState.Ok) then begin ; User clicked Ok.

            if ((sOptState.winXY ne sMainState.winX) or $
                (sOptState.winXY ne sMainState.winY)) then begin
               ; Resize window.
               sMainState.winX = sOptState.winXY
               sMainState.winY = sOptState.winXY
               sMainState.cubeOn = sOptState.cubeOn
               sMainState.axisOn = sOptState.axisOn

               SET_PLOT, 'Z'
               DEVICE, SET_RESOLUTION=[sMainState.winX, sMainState.winY], $
                       /Z_BUFFERING
               SET_PLOT, sMainState.screenDevice
               WDELETE, sMainState.pixWin
               WINDOW, /FREE, /PIXMAP, XSIZE=sMainState.winX, $
                  YSIZE=sMainState.winY, RETAIN=2
               sMainState.pixWin = !D.Window

               WIDGET_CONTROL, sMainState.wMainDraw, /DESTROY
               sMainState.wMainDraw = WIDGET_DRAW(sMainState.wDrawBase, $
                  XSIZE=sMainState.winX, YSIZE=sMainState.winY, RETAIN=2, $
                  UVALUE='wMainDraw', /BUTTON_EVENTS)
               WIDGET_CONTROL, sMainState.wMainDraw, GET_VALUE=mainWin
               sMainState.mainWin = mainWin
               WSET, sMainState.mainWin

               Viz3D_FillDepth, sMainState
               Viz3D_Redraw, sMainState
            endif else begin
               if (((sMainState.cubeOn eq 1) and (sOptState.cubeOn eq 0)) or $
                   ((sMainState.axisOn eq 1) and (sOptState.axisOn eq 0))) $
               then begin
                  sMainState.cubeOn = sOptState.cubeOn
                  sMainState.axisOn = sOptState.axisOn
                  Viz3D_Redraw, sMainState
               endif else begin
                  sMainState.cubeOn = sOptState.cubeOn
                  sMainState.axisOn = sOptState.axisOn
                  Viz3D_DrawData, sMainState
               endelse
            endelse
         endif
      end

      'wEraseBttn': begin ; Erase.
         ; Clean out the display list.
         for i=0, N_ELEMENTS(*(sMainState.hDisplayList)) -1 do begin
            hGraphic = (*(sMainState.hDisplayList))[i]
            IF (TAG_NAMES( *hGraphic, /STRUCTURE_NAME ) EQ 'VIZ3D_SURF') $
            THEN BEGIN
               PTR_FREE, (*hGraphic).hVertList
               PTR_FREE, (*hGraphic).hPolyList
               PTR_FREE, (*hGraphic).hColorList
            ENDIF
            IF (TAG_NAMES( *hGraphic, /STRUCTURE_NAME ) EQ $
               'VIZ3D_OBLIQ_PLANE') THEN BEGIN
               PTR_FREE, (*hGraphic).hVertPlane
               PTR_FREE, (*hGraphic).hVert3D
            ENDIF
            PTR_FREE, hGraphic
         endfor
         PTR_FREE, sMainState.hDisplayList
         sMainState.hDisplayList = PTR_NEW(/ALLOCATE_HEAP)
         if (sMainState.actMode eq 7) then begin
            ; View mode, redraw small slice and block windows.
            Viz3D_SliceShow, sMainState
            Viz3D_BlockShow, sMainState
         endif
         ; Repair "Delete" menu.
         delBttn = WIDGET_INFO(sMainState.wDeleteMenu, /CHILD)
         WHILE (WIDGET_INFO(delBttn, /VALID_ID)) DO BEGIN
            WIDGET_CONTROL, delBttn, /DESTROY
            delBttn = WIDGET_INFO(sMainState.wDeleteMenu, /CHILD)
         ENDWHILE
         ; Erase.
         SET_PLOT, 'Z'
         ERASE
         Viz3D_DrawData, sMainState
      end

      'wDeleteBttn': begin ; Delete a graphic.
         ; Determine which button in the menu was chosen
         widBttn = WIDGET_INFO(sMainState.wDeleteMenu, /Child)
         hGraphic = (*(sMainState.hDisplayList))[0]
         indx = 0
         WHILE (widBttn ne event.id) do begin
            widBttn = WIDGET_INFO(widBttn, /SIBLING)
            indx = indx +1
            hGraphic = (*(sMainState.hDisplayList))[indx]
         ENDWHILE
         ; Delete and redraw.
         IF (TAG_NAMES( *hGraphic, /STRUCTURE_NAME ) $
             EQ 'VIZ3D_SURF') THEN BEGIN
            PTR_FREE, (*(hGraphic)).hVertList
            PTR_FREE, (*(hGraphic)).hPolyList
            PTR_FREE, (*(hGraphic)).hColorList
         ENDIF
         IF (TAG_NAMES( *hGraphic, /STRUCTURE_NAME ) EQ $
            'VIZ3D_OBLIQ_PLANE') THEN BEGIN
            PTR_FREE, (*hGraphic).hVertPlane
            PTR_FREE, (*hGraphic).hVert3D
         ENDIF
         PTR_FREE, hGraphic
         listCount = N_ELEMENTS(*(sMainState.hDisplayList))
         if (listCount GT 1) then begin
           if (indx EQ 0) then begin
              *(sMainState.hDisplayList) = $
                 (*(sMainState.hDisplayList))[1:listCount-1]
           endif else if (indx EQ listCount-1) then begin
              *(sMainState.hDisplayList) = $
                 (*(sMainState.hDisplayList))[0:listCount-2]
           endif else begin
              ; general case
              *(sMainState.hDisplayList) = $
                 (*(sMainState.hDisplayList))[[indgen(indx), $
                              indgen(listCount-indx-1)+indx+1]]
           endelse
         ; the list is now empty, the listpointer itself points to undefined
         endif else begin
            PTR_FREE, sMainState.hDisplayList
            sMainState.hDisplayList = PTR_NEW(/ALLOCATE_HEAP)
         endelse
         WIDGET_CONTROL, widBttn, /DESTROY
         Viz3D_Redraw, sMainState
      end

      'wColResetBttn': begin ; Reset Colors.
         sMainState.sColorState = Viz3D_LoadColor(DIFFSHADE=20)
         WSET, sMainState.mainWin
         Viz3D_DrawData, sMainState
      end

      'wColDiffBttn': begin ; Differential shading.
         diffBase = WIDGET_BASE(TITLE='Slicer3', XPAD=1, YPAD=1, SPACE=1, $
            /COLUMN, UVALUE=event.top, GROUP_LEADER=event.top, /MODAL)
         wDiffSlid = WIDGET_SLIDER(diffBase, MIN=0, MAX=100, /DRAG, $
                        TITLE='Differential Shading (%)', $
                        VALUE=sMainState.sColorState.diffShade, $
                        UVALUE='wDiffSlid')
         wDiffCloseBttn = WIDGET_BUTTON(diffBase, VALUE='Close', $
                             UVALUE='wDiffCloseBttn')
         WIDGET_CONTROL, diffBase, /REALIZE
         WIDGET_CONTROL, event.top, SET_UVALUE=sMainState, /NO_COPY
         XMANAGER, 'Viz3D_Event', diffBase, $
            EVENT_HANDLER='Viz3D_DiffShadeEvent', /MODAL
         WIDGET_CONTROL, event.top, GET_UVALUE=sMainState, /NO_COPY
         Viz3D_DrawData, sMainState
      end
      'wColSliceBttn': begin ; Slice/block colors.
         saveP = !P
         saveX = !X
         saveY = !Y
         saveZ = !Z
         XLOADCT, BOTTOM=0, /USE_CURRENT, $
            NCOLORS=sMainState.sColorState.nColor, GROUP=event.top, /MODAL
         WIDGET_CONTROL, /HOURGLASS
         WSET, sMainState.mainWin
         !P = TEMPORARY(saveP)
         !X = TEMPORARY(saveX)
         !Y = TEMPORARY(saveY)
         !Z = TEMPORARY(saveZ)
         sViz3DColors = sMainState.sColorState
         TVLCT, r, g, b, /GET
         sViz3DColors.rSlice = r(0:sMainState.sColorState.nColor-1)
         sViz3DColors.gSlice = g(0:sMainState.sColorState.nColor-1)
         sViz3DColors.bSlice = b(0:sMainState.sColorState.nColor-1)
         Viz3D_DiffColor, sViz3DColors
         sMainState.sColorState = TEMPORARY(sViz3DColors)
         if (sMainState.sColorState.displayBits eq 24) then LOADCT, 0
         Viz3D_DrawData, sMainState
      end
      'wColSurfBttn': begin ; Surface colors.
         saveP = !P
         saveX = !X
         saveY = !Y
         saveZ = !Z
         XLOADCT, BOTTOM=(3*sMainState.sColorState.nColor), /USE_CURRENT, $
            NCOLORS=sMainState.sColorState.nColor, GROUP=event.top, /MODAL
         WIDGET_CONTROL, /HOURGLASS
         WSET, sMainState.mainWin
         !P = TEMPORARY(saveP)
         !X = TEMPORARY(saveX)
         !Y = TEMPORARY(saveY)
         !Z = TEMPORARY(saveZ)
         if (sMainState.sColorState.displayBits eq 24) then LOADCT, 0
         Viz3D_DrawData, sMainState
      end
      'wColProjBttn': begin ; Projection colors.
         saveP = !P
         saveX = !X
         saveY = !Y
         saveZ = !Z
         XLOADCT, BOTTOM=(4*sMainState.sColorState.nColor), /USE_CURRENT, $
            NCOLORS=sMainState.sColorState.nColor, GROUP=event.top, /MODAL
         WIDGET_CONTROL, /HOURGLASS
         WSET, sMainState.mainWin
         !P = TEMPORARY(saveP)
         !X = TEMPORARY(saveX)
         !Y = TEMPORARY(saveY)
         !Z = TEMPORARY(saveZ)
         if (sMainState.sColorState.displayBits eq 24) then LOADCT, 0
         Viz3D_DrawData, sMainState
      end

      'wAnimateBttn': begin ; Animate
         ; *** Future animation functionality.
      end

      'wHelpBttn': begin ; Help
         XDISPLAYFILE, FILEPATH("slicer3.txt", $
                SUBDIR=['help', 'widget']), $
                TITLE="Slicer3 Help", $
                GROUP=event.top, $
                /MODAL, $
                WIDTH=72, HEIGHT=24
      end

      'wDataDrop': begin ; Set data.
         if (event.index ne sMainState.curData) then begin
            sMainState.curData = event.index
            sMainState.sSurfState.curShade = event.index
            dName = *(sMainState.hDataName)
            sName = dName
            sName(event.index) = 'Light-source'
            WIDGET_CONTROL, sMainState.sSurfState.wSurfShadeDrop, $
               SET_VALUE=sName, /NO_COPY
            WIDGET_CONTROL, sMainState.sSurfState.wSurfShadeDrop, $
               SET_DROPLIST_SELECT=event.index
            Viz3D_SurfShow, sMainState, /CALC_HIST
            Viz3D_ThreshShow, sMainState
         endif
      end

      'wModeDrop': begin ; Set mode.
         case event.index of
            0: begin ; Slice mode.
               if (sMainState.actMode ne 0) then begin
                  WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=''
                  WIDGET_CONTROL, sMainState.wCurMapBase, MAP=0
                  WIDGET_CONTROL, sMainState.wSliceBase, /MAP
                  sMainState.wCurMapBase = sMainState.wSliceBase
                  sMainState.actMode = 0
                  WIDGET_CONTROL, sMainState.wDrawBase, $
                     EVENT_PRO='Viz3D_SliceEvent'
                  WIDGET_CONTROL, sMainState.wSaveSubsetBttn, SENSITIVE=0
               endif
            end
            1: begin ; Block mode.
               if (sMainState.actMode ne 1) then begin
                  WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=''
                  WIDGET_CONTROL, sMainState.wCurMapBase, MAP=0
                  WIDGET_CONTROL, sMainState.wBlockBase, /MAP
                  sMainState.wCurMapBase = sMainState.wBlockBase
                  sMainState.actMode = 1
                  WIDGET_CONTROL, sMainState.wDrawBase, $
                     EVENT_PRO='Viz3D_BlockEvent'
                  Viz3D_BlockShow, sMainState
                  Viz3D_BlockShow, sMainState, /DIRECT

                  WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=('(' + $
                     STRTRIM(STRING(sMainState.sBlockState.c1(0)),2) + $
                     ', ' + $
                     STRTRIM(STRING(sMainState.sBlockState.c1(1)),2) + $
                     ', ' + $
                     STRTRIM(STRING(sMainState.sBlockState.c1(2)),2) + $
                     '), (' + $
                     STRTRIM(STRING(sMainState.sBlockState.c2(0)),2) + $
                     ', ' + $
                     STRTRIM(STRING(sMainState.sBlockState.c2(1)),2) + $
                     ', ' + $
                     STRTRIM(STRING(sMainState.sBlockState.c2(2)),2) + ')')

                  sMainState.cleanupView = 1
                  WIDGET_CONTROL, sMainState.wSaveSubsetBttn, SENSITIVE=1
               endif
            end
            2: begin ; Surface mode.
               if (sMainState.actMode ne 2) then begin
                  WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=''
                  WIDGET_CONTROL, sMainState.wCurMapBase, MAP=0
                  WIDGET_CONTROL, sMainState.wSurfBase, /MAP
                  sMainState.wCurMapBase = sMainState.wSurfBase
                  sMainState.actMode = 2
                  WIDGET_CONTROL, sMainState.wDrawBase, $
                     EVENT_PRO='Viz3D_SurfEvent'
                  WIDGET_CONTROL, sMainState.wSaveSubsetBttn, SENSITIVE=0
               endif
            end
            3: begin ; Projection mode.
               if (sMainState.actMode ne 4) then begin
                  WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=''
                  WIDGET_CONTROL, sMainState.wCurMapBase, MAP=0
                  WIDGET_CONTROL, sMainState.wProjectBase, /MAP
                  sMainState.wCurMapBase = sMainState.wProjectBase
                  sMainState.actMode = 4
                  WIDGET_CONTROL, sMainState.wDrawBase, $
                     EVENT_PRO='Viz3D_ProjectEvent'
                  WIDGET_CONTROL, sMainState.wSaveSubsetBttn, SENSITIVE=0
               endif
            end
            4: begin ; Threshold mode.
               if (sMainState.actMode ne 3) then begin
                  WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=''
                  WIDGET_CONTROL, sMainState.wCurMapBase, MAP=0
                  WIDGET_CONTROL, sMainState.wThreshBase, /MAP
                  sMainState.wCurMapBase = sMainState.wThreshBase
                  sMainState.actMode = 3
                  WIDGET_CONTROL, sMainState.wDrawBase, $
                     EVENT_PRO='Viz3D_ThreshEvent'
                  WIDGET_CONTROL, sMainState.wSaveSubsetBttn, SENSITIVE=0
               endif
            end
            5: begin ; Profile mode.
               if (sMainState.actMode ne 5) then begin
                  WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=''
                  WIDGET_CONTROL, sMainState.wCurMapBase, MAP=0
                  WIDGET_CONTROL, sMainState.wProfileBase, /MAP
                  sMainState.wCurMapBase = sMainState.wProfileBase
                  sMainState.actMode = 5
                  WIDGET_CONTROL, sMainState.wDrawBase, $
                     EVENT_PRO='Viz3D_ProfileEvent'
                  Viz3D_ProfShow, sMainState
                  Viz3D_ProfDraw, sMainState
                  sMainState.cleanupView = 1
                  WIDGET_CONTROL, sMainState.wSaveSubsetBttn, SENSITIVE=0
               endif
            end
            6: begin ; Probe mode.
               if (sMainState.actMode ne 6) then begin
                  WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=''
                  WIDGET_CONTROL, sMainState.wCurMapBase, MAP=0
                  WIDGET_CONTROL, sMainState.wProbeBase, /MAP
                  sMainState.wCurMapBase = sMainState.wProbeBase
                  sMainState.actMode = 6
                  WIDGET_CONTROL, sMainState.wDrawBase, $
                     EVENT_PRO='Viz3D_ProbeEvent'
                  sMainState.cleanupView = 1
                  WIDGET_CONTROL, sMainState.wSaveSubsetBttn, SENSITIVE=0
                  Viz3D_ProbeDraw, sMainState
               endif
            end
            7: begin ; View mode.
               if (sMainState.actMode ne 7) then begin
                  WIDGET_CONTROL, sMainState.wStatText, SET_VALUE=''
                  WIDGET_CONTROL, sMainState.wCurMapBase, MAP=0
                  WIDGET_CONTROL, sMainState.wViewBase, /MAP
                  sMainState.wCurMapBase = sMainState.wViewBase
                  sMainState.actMode = 7
                  WIDGET_CONTROL, sMainState.wDrawBase, $
                     EVENT_PRO='Viz3D_ViewEvent'
                  WIDGET_CONTROL, sMainState.wSaveSubsetBttn, SENSITIVE=0
               endif
            end
         endcase
      end
      else:
   endcase

   WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY
   return

IO_BAD:

   if (NOT(PTR_VALID(*(sMainState.hDataList)))) then begin
      ans = DIALOG_MESSAGE(['Unrecoverable I/O error.', $
                            'Exiting Slicer3.'], /ERROR)
      PTR_FREE, sMainState.hMain
      WIDGET_CONTROL, event.top, /DESTROY
      return
   endif

   if (N_ELEMENTS(sMainState) gt 0L) then $
      WIDGET_CONTROL, mainID, SET_UVALUE=sMainState, /NO_COPY

   if (n_elements(fLun) gt 0L) then begin
      CLOSE, fLun
      FREE_LUN, fLun
   endif

   ans = DIALOG_MESSAGE('Unknown I/O error.', /ERROR)

end
;******************************************************************************


;******************************************************************************
; Procedure to clean up the main draw window if it needs it.

pro Viz3D_CleanView, sMainState, SKIP_BUFFER=skipBuffer

   if (not(KEYWORD_SET(skipBuffer))) then $
      Viz3D_CleanBuffer, sMainState

   if (sMainState.cleanupView) then begin
      WSET, sMainState.mainWin
      DEVICE, COPY=[0, 0, sMainState.winX, sMainState.winY, $
                    0, 0, sMainState.pixWin]
      sMainState.cleanupView = 0
   endif
end
;******************************************************************************


;******************************************************************************
; Procedure to clean up (redraw) the Z buffer (and other windows) if needed.

pro Viz3D_CleanBuffer, sMainState
   if (sMainState.cleanupBuffer) then begin
      Viz3D_FillDepth, sMainState
      Viz3D_SliceShow, sMainState
      Viz3D_BlockShow, sMainState
      Viz3D_SurfShow, sMainState
      Viz3D_ThreshShow, sMainState, /DYNAMIC
      Viz3D_ViewShow, sMainState
      Viz3D_ProfShow, sMainState
      Viz3D_Redraw, sMainState
   endif
end
;******************************************************************************


;******************************************************************************
; Cleanup handler for detached base.

pro Viz3D_KillDraw, wDrawBase

   WIDGET_CONTROL, wDrawBase, GET_UVALUE=wMainBase, /NO_COPY


   WIDGET_CONTROL, wMainBase, /DESTROY

;  if (WIDGET_INFO(wMainBase, /VALID_ID)) then begin
;     WIDGET_CONTROL, wMainBase, GET_UVALUE=sMainState, /NO_COPY
;     WIDGET_CONTROL, wMainBase, KILL_NOTIFY=''
;     WIDGET_CONTROL, wMainBase, /DESTROY
;     PTR_FREE, sMainState.hMain
;     WDELETE, sMainState.pixWin
;     WDELETE, sMainState.sProfState.profPix
;     WDELETE, sMainState.sSurfState.surfPix
;     WDELETE, sMainState.sThreshState.threshPix
;  endif

end
;******************************************************************************


;******************************************************************************
; Cleanup handler for main window.

pro Viz3D_KillMain, wMainBase

   WIDGET_CONTROL, wMainBase, GET_UVALUE=sMainState, /NO_COPY
   if (sMainState.detachBase) then begin
      IF WIDGET_INFO(sMainState.wDrawBase, /VALID) THEN BEGIN
         WIDGET_CONTROL, sMainState.wDrawBase, KILL_NOTIFY=''
         WIDGET_CONTROL, sMainState.wDrawBase, /DESTROY
      ENDIF
   end
   for i=0, N_ELEMENTS(*(sMainState.hDisplayList)) -1 do begin
      hGraphic = (*(sMainState.hDisplayList))[i]
      IF (TAG_NAMES(*hGraphic, /STRUCTURE_NAME ) $
            EQ 'VIZ3D_SURF') THEN BEGIN
         PTR_FREE, (*hGraphic).hVertList
         PTR_FREE, (*hGraphic).hPolyList
         PTR_FREE, (*hGraphic).hColorList
      ENDIF
      IF (TAG_NAMES( *hGraphic, /STRUCTURE_NAME ) EQ $
         'VIZ3D_OBLIQ_PLANE') THEN BEGIN
         PTR_FREE, (*hGraphic).hVertPlane
         PTR_FREE, (*hGraphic).hVert3D
      ENDIF
      PTR_FREE, hGraphic
;     PTR_FREE, (*(sMainState.hDisplayList))[i]
   endfor
   PTR_FREE, sMainState.hDisplayList
   FOR i=0, N_ELEMENTS(*(sMainState.sThreshState.hDataHist))-1 do begin
      PTR_FREE, (*(sMainState.sThreshState.hDataHist))[i]
   ENDFOR
   PTR_FREE, sMainState.sThreshState.hDataHist
   FOR i=0, N_ELEMENTS(*(sMainState.hDataList))-1 do begin
      PTR_FREE, (*(sMainState.hDataList))[i]
   ENDFOR
   PTR_FREE, sMainState.hDataList
   PTR_FREE, sMainState.hDataName
   PTR_FREE, sMainState.hFrontDepth
   PTR_FREE, sMainState.hBackDepth
   PTR_FREE, sMainState.hFrontImage
   PTR_FREE, sMainState.hBackImage
   PTR_FREE, sMainState.sThreshState.hMinData
   PTR_FREE, sMainState.sThreshState.hMaxData
   PTR_FREE, sMainState.sThreshState.hTransVal
   PTR_FREE, sMainState.sThreshState.hLowThresh
   PTR_FREE, sMainState.sThreshState.hHighThresh
   FOR i=0, N_ELEMENTS(*(sMainState.sSurfState.hRangeHist))-1 do begin
      PTR_FREE, (*(sMainState.sSurfState.hRangeHist))[i]
   ENDFOR
   PTR_FREE, sMainState.sSurfState.hRangeHist
   PTR_FREE, sMainState.sSurfState.hSurfThresh
   PTR_FREE, sMainState.hMain
   WDELETE, sMainState.pixWin
   WDELETE, sMainState.sProfState.profPix
   WDELETE, sMainState.sSurfState.surfPix
   WDELETE, sMainState.sThreshState.threshPix
end
;******************************************************************************


;******************************************************************************
; Main procedure (entry point).

pro Slicer3, hData3D, DATA_NAMES=dataNames, DETACH=detachView, $
             GROUP=groupLead, MODAL=runModal

   ; Use the "Colors" common block to save the current colortable.
   common colors, r, g, b, cur_red, cur_green, cur_blue
   IF (N_ELEMENTS(r) LE 0L) THEN LOADCT, 0
   saveR = r
   saveG = g
   saveB = b
   save_curR = cur_red
   save_curG = cur_green
   save_curB = cur_blue

   ; Check the incoming data (if any).

   maxHand = N_ELEMENTS(hData3D) - 1L
   if (maxHand ge 0L) then begin
      minData = DBLARR(N_ELEMENTS(hData3D))
      maxData = minData

      for i=0, maxHand do begin

         if NOT(PTR_VALID(hData3D(i))) then begin
            PRINT, "Slicer3: data pointer is not valid."
            return
         endif

         neData = N_ELEMENTS((*(hData3D(i))))
         if (neData le 0L) then begin
            PRINT, "Slicer3: data pointer does not contain valid data."
            return
         endif
         szData = SIZE((*(hData3D(i))))
         maxData(i) = MAX((*(hData3D(i))), MIN=minVal)
         minData(i) = minVal
         if (minData(i) eq maxData(i)) then maxData(i) = maxData(i) + 1.0D

         if (szData(0) ne 3L) then begin
            PRINT, "Slicer3: data does not have 3 dimensions."
            return
         endif else begin
            if (min(szData(1:3)) le 1L) then begin
               PRINT, $
                  "Slicer3: all 3 dimensions of data must be greater than 1."
               return
            endif
         endelse

         if (i eq 0) then sz1 = szData else begin
            diffDim = WHERE(szData(0:3) ne sz1(0:3))
            if (diffDim(0) ne (-1)) then begin
               PRINT, $
                  "Slicer3: all data sets must have the same dimensions."
               return
            endif
         endelse

      endfor
      defaultData = 0
   endif else begin
      data3D = BYTARR(2, 2, 2)
      szData = SIZE(data3D)
      minData = 0.0D
      maxData = 1.0D
      defaultData = 1
      dataNames = 'No Data'
      maxHand = 0
   endelse

   ; Main pointers.
   hMain = PTR_NEW(/ALLOCATE_HEAP)
   if (defaultData) then $
      hData3D = PTR_NEW(data3D, /NO_COPY)
   hDataList = PTR_NEW(hData3D)
   hDisplayList = PTR_NEW(/ALLOCATE_HEAP)

   ; Compute and store the histograms for each set of data.
   lDataHist = PTRARR(N_ELEMENTS(hData3D))
   lRangeHist = PTRARR(N_ELEMENTS(hData3D))
   for i=0, maxHand do begin
      data3D = *(hData3D(i))
      bSize = (maxData(i) - minData(i)) / 200.0
      n_samples = N_ELEMENTS((*(hData3D(i))))
      IF (n_samples GT 32767L) THEN BEGIN
         n_samples = 32767L + ROUND(SQRT(n_samples - 32767L))
         s_index = RANDOMU(s, 1)
         s_index = RANDOMU(s, n_samples)
         s_index = LONG(TEMPORARY(s_index) * FLOAT(n_samples))
         hist_data = FLOAT(((*(hData3D(i))))[TEMPORARY(s_index)])
         dataHist = $
            HISTOGRAM(TEMPORARY(hist_data), MIN=minData(i), MAX=maxData(i), $
               BINSIZE=bSize) > 1.0
      ENDIF ELSE BEGIN
         dataHist = $
            HISTOGRAM(FLOAT((*(hData3D(i)))), MIN=minData(i), MAX=maxData(i), $
               BINSIZE=bSize) > 1.0
      ENDELSE

      ; Histogram will be displayed as a LOG plot.
      dataHist = ALOG(FLOAT(TEMPORARY(dataHist))) > 0.0

      lDataHist[i] = PTR_NEW(dataHist)
      lRangeHist[i] = PTR_NEW(dataHist, /NO_COPY)
   endfor
   hDataHist = PTR_NEW(lDataHist, /NO_COPY)
   hRangeHist = PTR_NEW(lRangeHist, /NO_COPY)

   ; Store the name for each set of data.
   dName = STRARR(N_ELEMENTS(hData3D))
   for i=0, maxHand do $
      if (i lt N_ELEMENTS(dataNames)) then $
         dName(i) = STRTRIM(STRING(dataNames(i)), 2) $
      else $
         dName(i) = 'Data ' + STRTRIM(STRING(i), 2)
   hDataName = PTR_NEW(dName)


   ; Save the current system variables.
   saveOrder = !Order
   saveP = !P
   saveX = !X
   saveY = !Y
   saveZ = !Z
   saveWin = !D.Window

   ; Reset system variables to startup state.
   Viz3D_Reset

   ; Initialize parameters.

   DEVICE, GET_SCREEN_SIZE=screenSize
   winX = ((screenSize(0) / 2) < 640) > 320 ; View window size.
   winY = ((screenSize(1) / 2) < 640) > 320 ; View window size.
   winX = 8 * ((winX < winY) / 8)
   winY = winX

   detachBase = KEYWORD_SET(detachView)

   ; Create the widgets.

   if (N_ELEMENTS(groupLead) LE 0L) then groupLead = 0L
   if (WIDGET_INFO(groupLead, /VALID_ID)) then begin
      wMainBase = WIDGET_BASE(TITLE='IDL 3D Data Visualizer (Slicer3)', $
         APP_MBAR=wBarBase, /ROW, SPACE=1, XPAD=1, YPAD=1, $
         GROUP_LEADER=groupLead)
   endif else begin
      wMainBase = WIDGET_BASE(TITLE='IDL 3D Data Visualizer (Slicer3)', $
         APP_MBAR=wBarBase, /ROW, SPACE=1, XPAD=1, YPAD=1)
   endelse

   WIDGET_CONTROL, wMainBase, /MANAGED
   wFileMenu = WIDGET_BUTTON(wBarBase, VALUE='File', /MENU)
   wToolMenu = WIDGET_BUTTON(wBarBase, VALUE='Tools', /MENU)
   wHelpMenu = WIDGET_BUTTON(wBarBase, VALUE='About', /MENU, /HELP)
   wHelpBttn = WIDGET_BUTTON(wHelpMenu, VALUE='About Slicer', $
                             UVALUE='wHelpBttn')

   wLoadBttn = WIDGET_BUTTON(wFileMenu, VALUE='Load', UVALUE='wLoadBttn')

   wSaveMenu = WIDGET_BUTTON(wFileMenu, VALUE='Save', /MENU, /SEPARATOR)
   wSaveSubsetBttn = WIDGET_BUTTON(wSaveMenu, VALUE='Save Subset', $
                        UVALUE='wSaveSubsetBttn')
   WIDGET_CONTROL, wSaveSubsetBttn, SENSITIVE=0
   wSaveTiffBttn = WIDGET_BUTTON(wSaveMenu, VALUE='Save Tiff Image', $
                      UVALUE='wSaveTiffBttn')
   wSavePSBttn = WIDGET_BUTTON(wSaveMenu, VALUE='Save PS Image', $
                      UVALUE='wSavePSBttn')


   wQuitBttn = WIDGET_BUTTON(wFileMenu, VALUE='Quit', UVALUE='wQuitBttn', $
                             /SEPARATOR)

   wEraseBttn = WIDGET_BUTTON(wToolMenu, VALUE='Erase', $
                               UVALUE='wEraseBttn')
   wDeleteMenu = WIDGET_BUTTON(wToolMenu, VALUE='Delete', /MENU)

   wColorsMenu = WIDGET_BUTTON(wToolMenu, VALUE='Colors', /MENU)
   wColResetBttn = WIDGET_BUTTON(wColorsMenu, VALUE='Reset Colors', $
                     UVALUE='wColResetBttn')
   wColDiffBttn = WIDGET_BUTTON(wColorsMenu, VALUE='Differential Shading', $
                     UVALUE='wColDiffBttn')
   wColSliceBttn = WIDGET_BUTTON(wColorsMenu, VALUE='Slice/Block', $
                      UVALUE='wColSliceBttn')
   wColSurfBttn = WIDGET_BUTTON(wColorsMenu, VALUE='Surface', $
                     UVALUE='wColSurfBttn')
   wColProjBttn = WIDGET_BUTTON(wColorsMenu, VALUE='Projection', $
                     UVALUE='wColProjBttn')

;   *** Future animation functionality.
;   wAnimateBttn = WIDGET_BUTTON(wToolMenu, VALUE='Animate', $
;                     UVALUE='wAnimateBttn')

   wOptionsBttn = WIDGET_BUTTON(wToolMenu, VALUE='Options', $
                     UVALUE='wOptionsBttn')

   wModeBase = WIDGET_BASE(wMainBase, /COLUMN, SPACE=1, XPAD=1, YPAD=1)
   if ((detachBase) and NOT(KEYWORD_SET(runModal))) then $
      wDrawBase = WIDGET_BASE(TITLE='IDL 3-D Data Visualizer', $
         UVALUE=wMainBase, /COLUMN, SPACE=1, XPAD=1, YPAD=1, $
         EVENT_PRO='Viz3D_SliceEvent') $
   else $
      wDrawBase = WIDGET_BASE(wMainBase, UVALUE=wMainBase, /COLUMN, $
         SPACE=1, XPAD=1, YPAD=1, EVENT_PRO='Viz3D_SliceEvent')

   wStatText = WIDGET_TEXT(wDrawBase, SCR_XSIZE=winX-12, YSIZE=1, VALUE=' ')
   wMainDraw = WIDGET_DRAW(wDrawBase, XSIZE=winX, YSIZE=winY, $
      UVALUE='wMainDraw', RETAIN=2, /BUTTON_EVENTS)

   wDataModeBase = WIDGET_BASE(wModeBase, /ROW, /FRAME, $
                               SPACE=1, XPAD=1, YPAD=1)
   wDataLab = WIDGET_LABEL(wDataModeBase, VALUE='Data:')
   wDataDrop = WIDGET_DROPLIST(wDataModeBase, VALUE=dName, UVALUE='wDataDrop')
   if (N_ELEMENTS(dName) le 1L) then $
      WIDGET_CONTROL, wDataModeBase, SENSITIVE=0

   modes = ['Slice','Block','Surface','Projection', $
            'Threshold','Profile','Probe','View']
   wActModeBase = WIDGET_BASE(wModeBase, /ROW, /FRAME, $
                              SPACE=1, XPAD=1, YPAD=1)
   wDataLab = WIDGET_LABEL(wActModeBase, VALUE='Mode:')
   wModeDrop = WIDGET_DROPLIST(wActModeBase, VALUE=modes, UVALUE='wModeDrop')

   wMapBase = WIDGET_BASE(wModeBase, /FRAME, SPACE=1, XPAD=1, YPAD=1)

   wSliceBase = WIDGET_BASE(wMapBase, /COLUMN, SPACE=1, XPAD=1, YPAD=1, $
      EVENT_PRO='Viz3D_SliceEvent')
   wBlockBase = WIDGET_BASE(wMapBase, /COLUMN, SPACE=1, XPAD=1, YPAD=1, $
      MAP=0, EVENT_PRO='Viz3D_BlockEvent')
   wSurfBase = WIDGET_BASE(wMapBase, /COLUMN, SPACE=1, XPAD=1, YPAD=1, $
      MAP=0, EVENT_PRO='Viz3D_SurfEvent')
   wThreshBase = WIDGET_BASE(wMapBase, /COLUMN, SPACE=1, XPAD=1, YPAD=1, $
      MAP=0, EVENT_PRO='Viz3D_ThreshEvent')
   wProjectBase = WIDGET_BASE(wMapBase, /COLUMN, SPACE=1, XPAD=1, YPAD=1, $
      MAP=0, EVENT_PRO='Viz3D_ProjectEvent')
   wProfileBase = WIDGET_BASE(wMapBase, /COLUMN, SPACE=1, XPAD=1, YPAD=1, $
      MAP=0, EVENT_PRO='Viz3D_ProfileEvent')
   wProbeBase = WIDGET_BASE(wMapBase, /COLUMN, SPACE=1, XPAD=1, YPAD=1, $
      MAP=0, EVENT_PRO='Viz3D_ProbeEvent')
   wViewBase = WIDGET_BASE(wMapBase, /COLUMN, SPACE=1, XPAD=1, YPAD=1, $
      MAP=0, EVENT_PRO='Viz3D_ViewEvent')


   ; Slice base widgets.
   wSliceDraw = WIDGET_DRAW(wSliceBase, XSIZE=128, YSIZE=128, $
      RETAIN=2, UVALUE='wSliceDraw', /BUTTON_EVENTS)

   wBase = WIDGET_BASE(wSliceBase, /ROW, SPACE=1, XPAD=1, YPAD=1, $
              /EXCLUSIVE, /FRAME)
   wSliceDrawBttn = WIDGET_BUTTON(wBase, VALUE='Draw', $
      UVALUE='wSliceDrawBttn', /NO_RELEASE)
   wSliceExpoBttn = WIDGET_BUTTON(wBase, VALUE='Expose', $
      UVALUE='wSliceExpoBttn', /NO_RELEASE)
   WIDGET_CONTROL, wSliceDrawBttn, /SET_BUTTON

   wBase = WIDGET_BASE(wSliceBase, /COLUMN, SPACE=1, XPAD=1, YPAD=1, $
              /EXCLUSIVE, /FRAME)
   wSliceOrthoBttn = WIDGET_BUTTON(wBase, VALUE='Orthogonal', $
      UVALUE='wSliceOrthoBttn', /NO_RELEASE)
   wSliceObliqBttn = WIDGET_BUTTON(wBase, VALUE='Oblique', $
      UVALUE='wSliceObliqBttn', /NO_RELEASE)
   WIDGET_CONTROL, wSliceOrthoBttn, /SET_BUTTON

   wSliceMapBase = WIDGET_BASE(wSliceBase, SPACE=1, XPAD=1, YPAD=1, /FRAME)
   wSliceOrthoBase = WIDGET_BASE(wSliceMapBase, /ROW, SPACE=1, $
                        XPAD=1, YPAD=1, /EXCLUSIVE)
   wSliceXBttn = WIDGET_BUTTON(wSliceOrthoBase, VALUE='X', $
      UVALUE='wSliceXBttn', /NO_RELEASE)
   wSliceYBttn = WIDGET_BUTTON(wSliceOrthoBase, VALUE='Y', $
      UVALUE='wSliceYBttn', /NO_RELEASE)
   wSliceZBttn = WIDGET_BUTTON(wSliceOrthoBase, VALUE='Z', $
      UVALUE='wSliceZBttn', /NO_RELEASE)
   WIDGET_CONTROL, wSliceXBttn, /SET_BUTTON

   wSliceObliqBase = WIDGET_BASE(wSliceMapBase, /Column, SPACE=1, $
                        XPAD=1, YPAD=1, MAP=0)
   wSliceObliqBase1 = WIDGET_BASE(wSliceObliqBase, /ROW, SPACE=1, $
                        XPAD=1, YPAD=1, /EXCLUSIVE)
   wSliceNormalBttn = WIDGET_BUTTON(wSliceObliqBase1, VALUE='Normal', $
      UVALUE='wSliceNormalBttn', /NO_RELEASE)
   wSliceCenterBttn = WIDGET_BUTTON(wSliceObliqBase1, VALUE='Center', $
      UVALUE='wSliceCenterBttn', /NO_RELEASE)
   WIDGET_CONTROL, wSliceNormalBttn, /SET_BUTTON
   wSliceObliqGoBttn = WIDGET_BUTTON(wSliceObliqBase, VALUE='Display', $
                                     UVALUE='wSliceObliqGoBttn')


   ; Block base widgets.
   wBlockDraw = WIDGET_DRAW(wBlockBase, XSIZE=128, YSIZE=128, $
      RETAIN=2, UVALUE='wBlockDraw', /BUTTON_EVENTS)
   wBlockModeBase = WIDGET_BASE(wBlockBase, /ROW, /EXCLUSIVE, /FRAME, $
                       XPAD=1, YPAD=1, SPACE=1)
   wBlockAddBttn = WIDGET_BUTTON(wBlockModeBase, VALUE='Add', $
                      UVALUE='wBlockAddBttn', /NO_RELEASE)
   wBlockSubBttn = WIDGET_BUTTON(wBlockModeBase, VALUE='Subtract', $
                      UVALUE='wBlockSubBttn', /NO_RELEASE)
   WIDGET_CONTROL, wBlockAddBttn, /SET_BUTTON
   wBlockGoBttn = WIDGET_BUTTON(wBlockBase, VALUE='Display', $
                                UVALUE='wBlockGoBttn')


   ; Iso-surface base widgets.
   wSurfDrawBase = WIDGET_BASE(wSurfBase, /COLUMN, SPACE=1, $
                               XPAD=1, YPAD=1, /FRAME)
   wSurfLabel = WIDGET_LABEL(wSurfDrawBase, VALUE='Surface Threshold', $
                             /ALIGN_LEFT)
   wSurfDraw = WIDGET_DRAW(wSurfDrawBase, XSIZE=192, YSIZE=128, $
      RETAIN=2, UVALUE='wSurfDraw', /BUTTON_EVENTS)
   wSurfText = WIDGET_TEXT(wSurfDrawBase, XSIZE=8, YSIZE=1, /FRAME, $
                           UVALUE='wSurfText', /EDITABLE)
   wSurfSideBase = WIDGET_BASE(wSurfBase, /ROW, SPACE=1, /EXCLUSIVE, $
                               XPAD=1, YPAD=1, /FRAME)
   wSurfLowBttn = WIDGET_BUTTON(wSurfSideBase, VALUE='Low', $
                                UVALUE='wSurfLowBttn', /NO_RELEASE)
   wSurfHighBttn = WIDGET_BUTTON(wSurfSideBase, VALUE='High', $
                                 UVALUE='wSurfHighBttn', /NO_RELEASE)
   WIDGET_CONTROL, wSurfLowBttn, /SET_BUTTON
   wSurfShadeBase = WIDGET_BASE(wSurfBase, /COLUMN, $
                       SPACE=1, XPAD=1, YPAD=1, /FRAME)
   wSurfShadeLab = WIDGET_LABEL(wSurfShadeBase, VALUE='Shading:')
   sName = dName
   sName(0) = 'Light-source'
   wSurfShadeDrop = WIDGET_DROPLIST(wSurfShadeBase, VALUE=sName, $
      UVALUE='wSurfShadeDrop')
   if (N_ELEMENTS(sName) LE 1L) then WIDGET_CONTROL, wSurfShadeBase, $
      SENSITIVE=0
   wSurfGoBttn = WIDGET_BUTTON(wSurfBase, VALUE='Display', $
                                     UVALUE='wSurfGoBttn')


   ; Projection base widgets.
   wProjTypeBase = WIDGET_BASE(wProjectBase, /COLUMN, SPACE=1, $
                               XPAD=1, YPAD=1, /FRAME)
   wProjTypeLab = WIDGET_LABEL(wProjTypeBase, VALUE='Projection Type', $
                     /ALIGN_LEFT)
   wProjTypeBase1 = WIDGET_BASE(wProjTypeBase, /ROW, /EXCLUSIVE, SPACE=1, $
                               XPAD=1, YPAD=1)
   wProjMaxBttn = WIDGET_BUTTON(wProjTypeBase1, VALUE='Max', $
                                UVALUE='wProjMaxBttn', /NO_RELEASE)
   wProjAvgBttn = WIDGET_BUTTON(wProjTypeBase1, VALUE='Avg', $
                                UVALUE='wProjAvgBttn', /NO_RELEASE)
   WIDGET_CONTROL, wProjMaxBttn, /SET_BUTTON
   wProjResoBase = WIDGET_BASE(wProjectBase, /COLUMN, SPACE=1, $
                               XPAD=1, YPAD=1, /FRAME)
   wProjResoLab = WIDGET_LABEL(wProjResoBase, VALUE='Resolution', $
                     /ALIGN_LEFT)
   wProjResoBase1 = WIDGET_BASE(wProjResoBase, /ROW, /EXCLUSIVE, $
                                SPACE=1, XPAD=1, YPAD=1)
   wProjLowBttn = WIDGET_BUTTON(wProjResoBase1, VALUE='Low', $
                                UVALUE='wProjLowBttn', /NO_RELEASE)
   wProjMedBttn = WIDGET_BUTTON(wProjResoBase1, VALUE='Med', $
                                UVALUE='wProjMedBttn', /NO_RELEASE)
   wProjHighBttn = WIDGET_BUTTON(wProjResoBase1, VALUE='High', $
                                 UVALUE='wProjHighBttn', /NO_RELEASE)
   WIDGET_CONTROL, wProjMedBttn, /SET_BUTTON
   wProjDQSlid = WIDGET_SLIDER(wProjectBase, MIN=0, MAX=100, VALUE=100, $
                             UVALUE='wProjDQSlid', TITLE='Depth Queue %')
   wProjGoBttn = WIDGET_BUTTON(wProjectBase, VALUE='Display', $
                               UVALUE='wProjGoBttn')


   ; Profile base widgets.
   profWinY = 256
   wProfDraw = WIDGET_DRAW(wProfileBase, XSIZE=192, YSIZE=profWinY, $
      RETAIN=2, UVALUE='wProfDraw', /BUTTON_EVENTS)
   wProfTypeBase = WIDGET_BASE(wProfileBase, /COLUMN, /EXCLUSIVE, $
                               XPAD=1, YPAD=1, SPACE=1, /FRAME)
   wProfOrthoBttn = WIDGET_BUTTON(wProfTypeBase, Value='Orthogonal', $
      UVALUE='wProfOrthoBttn', /ALIGN_LEFT)
   wProfObliqBttn = WIDGET_BUTTON(wProfTypeBase, Value='Oblique', $
      UVALUE='wProfObliqBttn', /ALIGN_LEFT)
   WIDGET_CONTROL, wProfOrthoBttn, /SET_BUTTON


   ; Probe base widgets.
   wProbeXBase = WIDGET_BASE(wProbeBase, /ROW, XPAD=1, YPAD=1, SPACE=1)
   wProbeXLab = WIDGET_LABEL(wProbeXBase, VALUE='X')
   wProbeXText = WIDGET_TEXT(wProbeXBase, UVALUE='wProbeText', /EDITABLE, $
      VALUE=STRTRIM(STRING(FLOAT(szData(1)-1)/2.0)), XSIZE=12)
   wProbeYBase = WIDGET_BASE(wProbeBase, /ROW, XPAD=1, YPAD=1, SPACE=1)
   wProbeYLab = WIDGET_LABEL(wProbeYBase, VALUE='Y')
   wProbeYText = WIDGET_TEXT(wProbeYBase, UVALUE='wProbeText', /EDITABLE, $
      VALUE=STRTRIM(STRING(FLOAT(szData(2)-1)/2.0)), XSIZE=12)
   wProbeZBase = WIDGET_BASE(wProbeBase, /ROW, XPAD=1, YPAD=1, SPACE=1)
   wProbeZLab = WIDGET_LABEL(wProbeZBase, VALUE='Z')
   wProbeZText = WIDGET_TEXT(wProbeZBase, UVALUE='wProbeText', /EDITABLE, $
      VALUE=STRTRIM(STRING(FLOAT(szData(3)-1)/2.0)), XSIZE=12)


   ; View base widgets.
   ang1 = ( 30)
   ang2 = (-60)
   ang3 = (  0)
   wViewDrawBase = WIDGET_BASE(wViewBase, /ROW, SPACE=1, XPAD=1, YPAD=1)
   wViewDraw = WIDGET_DRAW(wViewDrawBase, XSIZE=128, YSIZE=128, $
      RETAIN=2, UVALUE='wViewDraw', /BUTTON_EVENTS)
   wViewDispBttn = WIDGET_BUTTON(wViewDrawBase, VALUE='Display', $
                      UVALUE='wViewDispBttn', /ALIGN_BOTTOM)

   wRotBase = WIDGET_BASE(wViewBase, /COLUMN, /FRAME, $
                 SPACE=1, XPAD=1, YPAD=1)
   wRotLab = WIDGET_LABEL(wRotBase, VALUE='Rotations:')

   wRot1Base = WIDGET_BASE(wRotBase, /ROW, SPACE=1, XPAD=1, YPAD=1)
   wRot1Slid = WIDGET_SLIDER(wRot1Base, Min=(-180), Max=(180), $
      VALUE=ang1, UVALUE='wRot1Slid', XSIZE=136, /DRAG)
   wRot1Drop = WIDGET_DROPLIST(wRot1Base, VALUE=['X','Y','Z'], $
                  UVALUE='wRot1Drop')
   WIDGET_CONTROL, wRot1Drop, SET_DROPLIST_SELECT=2

   wRot2Base = WIDGET_BASE(wRotBase, /ROW, SPACE=1, XPAD=1, YPAD=1)
   wRot2Slid = WIDGET_SLIDER(wRot2Base, Min=(-180), Max=(180), $
      VALUE=ang2, UVALUE='wRot2Slid', XSIZE=136, /DRAG)
   wRot2Drop = WIDGET_DROPLIST(wRot2Base, VALUE=['X','Y','Z'], $
                  UVALUE='wRot2Drop')
   WIDGET_CONTROL, wRot2Drop, SET_DROPLIST_SELECT=0

;   wRot3Base = WIDGET_BASE(wRotBase, /ROW, SPACE=1, XPAD=1, YPAD=1)
;   wid = WIDGET_LABEL(wRot3Base, VALUE='3')
;   wRot3Slid = WIDGET_SLIDER(wRot3Base, Min=(-180), Max=(180), $
;      VALUE=ang3, UVALUE='wRot3Slid', XSIZE=136, /DRAG)
;   wRot3Drop = WIDGET_DROPLIST(wRot3Base, VALUE=['X','Y','Z'], $
;                  UVALUE='wRot3Drop')
;   WIDGET_CONTROL, wRot3Drop, SET_DROPLIST_SELECT=0

   wZoomBase = WIDGET_BASE(wViewBase, /ROW, SPACE=1, XPAD=1, YPAD=1, $
                           /FRAME)
   wZoomSlid = WIDGET_SLIDER(wZoomBase, Min=50, Max=100, VALUE=50, $
                  XSIZE=136, UVALUE='wZoomSlid', /DRAG)
   wid = WIDGET_LABEL(wZoomBase, VALUE='Zoom %')

   wZFacBase = WIDGET_BASE(wViewBase, /ROW, SPACE=1, XPAD=1, YPAD=1, $
                           /FRAME)
   wZFacSlid = WIDGET_SLIDER(wZFacBase, Min=50, Max=200, VALUE=100, $
                  XSIZE=136, UVALUE='wZFacSlid', /DRAG)
   wid = WIDGET_LABEL(wZfacBase, VALUE='Z %')

;   wPerspBase = WIDGET_BASE(wViewBase, /ROW, SPACE=1, XPAD=1, YPAD=1, $
;                            /FRAME)
;   wPerspSlid = WIDGET_SLIDER(wPerspBase, Min=0, Max=10, VALUE=0, $
;                   XSIZE=136, UVALUE='wPerspSlid', /DRAG)
;   wid = WIDGET_LABEL(wPerspBase, VALUE='Persp.')


   ; Threshold base widgets.
   wThreshLabel = WIDGET_LABEL(wThreshBase, VALUE='Data Threshold', $
                               /ALIGN_LEFT)
   wThreshDraw = WIDGET_DRAW(wThreshBase, XSIZE=192, YSIZE=128, $
      RETAIN=2, UVALUE='wThreshDraw', /BUTTON_EVENTS)
   wThreshLowBase = WIDGET_BASE(wThreshBase, /ROW, SPACE=1, $
                                XPAD=1, YPAD=1)
   wThreshLowText = WIDGET_TEXT(wThreshLowBase, XSIZE=8, YSIZE=1, /FRAME, $
                                UVALUE='wThreshLowText', /EDITABLE)
   wThreshLowLabel = WIDGET_LABEL(wThreshLowBase, VALUE='Min')
   wThreshHighBase = WIDGET_BASE(wThreshBase, /ROW, SPACE=1, $
                                 XPAD=1, YPAD=1)
   wThreshHighText = WIDGET_TEXT(wThreshHighBase, XSIZE=8, YSIZE=1, /FRAME, $
                                UVALUE='wThreshHighText', /EDITABLE)
   wThreshHighLabel = WIDGET_LABEL(wThreshHighBase, VALUE='Max')
   wThreshTransBase = WIDGET_BASE(wThreshBase, /ROW, SPACE=1, $
                                  XPAD=1, YPAD=1)
   wThreshTransText = WIDGET_TEXT(wThreshTransBase, XSIZE=8, YSIZE=1, $
                         /FRAME, UVALUE='wThreshTransText', /EDITABLE)
   wThreshTransLabel = WIDGET_LABEL(wThreshTransBase, VALUE='Transp.')

   ; Set the color mode.
   if (!D.N_Colors GT 256L) then DEVICE, DECOMPOSED=0

   ; Create the pixmap used for profiling.
   WINDOW, /FREE, /PIXMAP, XSIZE=192, YSIZE=profWinY, RETAIN=2
   profPix = !D.Window

   ; Create the pixmap used for the surface histogram.
   WINDOW, /FREE, /PIXMAP, XSIZE=192, YSIZE=128, RETAIN=2
   surfPix = !D.Window

   ; Create the pixmap used for the threshold histogram.
   WINDOW, /FREE, /PIXMAP, XSIZE=192, YSIZE=128, RETAIN=2
   threshPix = !D.Window

   ; Create the pixmap used for dynamic screen updating.
   WINDOW, /FREE, /PIXMAP, XSIZE=winX, YSIZE=winY, RETAIN=2
   XYOUTS, 0, 0, '.', /DEVICE
   ERASE
   pixWin = !D.Window

   ; Initialize the Z buffer.
   screenDevice = !D.Name
   SET_PLOT, 'Z'
   DEVICE, /Z_BUFFERING, SET_RESOLUTION=[winX, winY]
   SET_PLOT, screenDevice

   if ((detachBase) and NOT(KEYWORD_SET(runModal))) then $
      WIDGET_CONTROL, wDrawBase, /REALIZE
   WIDGET_CONTROL, wMainBase, /REALIZE

   WIDGET_CONTROL, /HOURGLASS

   WIDGET_CONTROL, wMainDraw, GET_VALUE=mainWin
   WIDGET_CONTROL, wSliceDraw, GET_VALUE=sliceWin
   WIDGET_CONTROL, wBlockDraw, GET_VALUE=blockWin
   WIDGET_CONTROL, wSurfDraw, GET_VALUE=SurfWin
   WIDGET_CONTROL, wThreshDraw, GET_VALUE=threshWin
   WIDGET_CONTROL, wProfDraw, GET_VALUE=profWin
   WIDGET_CONTROL, wViewDraw, GET_VALUE=viewWin
   WSET, mainWin

   ; Set up the color state.
   sColorState = Viz3D_LoadColor(DIFFSHADE=20)

   ; View state.
   sViewState = Viz3D_View(viewWin, ANG1=ang1, ANG2=ang2, ANG3=ang3, $
      DIR1='Z', DIR2='X', DIR3='Y', ZOOM=0.75, XMAX=szData(1)-1L, $
      YMAX=szData(2)-1L, ZMAX=szData(3)-1L)
   WIDGET_CONTROL, wZoomSlid, SET_VALUE=ROUND(100.0 * sViewState.zoomFac)

   ; Threshold state.
   hMinData = PTR_NEW(minData)
   hMaxData = PTR_NEW(maxData)
   hLowThresh = PTR_NEW(minData)
   hHighThresh = PTR_NEW(maxData)
   hTransVal = PTR_NEW(minData)
   sThreshState = {SThreshState, $
                   histColor:7, lowColor:2, highColor:5, $
                   axisColor:3, transColor:6, $
                   wThreshDraw:wThreshDraw, $
                   wThreshLowText:wThreshLowText, $
                   wThreshHighText:wThreshHighText, $
                   wThreshTransText:wThreshTransText, $
                   threshWin:threshWin, hDataHist:hDataHist, $
                   threshPix:threshPix, $
                   hMinData:hMinData, hMaxData:hMaxData, $
                   hTransVal:hTransVal, $
                   hLowThresh:hLowThresh, $
                   hHighThresh:hHighThresh, $
                   moveType:0}
   WIDGET_CONTROL, wThreshLowText, $
      SET_VALUE=STRTRIM(STRING(minData(0)),2)
   WIDGET_CONTROL, wThreshHighText, $
      SET_VALUE=STRTRIM(STRING(maxData(0)),2)

   ; Slice state.
   sSliceState = $
      {SSliceState, sliceWin:sliceWin, $
       wSliceOrthoBase:wSliceOrthoBase, wSliceObliqBase:wSliceObliqBase, $
       wSliceXBttn:wSliceXBttn, wSliceYBttn:wSliceYBttn, $
       wSliceZBttn:wSliceZBttn, drawMode:0, planeMode:1, orthoDir:1, $
       orthoPos:0.5, obliqNormal:[0.0,0.0,1.0], obliqCenter:[0.5,0.5,0.5], $
       obliqMove:0, edgeColor:6, fillColor:1, normColor:5}

   ; Block state.
   sBlockState = $
      {SBlockState, blockWin:blockWin, $
       blockMode:1, c1:[szData(1)-1, szData(2)-1, szData(3)-1]/3, $
       c2:(2*[szData(1)-1, szData(2)-1, szData(3)-1])/3, $
       blockColor:6, c1Color:4, c2Color:1, frontBack:1}

   ; Iso-surface state.
   hSurfThresh = PTR_NEW((TEMPORARY(minData) + TEMPORARY(maxData))/2.0D)
   sSurfState = $
      {SSurfState, surfWin:surfWin, wSurfShadeBase:wSurfShadeBase, $
       surfPix:surfPix, wSurfShadeDrop:wSurfShadeDrop, curShade:0, $
       histColor:7, tColor:2, axisColor:3, hRangeHist:hRangeHist, $
       surfSide:0, hSurfThresh:hSurfThresh, wSurfText:wSurfText}

   ; Projection state.
   sProjState = $
      {SProjState, projType:0, projReso:1, depthQ:100}

   ; Profile state.
   sProfState = $
      {SProfState, profType:0, profDir:1, frontBack:1, onFace:1, $
       x1Ortho:FLOAT(szData(1)-1)/2.0, x2Ortho:FLOAT(szData(1)-1)/2.0, $
       y1Ortho:FLOAT(szData(2)-1)/2.0, y2Ortho:FLOAT(szData(2)-1)/2.0, $
       z1Ortho:0.0, z2Ortho:FLOAT(szData(3)-1), $
       x1Obliq:FLOAT(szData(1)-1)/2.0, x2Obliq:FLOAT(szData(1)-1)/2.0, $
       y1Obliq:FLOAT(szData(2)-1)/2.0, y2Obliq:FLOAT(szData(2)-1)/2.0, $
       z1Obliq:0.0, z2Obliq:FLOAT(szData(3)-1), $
       lineColor:6, axisColor:3, markColor1:4, markColor2:1, $
       profWin:profWin, profPix:profPix, profWinY:profWinY}

   ; Profile state.
   sProbeState = $
      {SProbeState, wProbeXText:wProbeXText, wProbeYText:wProbeYText, $
                    wProbeZText:wProbeZText, x:FLOAT(szData(1)-1)/2.0, $
                    y:FLOAT(szData(2)-1)/2.0, z:FLOAT(szData(3)-1)/2.0, $
                    lineColor:1}

   ; Main state.
   sMainState = {SVizState, curData:0, $
       szData:szData, detachBase:detachBase, winX:winX, winY:winY, $
       actMode:0, wMainBase:wMainBase, wDrawBase:wDrawBase, $
       wMainDraw:wMainDraw, wSliceBase:wSliceBase, wBlockBase:wBlockBase, $
       wSurfBase:wSurfBase, wThreshBase:wThreshBase, $
       wProjectBase:wProjectBase, wProfileBase:wProfileBase, $
       wProbeBase:wProbeBase, wViewBase:wViewBase, wCurMapBase:wSliceBase, $
       wStatText:wStatText, wSaveSubsetBttn:wSaveSubsetBttn, $
       wDeleteMenu:wDeleteMenu, wDataModeBase:wDataModeBase, $
       wDataDrop:wDataDrop, mainWin:mainWin, pixWin:pixWin, interpMode:1, $
       hMain:hMain, hDisplayList:hDisplayList, hDataList:hDataList, $
       hDataName:hDataName, cubeOn:1, axisOn:1, $
       cubeColor:3, axisColor:2, backColor:0, screenDevice:screenDevice, $
       cleanupView:0, cleanupBuffer:0, $
       hFrontDepth:PTR_NEW(Intarr(winX,winY)), $
       hBackDepth:PTR_NEW(Intarr(winX,winY)), $
       hFrontImage:PTR_NEW(Bytarr(winX,winY)), $
       hBackImage:PTR_NEW(Bytarr(winX,winY)), $
       sColorState:Temporary(sColorState), $
       sSliceState:Temporary(sSliceState), $
       sBlockState:Temporary(sBlockState), $
       sSurfState:Temporary(sSurfState), $
       sProjState:Temporary(sProjState), $
       sViewState:Temporary(sViewState), $
       sProfState:Temporary(sProfState), $
       sProbeState:Temporary(sProbeState), $
       sThreshState:Temporary(sThreshState) $
    }

   ; Fill the depth and image buffers.
   Viz3D_FillDepth, sMainState

   ; Draw the contents of the small windows.
   Viz3D_SliceShow, sMainState
   Viz3D_SurfShow, sMainState, /CALC_HIST
   Viz3D_ThreshShow, sMainState
   Viz3D_ViewShow, sMainState

   ; Draw cube in main window.
   SET_PLOT, 'Z'
   ERASE
   SET_PLOT, screenDevice
   Viz3D_DrawData, sMainState

   ; Put the main state in the top level base.
   WIDGET_CONTROL, wMainBase, SET_UVALUE=sMainState, /NO_COPY

   ; Turn off the hourglass cursor.
   event = WIDGET_EVENT(wMainBase, /NOWAIT)

   if (KEYWORD_SET(runModal)) then begin

      ; Start event processing.

      XMANAGER, 'Slicer3', wMainBase, EVENT_HANDLER='Viz3D_Event', $
         /MODAL, CLEANUP='Viz3D_KillMain'


      ; All done, cleanup.


      ; Restore the system variables to their previous state.
      !Order = saveOrder
      !P = TEMPORARY(saveP)
      !X = TEMPORARY(saveX)
      !Y = TEMPORARY(saveY)
      !Z = TEMPORARY(saveZ)

      ; Restore the previous color table.
      r = TEMPORARY(saveR)
      g = TEMPORARY(saveG)
      b = TEMPORARY(saveB)
      cur_red = TEMPORARY(save_curR)
      cur_green = TEMPORARY(save_curG)
      cur_blue = TEMPORARY(save_curB)

      ; Make the previous IDL window active.
      if (saveWin ge 0L) then WSET, saveWin

   endif else begin

      ; Start event processing.
      if (detachBase) then $
         XMANAGER, 'Slicer3', wDrawBase, /JUST_REG, $
            EVENT_HANDLER='Viz3D_SliceEvent', CLEANUP='Viz3D_KillDraw'

      XMANAGER, 'Slicer3', wMainBase, EVENT_HANDLER='Viz3D_Event', $
         CLEANUP='Viz3D_KillMain'

   endelse

end
;******************************************************************************

