to start in a terminal:
$idl $AMRVAC_DIR/tools/.idlrc

.r getpict  ; get files
.r plotfunc ; plot functions
.r animate  ; get files, plot images and make a animation
.r getav    ; get function value or mean value at a space point or spot (functional in 1D; under developing)
.r getmax   ; get function maximum in a given domain (functional in 1D; under developing)

doanimate=0 ; set value 0 to not creat gui animation window for saving system memory when drawing large amount of pictures(e.g., more than 200 pictures) 
iextrapol=1 ; set value 1 to extrapolate amr data to one big finest grid; 0 do not
bottomline=1 ; set value 0 to display no bottom line of "it time= "; 1 to display only "time= "
headerline=0 ; 1 to display headerline; 0 do not display headerline
savemovie='y' ; 'y' to save plotted images; 'n' do not save
autorange='n' & noautorange=1 ; 'n' and 1 to manually set ranges of functions to plot; 
                              ;  'y' and '0' to automatically set range to be maximum and minimum of the functions.
domainset=0 ; 0 to set space domain manually
vertical=1 ; 1 to plot functions in a vertical sequence; 0 to plot functions in a matrix distribution
islope=1   ; use linear restruction to extraplote
filename=' ' ; set file to be load in
func='T rho pth' ; set functions to be T rho pth
plotmode='plotn_amr' ; set plotmodes
xtitle1d='string' ; set title for x axis in 1D.

Example 1, single snapshot:
$idl $AMRVAC_DIR/tools/.idlrc
.....
IDL>.r getpict
% Compiled module: $MAIN$.
filename(s)   ? datamr/yourfile0010.out
filetype(s)   = bin_amr
npictinfile(s)=      -1
npict? 1
.....
extrapolate to single grid (0/1/2)       ? 0
IDL>.r plotfunc
% Compiled module: $MAIN$.
physics (e.g. mhd12)      = hd11 
======= CURRENT PLOTTING PARAMETERS ================
ax,az=  30, 30, contourlevel=100, velvector= 200, velspeed (0..5)= 5
multiplot= 0 (default), axistype (coord/cells)=coord
bottomline=2, headerline=2
======= PLOTTING PARAMETERS =========================
wnames                     = rho m1 e Te cond rhog H dradcool radcool criteria
func(s) (e.g. rho p m1;m2 r*m1 -T) ? 

.....

Example 2, series of snapshots:
$doidlcat datamr/yourfile 0 10 1
.....
 all will be gathered in  datamr/yourfileall.out
$idl $AMRVAC_DIR/tools/.idlrc
.....
IDL>.r getpict
% Compiled module: $MAIN$.
filename(s)   ? datamr/yourfileall.out
filetype(s)   = bin_amr
npictinfile(s)=      -1
npict? 1
.....
extrapolate to single grid (0/1/2)       ? 0

IDL> .r animate
% Compiled module: $MAIN$.
======= CURRENT ANIMATION PARAMETERS ================
firstpict=   1, dpict=   1, npictmax=3000, savemovie (y/n)=n
ax,az=  30, 30, contourlevel=100, velvector= 200, velspeed (0..5)= 5
multiplot= 0 (default), axistype (coord/cells)=coord
bottomline=2, headerline=2
======= FILE DESCRIPTION ============================
filename(s)   = datamr/relaxaall.out
filetype(s)   = bin_amr
npictinfile(s)=      -1
headline                  =prominence
variables                 =X rho m1 e Te cond rhog H dradcool radcool criteria gamma mu g1 p1 p2 p3 p4 p5 p  (ndim= 1, nw=10)
physics (e.g. mhd12)      = hd11
vertical       0
======= PLOTTING PARAMETERS =========================
No transformations allowed for AMR files
setting transform to amr
extrapolate to single grid (0/1/2)       ? 0
NOAUTODOMAIN    INT       =        0
func(s) (e.g. rho p m1;m2 r*m1 -T) ? Te logn v1 logpth
