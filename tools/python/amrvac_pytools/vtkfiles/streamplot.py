"""
Streamline plotting like Mathematica.
Copyright (c) 2011 Tom Flannaghan.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

version = '4'

import numpy
import pylab
import matplotlib
import matplotlib.patches as mpp
import matplotlib.colors
import matplotlib.cm

def streamplot(x, y, u, v, x_0=None, y_0=None, density=1, linewidth=1,
               color='k', cmap=None, norm=None, vmax=None, vmin=None,
               arrowsize=1, INTEGRATOR='RK4',alpha=1.):
    """Draws streamlines of a vector flow.

    * x and y are 1d arrays defining an *evenly spaced* grid.
    * u and v are 2d arrays (shape [y,x]) giving velocities.
    * density controls the closeness of the streamlines. For different
      densities in each direction, use a tuple or list [densityx, densityy].
    * linewidth is either a number (uniform lines) or a 2d array
      (variable linewidth).
    * color is either a color code (of any kind) or a 2d array. This is
      then transformed into color by the cmap, norm, vmin and vmax args.
      A value of None gives the default for each.

    INTEGRATOR is experimental. Currently, RK4 should be used.
      """

    ## Sanity checks.
    assert len(x.shape)==1
    assert len(y.shape)==1
    assert u.shape == (len(y), len(x))
    assert v.shape == (len(y), len(x))
    if type(linewidth) == numpy.ndarray:
        assert linewidth.shape == (len(y), len(x))
    if type(color) == numpy.ndarray:
        assert color.shape == (len(y), len(x))

    ## Set up some constants - size of the grid used.
    NGX = len(x)
    NGY = len(y)
    ## Constants used to convert between grid index coords and user coords.
    DX = x[1]-x[0]
    DY = y[1]-y[0]
    XOFF = x[0]
    YOFF = y[0]

    ## Now rescale velocity onto axes-coordinates
    u = u / (x[-1]-x[0])
    v = v / (y[-1]-y[0])
    speed = numpy.sqrt(u*u+v*v)
    ## s (path length) will now be in axes-coordinates, but we must
    ## rescale u for integrations.
    u *= NGX
    v *= NGY
    ## Now u and v in grid-coordinates.

    ## Blank array: This is the heart of the algorithm. It begins life
    ## zeroed, but is set to one when a streamline passes through each
    ## box. Then streamlines are only allowed to pass through zeroed
    ## boxes. The lower resolution of this grid determines the
    ## approximate spacing between trajectories.
#    if (x_0 == None and y_0 == None):
    if type(density) == float or type(density) == int:
        assert density > 0
        NBX = int(30*density)
        NBY = int(30*density)
    else:
        assert len(density) > 0
        NBX = int(30*density[0])
        NBY = int(30*density[1])
#    else:
#        NBX = NGX
#        NBY = NGY
    blank = numpy.zeros((NBY,NBX))
    
    ## Constants for conversion between grid-index space and
    ## blank-index space
    bx_spacing = NGX/float(NBX-1)
    by_spacing = NGY/float(NBY-1)

    def blank_pos(xi, yi):
        ## Takes grid space coords and returns nearest space in
        ## the blank array.
        return int((xi / bx_spacing) + 0.5), \
               int((yi / by_spacing) + 0.5)

    def value_at(a, xi, yi):
        ## Linear interpolation - nice and quick because we are
        ## working in grid-index coordinates.
        if type(xi) == numpy.ndarray:
            x = xi.astype(numpy.int)
            y = yi.astype(numpy.int)
        else:
            x = numpy.int(xi)
            y = numpy.int(yi)
        a00 = a[y,x]
        a01 = a[y,x+1]
        a10 = a[y+1,x]
        a11 = a[y+1,x+1]
        xt = xi - x
        yt = yi - y
        a0 = a00*(1-xt) + a01*xt
        a1 = a10*(1-xt) + a11*xt
        return a0*(1-yt) + a1*yt

    def rk4_integrate(x0, y0):
        ## This function does RK4 forward and back trajectories from
        ## the initial conditions, with the odd 'blank array'
        ## termination conditions. TODO tidy the integration loops.
        
        def f(xi, yi):
            if (value_at(speed, xi, yi) != 0.):
                dt_ds = 1./value_at(speed, xi, yi)
            else:
                dt_ds = 0.
            ui = value_at(u, xi, yi)
            vi = value_at(v, xi, yi)
            return ui*dt_ds, vi*dt_ds

        def g(xi, yi):
            if (value_at(speed, xi, yi) != 0.):
                dt_ds = 1./value_at(speed, xi, yi)
            else:
                dt_ds = 0.
            ui = value_at(u, xi, yi)
            vi = value_at(v, xi, yi)
            return -ui*dt_ds, -vi*dt_ds

        check = lambda xi, yi: xi>=0 and xi<NGX-1 and yi>=0 and yi<NGY-1 and xi==xi and yi==yi

        bx_changes = []
        by_changes = []

        ## Integrator function
        def rk4(x0, y0, f):
            ds = 0.01 #min(1./NGX, 1./NGY, 0.01)
            stotal = 0
            xi = x0
            yi = y0
            xb, yb = blank_pos(xi, yi)
            xf_traj = []
            yf_traj = []
            while check(xi, yi):
                # Time step. First save the point.
                xf_traj.append(xi)
                yf_traj.append(yi)
                # Next, advance one using RK4
                try:
                    k1x, k1y = f(xi, yi)
                    k2x, k2y = f(xi + .5*ds*k1x, yi + .5*ds*k1y)
                    k3x, k3y = f(xi + .5*ds*k2x, yi + .5*ds*k2y)
                    k4x, k4y = f(xi + ds*k3x, yi + ds*k3y)
                except IndexError:
                    # Out of the domain on one of the intermediate steps
                    break
                xi += ds*(k1x+2*k2x+2*k3x+k4x) / 6.
                yi += ds*(k1y+2*k2y+2*k3y+k4y) / 6.
                # Final position might be out of the domain
                if not check(xi, yi): break
                stotal += ds
                # Next, if s gets to thres, check blank.
                new_xb, new_yb = blank_pos(xi, yi)
                if new_xb != xb or new_yb != yb:
                    # New square, so check and colour. Quit if required.
                    if blank[new_yb,new_xb] == 0:
                        blank[new_yb,new_xb] = 1
                        bx_changes.append(new_xb)
                        by_changes.append(new_yb)
                        xb = new_xb
                        yb = new_yb
                    else:
                        break
                if stotal > 2:
                    break
            return stotal, xf_traj, yf_traj
        
        ## Alternative Integrator function

        ## RK45 does not really help in it's current state. The
        ## resulting trajectories are accurate but low-resolution in
        ## regions of high curvature and thus fairly ugly. Maybe a
        ## curvature based cap on the maximum ds permitted is the way
        ## forward.
    
        def rk45(x0, y0, f):
            maxerror = 0.001
            maxds = 0.03
            ds = 0.03
            stotal = 0
            xi = x0
            yi = y0
            xb, yb = blank_pos(xi, yi)
            xf_traj = []
            yf_traj = []
            while check(xi, yi):
                # Time step. First save the point.
                xf_traj.append(xi)
                yf_traj.append(yi)
                # Next, advance one using RK45
                try:
                    k1x, k1y = f(xi, yi)
                    k2x, k2y = f(xi + .25*ds*k1x,
                                 yi + .25*ds*k1y)
                    k3x, k3y = f(xi + 3./32*ds*k1x + 9./32*ds*k2x,
                                 yi + 3./32*ds*k1y + 9./32*ds*k2y)
                    k4x, k4y = f(xi + 1932./2197*ds*k1x - 7200./2197*ds*k2x + 7296./2197*ds*k3x,
                                 yi + 1932./2197*ds*k1y - 7200./2197*ds*k2y + 7296./2197*ds*k3y)
                    k5x, k5y = f(xi + 439./216*ds*k1x - 8*ds*k2x + 3680./513*ds*k3x - 845./4104*ds*k4x,
                                 yi + 439./216*ds*k1y - 8*ds*k2y + 3680./513*ds*k3y - 845./4104*ds*k4y)
                    k6x, k6y = f(xi - 8./27*ds*k1x + 2*ds*k2x - 3544./2565*ds*k3x + 1859./4104*ds*k4x - 11./40*ds*k5x,
                                 yi - 8./27*ds*k1y + 2*ds*k2y - 3544./2565*ds*k3y + 1859./4104*ds*k4y - 11./40*ds*k5y)
                    
                except IndexError:
                    # Out of the domain on one of the intermediate steps
                    break
                dx4 = ds*(25./216*k1x + 1408./2565*k3x + 2197./4104*k4x - 1./5*k5x)
                dy4 = ds*(25./216*k1y + 1408./2565*k3y + 2197./4104*k4y - 1./5*k5y)
                dx5 = ds*(16./135*k1x + 6656./12825*k3x + 28561./56430*k4x - 9./50*k5x + 2./55*k6x)
                dy5 = ds*(16./135*k1y + 6656./12825*k3y + 28561./56430*k4y - 9./50*k5y + 2./55*k6y)

                ## Error is normalized to the axes coordinates (it's a distance)
                error = numpy.sqrt(((dx5-dx4)/NGX)**2 + ((dy5-dy4)/NGY)**2)
                if error < maxerror:
                    # Step is within tolerance so continue
                    xi += dx5
                    yi += dy5
                    # Final position might be out of the domain
                    if not check(xi, yi): break
                    stotal += ds
                    # Next, if s gets to thres, check blank.
                    new_xb, new_yb = blank_pos(xi, yi)
                    if new_xb != xb or new_yb != yb:
                        # New square, so check and colour. Quit if required.
                        if blank[new_yb,new_xb] == 0:
                            blank[new_yb,new_xb] = 1
                            bx_changes.append(new_xb)
                            by_changes.append(new_yb)
                            xb = new_xb
                            yb = new_yb
                        else:
                            break
                    if stotal > 2:
                        break
                # Modify ds for the next iteration.
                if len(xf_traj) > 2:
                    ## hacky curvature dependance:
                    v1 = numpy.array((xf_traj[-1]-xf_traj[-2], yf_traj[-1]-yf_traj[-2]))
                    v2 = numpy.array((xf_traj[-2]-xf_traj[-3], yf_traj[-2]-yf_traj[-3]))
                    costheta = (v1/numpy.sqrt((v1**2).sum()) * v2/numpy.sqrt((v2**2).sum())).sum()
                    if costheta < .8:
                        ds = .01
                        continue
                ds = min(maxds, 0.85*ds*(maxerror/error)**.2)
            return stotal, xf_traj, yf_traj

        ## Forward and backward trajectories
        if INTEGRATOR == 'RK4':
            integrator = rk4
        elif INTEGRATOR == 'RK45':
            integrator = rk45
        
        sf, xf_traj, yf_traj = integrator(x0, y0, f)
        sb, xb_traj, yb_traj = integrator(x0, y0, g)
        stotal = sf + sb
        x_traj = xb_traj[::-1] + xf_traj[1:]
        y_traj = yb_traj[::-1] + yf_traj[1:]

        ## Tests to check length of traj. Remember, s in units of axes.
        if len(x_traj) < 1: return None
        if stotal > .2:
            initxb, inityb = blank_pos(x0, y0)
            blank[inityb, initxb] = 1
            return x_traj, y_traj
        else:
            for xb, yb in zip(bx_changes, by_changes):
                blank[yb, xb] = 0
            return None

    ## A quick function for integrating trajectories if blank==0.
    trajectories = []
    def traj(xb, yb):
        if xb < 0 or xb >= NBX or yb < 0 or yb >= NBY:
            return
        if blank[yb, xb] == 0:
            t = rk4_integrate(xb*bx_spacing, yb*by_spacing)
            if t is not None:
                trajectories.append(t)

    ## Now we build up the trajectory set. I've found it best to look
    ## for blank==0 along the edges first, and work inwards.
    if (x_0 is None and y_0 is None):
        for indent in range(int((max(NBX,NBY))/2)):
            for xi in range(int(max(NBX,NBY)-2*indent)):
                traj(xi+indent, indent)
                traj(xi+indent, NBY-1-indent)
                traj(indent, xi+indent)
                traj(NBX-1-indent, xi+indent)
    else:
#        for myx in zip(x_0,y_0):
#            traj(myx[0]/DX,myx[1]/DY)
#        print 'starting point:', x_0,y_0, (x_0-XOFF)/DX, (y_0-YOFF)/DY
        traj((x_0-XOFF)/DX/bx_spacing,(y_0-YOFF)/DY/by_spacing)
#        print trajectories

    ## PLOTTING HERE.
    #pylab.pcolormesh(numpy.linspace(x.min(), x.max(), NBX+1),
    #                 numpy.linspace(y.min(), y.max(), NBY+1), blank)
    
    # Load up the defaults - needed to get the color right.
    if type(color) == numpy.ndarray:
        if vmin is None: vmin = color.min()
        if vmax is None: vmax = color.max()
        if norm is None: norm = matplotlib.colors.normalize
        if cmap is None: cmap = matplotlib.cm.get_cmap(
            matplotlib.rcParams['image.cmap'])

    for t in trajectories:
        # Finally apply the rescale to adjust back to user-coords from
        # grid-index coordinates.
        tx = numpy.array(t[0])*DX+XOFF
        ty = numpy.array(t[1])*DY+YOFF

        tgx = numpy.array(t[0])
        tgy = numpy.array(t[1])
        
        points = numpy.array([tx, ty]).T.reshape(-1,1,2)
        segments = numpy.concatenate([points[:-1], points[1:]], axis=1)

        args = {}
        if type(linewidth) == numpy.ndarray:
            args['linewidth'] = value_at(linewidth, tgx, tgy)[:-1]
            arrowlinewidth = args['linewidth'][len(tgx)/2]
        else:
            args['linewidth'] = linewidth
            arrowlinewidth = linewidth
            
        if type(color) == numpy.ndarray:            
            args['color'] = cmap(norm(vmin=vmin,vmax=vmax)
                                 (value_at(color, tgx, tgy)[:-1]))
            arrowcolor = args['color'][len(tgx)/2]
        else:
            args['color'] = color
            arrowcolor = color
        args['alpha'] = alpha
        
        lc = matplotlib.collections.LineCollection\
             (segments, **args)
        pylab.gca().add_collection(lc)
            
        ## Add arrows half way along each trajectory.
        n = int(len(tx)/2)
        p = mpp.FancyArrowPatch((tx[n],ty[n]), (tx[n+1],ty[n+1]),
                                arrowstyle='->', lw=arrowlinewidth,
                                mutation_scale=20*arrowsize, color=arrowcolor, alpha=alpha)
        pylab.gca().add_patch(p)

    pylab.xlim(x.min(), x.max())
    pylab.ylim(y.min(), y.max())    
    return

def test():
    pylab.figure(1)
    x = numpy.linspace(-3,3,100)
    y = numpy.linspace(-3,3,100)
    u = -1-x**2+y[:,numpy.newaxis]
    v = 1+x-y[:,numpy.newaxis]**2
    speed = numpy.sqrt(u*u + v*v)
    pylab.subplot(121)
    streamplot(x, y, u, v, density=1, INTEGRATOR='RK4', color='b')
    pylab.subplot(122)
    streamplot(x, y, u, v, density=(1,1), INTEGRATOR='RK4', color=u,
               linewidth=5*speed/speed.max())
    pylab.show()

if __name__ == '__main__':
    test()
