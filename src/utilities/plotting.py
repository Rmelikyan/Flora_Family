# -*- coding: utf-8 -*-
import math
from itertools import cycle

from src.data_collection.Bennuid import Bennuid
from src.utilities.constants import m2au


def OrbitPlot(particles, figsize=None, slices=0, xlim=None, ylim=None, unitlabel=None, color=False, periastron=False, orbit_type="trail", lw=1., primary=None, Narc=128):
    """
    Convenience function for plotting instantaneous orbits.

    Parameters
    ----------
    figsize         : tuple of float, optional
        Tuple defining the figure size (default: (5,5))
    fancy           : bool (default: False)
        Changes various settings to create a fancy looking plot
    slices          : float, optional
        Default is 0, showing the orbits in the xy plane only. Set to a value between 0 and 1 to create three plots, showing the orbits from different directions. The value corresponds to the size of the additional plots relative to the main plot.
    xlim            : tuple of float, optional           
        Limits for x axes (default: None = automatically determined)
    ylim            : tuple of float, optional           
        Limits for y axes (default: None = automatically determined)
    unitlabel       : str, optional          
        String describing the units, shown on axis labels (default: None)
    color           : bool, str or list, optional            
        By default plots are black and white. If set to True, plots use a color cycle. If a string or list of strings, e.g. ['red', 'cyan'], will cycle between passed colors.
    periastron  : bool, optional            
        Draw a marker at periastron (default: False)
    orbit_type       : str, optional
        This argument determines the type of orbit show. By default it shows the orbit as a trailing and fading line ("trail"). Other object are: "solid", None.
    lw              : float, optional           
        Linewidth used in plots (default: 1.)
    plotparticles   : list, optional
        List of particles to plot. Can be a list of any valid keys for accessing sim.particles, i.e., integer indices or hashes (default: plot all particles)
    primary         : rebound.Particle, optional
        Primary to use for the osculating orbit (default: Jacobi center of mass)
    Narc            : int, optional
        Number of points used in an orbit. Increase this number for highly eccentric orbits. (default: 128)

    Returns
    -------
    fig, ax_main, (ax_top, ax_right)
        The function return the matplotlib figure as well as the axes (three axes if slices>0.)

    Examples
    --------
    The following example illustrates a typical use case.

    >>> sim = rebound.Simulation()
    >>> sim.add(m=1)
    >>> sim.add(a=1)
    >>> fig, ax_main = rebound.OrbitPlot(sim)
    >>> fig.savefig("image.png") # save figure to file
    >>> fig.show() # show figure on screen

    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib import gridspec
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        import numpy as np
    except:
        raise ImportError("Error importing matplotlib and/or numpy. Plotting functions not available. If running from within a jupyter notebook, try calling '%matplotlib inline' beforehand.")
    if unitlabel is not None:
        unitlabel = " " + unitlabel
    else:
        unitlabel = ""
    if figsize is None:
        if slices>0.:
            figsize = (8,8)
        else:
            figsize = (5,5)
    
    fig = plt.figure(figsize=figsize)
    ax_main = plt.subplot(111,aspect="equal")
    ax_main.set_xlabel("x"+unitlabel)
    ax_main.set_ylabel("y"+unitlabel)

    if slices>0.:
        divider = make_axes_locatable(ax_main)
        divider.set_aspect(True)
        ax_top   = divider.append_axes("top",  size="%.2f%%"%(100.*slices), sharex=ax_main)
        ax_top.set_aspect('equal', adjustable='datalim')
        ax_right = divider.append_axes("right", size="%.2f%%"%(100.*slices), sharey=ax_main)
        ax_right.set_aspect('equal', adjustable='datalim')
      
        plt.setp(ax_top.get_xticklabels(), visible=False)
        plt.setp(ax_top.get_xticklines(), visible=False)
        ax_top.set_ylabel("z"+unitlabel)
        
        plt.setp(ax_right.get_yticklabels(), visible=False)
        plt.setp(ax_right.get_yticklines(), visible=False)
        ax_right.set_xlabel("z"+unitlabel)

    OrbitPlotOneSlice(particles, ax_main, Narc = Narc, color = color, periastron = periastron, orbit_type = orbit_type, lw = lw, axes = "xy", primary = primary)

    if slices>0.:
        OrbitPlotOneSlice(particles, ax_right, Narc = Narc, color = color, periastron = periastron, orbit_type = orbit_type, lw = lw, axes = "zy", primary = primary)
        OrbitPlotOneSlice(particles, ax_top, Narc = Narc, color = color, periastron = periastron, orbit_type = orbit_type, lw = lw, axes = "xz", primary = primary)

    if xlim is not None:
        ax_main.set_xlim(xlim)
    if ylim is not None:
        ax_main.set_ylim(ylim)

    if slices>0.:
        return fig, ax_main, ax_top, ax_right
    else:
        return fig, ax_main

def get_color(color):
    """
    Takes a string for a color name defined in matplotlib and returns of a 3-tuple of RGB values.
    Will simply return passed value if it's a tuple of length three.

    Parameters
    ----------
    color   : str
        Name of matplotlib color to calculate RGB values for.
    """

    if isinstance(color, tuple) and len(color) == 3: # already a tuple of RGB values
        return color

    try:
        import matplotlib.colors as mplcolors
    except:
        raise ImportError("Error importing matplotlib. If running from within a jupyter notebook, try calling '%matplotlib inline' beforehand.")
   
    try:
        hexcolor = mplcolors.cnames[color]
    except KeyError:
        raise AttributeError("Color not recognized in matplotlib.")

    hexcolor = hexcolor.lstrip('#')
    lv = len(hexcolor)
    return tuple(int(hexcolor[i:i + lv // 3], 16)/255. for i in range(0, lv, lv // 3)) # tuple of rgb values

def fading_line(x, y, color='black', alpha_initial=1., alpha_final=0., glow=False, **kwargs):
    """
    Returns a matplotlib LineCollection connecting the points in the x and y lists, with a single color and alpha varying from alpha_initial to alpha_final along the line.
    Can pass any kwargs you can pass to LineCollection, like linewidgth.

    Parameters
    ----------
    x       : list or array of floats for the positions on the (plot's) x axis
    y       : list or array of floats for the positions on the (plot's) y axis
    color   : matplotlib color for the line. Can also pass a 3-tuple of RGB values (default: 'black')
    alpha_initial:  Limiting value of alpha to use at the beginning of the arrays.
    alpha_final:    Limiting value of alpha to use at the end of the arrays.
    """
    try:
        from matplotlib.collections import LineCollection
        from matplotlib.colors import LinearSegmentedColormap
        import numpy as np
    except:
        raise ImportError("Error importing matplotlib and/or numpy. Plotting functions not available. If running from within a jupyter notebook, try calling '%matplotlib inline' beforehand.")


    if "lw" not in kwargs:
        kwargs["lw"] = 1
    lw = kwargs["lw"]

    if glow:
        kwargs["lw"] = 1*lw
        fl1 = fading_line(x, y, color, alpha_initial, alpha_final, glow=False, **kwargs)
        kwargs["lw"] = 2*lw
        alpha_initial *= 0.5
        alpha_final *= 0.5
        fl2 = fading_line(x, y, color, alpha_initial, alpha_final, glow=False, **kwargs)
        kwargs["lw"] = 6*lw
        alpha_initial *= 0.5
        alpha_final *= 0.5
        fl3 = fading_line(x, y, color, alpha_initial, alpha_final, glow=False, **kwargs)
        return [fl3,fl2,fl1]

    color = get_color(color)
    cdict = {'red': ((0.,color[0],color[0]),(1.,color[0],color[0])),
             'green': ((0.,color[1],color[1]),(1.,color[1],color[1])),
             'blue': ((0.,color[2],color[2]),(1.,color[2],color[2])),
             'alpha': ((0.,alpha_initial, alpha_initial), (1., alpha_final, alpha_final))}
    
    Npts = len(x)
    if len(y) != Npts:
        raise AttributeError("x and y must have same dimension.")
   
    segments = np.zeros((Npts-1,2,2))
    segments[0][0] = [x[0], y[0]]
    for i in range(1,Npts-1):
        pt = [x[i], y[i]]
        segments[i-1][1] = pt
        segments[i][0] = pt 
    segments[-1][1] = [x[-1], y[-1]]

    individual_cm = LinearSegmentedColormap('indv1', cdict)
    lc = LineCollection(segments, cmap=individual_cm, **kwargs)
    lc.set_array(np.linspace(0.,1.,len(segments)))
    return lc

def OrbitPlotOneSlice(particles, ax, Narc=128, color=False, periastron=False, orbit_type="trail", lw=1., axes="xy", primary = None):
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
    # import random

    #ax.set_aspect("equal")

    if color:
        if color == True:
            colors = [(1.,0.,0.),(0.,0.75,0.75),(0.75,0.,0.75),(0.75, 0.75, 0,),(0., 0., 0.),(0., 0., 1.),(0., 0.5, 0.)]
        if isinstance(color, str):
            colors = [get_color(color)]
        if isinstance(color, list):
            colors = []
            for c in color:
                colors.append(get_color(c))
    else:
        colors = ["black"]
    coloriterator = cycle(colors)
   
    if primary is None:
        prim = Bennuid([0,0,0,0,0,0], 'SI', 'Helio', '1788')
    else:
        prim = primary

    ax.scatter(getattr(prim,axes[0]),getattr(prim,axes[1]), marker="*", s=35*lw, facecolor="black", edgecolor=None, zorder=3)
    
    proj = {}
    for p in particles:

        if p.a is None:
            p.calcElements()

        colori = next(coloriterator)

        ax.scatter(getattr(p,axes[0]), getattr(p,axes[1]), s=25*lw, facecolor="black", edgecolor=None, zorder=3)
       
        if orbit_type is not None:
            alpha_final = 0. if orbit_type=="trail" else 1. # fade to 0 with trail

            hyperbolic = bool(p.a < 0) # Boolean for whether orbit is hyperbolic
            if hyperbolic is False:
                pts = np.array(sample_orbit(p, Npts=Narc+1, useTrueAnomaly=True))
                proj['x'],proj['y'],proj['z'] = [pts[:,i] for i in range(3)]
                lc = fading_line(proj[axes[0]], proj[axes[1]], colori, alpha_final=alpha_final, lw=lw)
                if type(lc) is list:
                    for l in lc:
                        ax.add_collection(l)
                else:
                    ax.add_collection(lc)

            else:
                pts = np.array(sample_orbit(p, Npts=Narc+1, useTrueAnomaly=False))
                # true anomaly stays close to limiting value and switches quickly at pericenter for hyperbolic orbit, so use mean anomaly
                proj['x'],proj['y'],proj['z'] = [pts[:,i] for i in range(3)]
                lc = fading_line(proj[axes[0]], proj[axes[1]], colori, alpha_final=alpha_final, lw=lw)
                if type(lc) is list:
                    for l in lc:
                        ax.add_collection(l)
                else:
                    ax.add_collection(lc)
          
                alpha = 0.2 if orbit_type=="trail" else 1.
                pts = np.array(sample_orbit(p, Npts=Narc+1, useTrueAnomaly=False))
                proj['x'],proj['y'],proj['z'] = [pts[:,i] for i in range(3)]
                lc = fading_line(proj[axes[0]], proj[axes[1]], colori, alpha_initial=alpha, alpha_final=alpha, lw=lw)
                if type(lc) is list:
                    for l in lc:
                        ax.add_collection(l)
                else:
                    ax.add_collection(lc)

        if periastron:
            newp = KeplerToCartesian(a=p.a, f=0., inc=p.inc, omega=p.omega, Omega=p.Omega, e=p.e)
            newp = {'x': newp[0], 'y':newp[1], 'z':newp[2]}
            ax.plot([getattr(prim,axes[0]), getattr(newp,axes[0])], [getattr(prim,axes[1]), getattr(newp,axes[1])], linestyle="dotted", c=colori, zorder=1, lw=lw)
            ax.scatter([getattr(newp,axes[0])],[getattr(newp,axes[1])], marker="o", s=5.*lw, facecolor="none", edgecolor=colori, zorder=1)


def sample_orbit(particle, Npts=100, trailing=True, timespan=None, useTrueAnomaly=None):
    """
    Returns a nested list of xyz positions along the osculating orbit of the particle.

    Parameters
    ----------
    Npts    : int, optional  
        Number of points along the orbit to return  (default: 100)
    trailing: bool, optional
        Whether to return points stepping backwards in time (True) or forwards (False). (default: True)
    timespan: float, optional    
        Return points (for the osculating orbit) from the current position to timespan (forwards or backwards in time depending on trailing keyword). 
        Defaults to the orbital period for bound orbits, and to the rough time it takes the orbit to move by the current distance from the primary for a hyperbolic orbit. Implementation currently only supports this option if useTrueAnomaly=False.
    useTrueAnomaly: bool, optional
        Will sample equally spaced points in true anomaly if True, otherwise in mean anomaly.
        Latter might be better for hyperbolic orbits, where true anomaly can stay near the limiting value for a long time, and then switch abruptly at pericenter. (Default: True for bound orbits, False for unbound orbits)
    """
    pts = []
    if particle.a is None:
        particle.calcElements()  # Calc elements for bennuid particle

    if timespan is None:
        if particle.a < 0.: # hyperbolic orbit
            raise Exception("Not Implemented")
            # timespan = 2*math.pi*particle.d/particle.v # rough time to cross display box
        else:
            timespan = particle.P
    
    lim_phase = abs(particle.n)*timespan # n is negative for hyperbolic orbits

    if trailing is True:
        lim_phase *= -1 # sample phase backwards from current value
    phase = [lim_phase*i/(Npts-1) for i in range(Npts)]

    if useTrueAnomaly is None:
        if particle.a <0.:
            raise Exception("Not Implemented")
            # useTrueAnomaly = False
        else:
            useTrueAnomaly = True

    for ph in phase:
        newp = KeplerToCartesian(a=particle.a, f=particle.f+ph, inc=particle.inc, omega=particle.omega, Omega=particle.Omega, e=particle.e)
        pts.append(newp)
    
    return pts

def KeplerToCartesian(a, f, inc, omega, Omega, e):
    ''' From Keplarian elements return cartesian points [x,y,z]'''
    r = a*(1-e**2)/(1+e*math.cos(f))
    v = omega + f
    x = r*(math.cos(Omega)*math.cos(v)- math.sin(Omega)*math.sin(v)*math.cos(inc))
    y = r*(math.sin(Omega)*math.cos(v) + math.cos(Omega)*math.sin(v)*math.cos(inc))
    z = r*(math.sin(inc)*math.sin(v))
    return [x, y, z]

def optimizeGif(path):
    print("TODO: optimizeGif")
    return 0
    #from pygifsicle import optimize
    #import moviepy.editor as mp
    #optimize(path)
    #clip = mp.VideoFileClip(path)
    #clip.write_videofile(path[:-4] + ".mp4")
