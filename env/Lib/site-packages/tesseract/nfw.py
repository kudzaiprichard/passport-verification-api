#!/usr/bin/python
"""
nfw
===

This module provides methods for determining various aspects relating to NFW
profiles.

Attributes:
    _default_delta (float): Virial overdensity factor that should be used any 
        time it is not provided. It is specified by the `default-rhoc` option
        in the `nfw-options` section of the tesseract config file.
    _default_rhoc (float): Critical density that should be used any time it is 
        not provided. it is specified by the `default-rhoc` option in the 
        `nfw-options` section of the tesseract config file. Some values of note
        include:

            *1.1845e2 ``Msol/kpc**3`` - Value assumed by langmm's buildgal.
            *1.4139e2 ``Msol/kpc**3`` - Value assumed by khb's buildgal.
            *1.4624e2 ``Msol/kpc**3`` - Value derived from the Gadget2 cosmology.

"""

# Standard packages
import numpy as np
import copy

# Sister modules
from . import util
from . import config_parser

# ------------------------------------------------------------------------------
# DEFAULT VALUES FOR COSMOLOGY
_default_delta = config_parser.getfloat('nfw-options','default-delta')
_default_rhoc = config_parser.getfloat('nfw-options','default-rhoc')

# ------------------------------------------------------------------------------
# UTILITIES RELATING TO DETERMINING VIRIAL QUANTITIES
def calc_virial(r,rho=None,menc=None,delta=None,rhoc=None):
    """Calculates the virial radius and mass of a halo from particle 
    information. Either rho or menc must be provided. If both are provided, 
    rho is used. If the whole halo is more dense than the virial overdensity, 
    the maximum radius is used.
    
    Args:
        r (np.ndarray): (N,) Particle radii.
        rho (Optional[np.ndarray]): (N,) Mean enclosed density at each r.
        menc (Optional[np.ndarray]): (N,) Enclosed mass at each r.
        delta (Optional[float]): Factor determining virial overdensity 
            (default = `_default_delta`)
        rhoc (Optional[float]): Critical density in the same units as r and 
            menc or rho (default = `_default_rhoc`)

    Returns:
        rvir (float): Virial radius.
        mvir (float): Virial mass.

    Raises:
        ValueError: If neither rho or menc is provided.
        Exception: If no part of the halo is not denser than the virial 
            overdensity.
        

    .. todo:: Raise exception for halo that is entirely denser than viral?

    """
    # Set defaults
    if delta is None: delta = _default_delta
    if rhoc is None: rhoc = _default_rhoc
    delrho = delta*rhoc
    # Calculated rho if not provided
    if rho is None:
        if menc is None:
            raise ValueError('Both rho and menc are None. Provided at least one.')
        rho = menc/util.sphvol(r)
    # Find where density exceeds virial overdensity
    idxvir = (rho>=delrho)
    if not np.any(idxvir):
        raise Exception('Density does not exceed virial overdensity. delta*rhoc={}, min(rho)={}, max(rho)={}'.format(delrho,min(rho),max(rho)))
    if np.all(idxvir):
        pass
        #raise Exception('All densities exceed virial overdensity. delta*rhoc={}, min(rho)={}, max(rho)={}'.format(delrho,min(rho),max(rho)))
    # Find this radius and mass starting from radius
    rvir = np.max(r[idxvir])
    mvir = delrho*util.sphvol(rvir)
    # Return
    return rvir,mvir
def calc_rvir(*args,**kwargs):
    """Calculates virial radius. See `nfw.calc_virial` for details.

    Returns:
        rvir (float): Virial radius.

    """
    return get_virial(*args,**kwargs)[0]
def calc_mvir(*args,**kwargs):
    """Calculates virial mass. See `nfw.calc_virial` for details.

    Returns:
        mvir (float): Virial mass.

    """
    return get_virial(*args,**kwargs)[1]

def calc_rhalf(r,menc,mvir=None,delta=None,rhoc=None):
    """Calculate half mass radius from particle information.
    
    Args:
        r (np.ndarray): (N,) Particle radii.
        menc (np.ndarray): (N,) Enclosed mass at each radius.
        mvir (Optional[float]): Virial mass (provide to save time).
        delta (Optional[float]): Factor determining virial overdensity.
            (default = `_default_delta`)
        rhoc (Optional[float]): Critical density in the same units as r and 
            menc. (default = `_default_delta`)

    Returns:
        float: Radius enclosing half the Virial mass.

    Raises:
        Exception: If all masses are greater than the half mass.
        Exception: If all masses are less than the half mass.

    """
    # Get enclosed mass and virial mass
    if mvir is None:
        mvir = calc_mvir(r,menc=menc,delta=delta,rhoc=rhoc)
    # Half mass
    idx_over = (menc>=(mvir/2.))
    if np.all(idx_over):
        raise Exception('Masses are all greater than the half mass. '+
                        'Mhalf={}, min(Menc)={}, max(Menc)={}'.format(mvir/2.,min(menc),max(menc)))
    if not np.any(idx_over):
        raise Exception('Masses are all smaller than the half mass. '+
                        'Mhalf={}, min(Menc)={}, max(Menc)={}'.format(mvir/2.,min(menc),max(menc)))
    rhalf_min = np.max(r[np.logical_not(idx_over)])
    rhalf_max = np.min(r[idx_over])
    rhalf = (rhalf_min+rhalf_max)/2.
    return rhalf
    
# ------------------------------------------------------------------------------
# FITTING AND SUCH
def _nfwgc(c):
    """Utility function for NFW g(r_rs).
    
    Args:
        c (float): Concentraton.

    Returns:
        float: np.log(1+c) - c/(1+c)

    """
    return (np.log(1.+c)-c/(1.+c))

def profile(r,ms=1.,rs=1.,G=1.,method='mass'):
    """Return NFW profile at provided radii for the given parameters.

    Args:
        r (np.ndarray): (N,) Radii.
        ms (Optional[float]): Scale mass. (default = 1)
        rs (Optional[float]): Scale radius. (default = 1)
        G (Optional[float]): Gravitational constant in correct units. 
            Only used for `method = 'pot'`. (default = 1) 
        method (Optional[str]): Type of profile that should be returned. 
            (default = 'mass') Valid values include:
                'rho' : density profile
                'mass': cummulative mass profile
                'pot' : potential profile
    
    Returns:
        np.ndarray: (N,) Density, enclosed mass, or potential profile at
            each radius. 

    Raises:
        ValueError: If `method` is not a recognized option.

    """
    if ms is None: ms = 1.
    if rs is None: rs = 1.
    if G is None: G = 1.
    # Set dimensionless parameter
    rho0=ms/(4.*np.pi*(rs**3.))
    s=r/rs
    # Calculate profile
    if   method=='rho' : out=rho0/(s*((1.0+s)**2.0))
    elif method=='mass': out=ms*(np.log(1.+s)-s/(1.+s))
    elif method=='pot' : out=-(G*ms/rs)*(np.log(1.+s)/s)
    else: raise ValueError('Invalid profile method: {}'.format(method))
    # Return
    return out

def fit(r,y,rs=1.,ms=1.,G=1.,method='mass',
        plotflag=False,plotfile=None,**kwargs):
    """Fit an NFW profile to provided data.

    Args:
        r (np.ndarray): (N,) Radii.
        y (np.ndarray): (N,) Variable to be fit (depends on method).
        rs (Optional[float]): Initial guess at scale radius. (default = 1)
        ms (Optional[float]): Initial guess at scale mass. (default = 1)
        G (Optional[float]): Gravitational constant in correct units. Only 
            required for `method = 'pot'`. (default = 1)
        method (Optional[str]): Type of profile that should be fit to y. 
            (default = 'mass') Valid values include:
                'rho' : density profile
                'mass': cummulative mass profile
                'pot' : potential profile
        plotflag (Optional[bool]): When True, the fit is plotted. 
            (default = False)
        plotfile (Optional[str]): Path to file where plot should be saved.
        **kwargs: Additional keywords are passed to the minimize/leastsq 
            function from the scipy.optimize module, depending on if minimize 
            exists or not.

    Returns:
        ms (float): Scale mass.
        rs (float): Scale radius.

    Raises:
        Exception: If fit is not successful.
        
    """
    # Set initial guess and bounds
    p0 = [ms,rs]
    bounds = [(0.,None),(0.,None)]
    # Create fit functions
    fitfunc = lambda p,ri: profile(ri,ms=p[0],rs=p[1],method='mass')
    errfunc = lambda p: fitfunc(p,r)-y
    # Fit
    p1,success = util.myleastsq(errfunc,p0,bounds=bounds,**kwargs)
    if success not in [1,2,3,4]:
        raise Exception('Fit was unsuccessful with code {}'.format(success))
    # Plot
    if plotflag:
        plotfit(r,y,ms=p1[0],rs=p1[1],G=G,method=method,plotfile=plotfile)
    # Return
    return p1

def plotfit(r,ydat,rs=None,ms=None,G=1.,method='mass',
            plotfile=None,residuals=True,label='fit',color='k',
            axs=None,axs_res=None,**kwargs):
    """Plot a fit to an NFW profile.

    Args:
        r (np.ndarray): (N,) Radii.
        ydat (np.ndarray): (N,) Measured values at r. (depends on method)
        rs (Optional[float]): Scale radius. If not provided, it is fit.
        ms (Optional[float]): Scale mass. If not provided, it is fit.
        G (Optional[float]): Gravitational constant in the correct units. Only used if
            `method = 'pot'`. (default = 1)
        method (Optional[str]): Description of what profile ydat contains.
            (default = 'mass') Valid values include:
                'rho' : density profile
                'mass': cummulative mass profile
                'pot' : potential profile
        plotfile (Optional[str]): Full path to file where plot should be saved.
            If not provided, the plot is displayed instead.
        residuals (Optional[bool]): If True, fit residuals are plotted in 
            axs_res. (default = False)
        label (Optional[str]): Label to give the profile lines.  
            (default = 'fit')
        color (Optional[tuple,str]): Color for profile lines.
        axs (Optional[plt.Axes]): Existing axes that profile should be plotted 
            in. If not provided, one is created.
        axs_res (Optional[plt.Axes]): Existing axes that profile should be 
            plotted in. If not provided, one is created.
        **kwargs: Additional keyword arguments are passed to `nfw.fit`.

    """
    import matplotlib.pyplot as plt
    labelx = -0.1
    # Determine which points you should plot
    Npmax = 99000 # maximum number of plot points that should be used
    Np = len(r)
    if (Np/Npmax)>1:
        idxplt = np.arange(0,Np,Np/Npmax)
    else:
        idxplt = slice(0,Np)
    # Fit parameters if not provided & get profile
    if rs is None or ms is None:
        ms,rs = fit(r,ydat,G=G,method=method,**kwargs)
    yfit = profile(r,ms=ms,rs=rs,G=G,method=method)
    yres = (ydat-yfit)/ydat
    # Get bounds
    xlim = (r.min(),r.max())
    ylim = (ydat.min(),ydat.max())
    rmax = np.max(np.abs(yres))
    rlim = (-rmax,rmax)
    # Set up axes
    if axs is None:
        plt.clf()
        fig = plt.figure(figsize=(10,5))
        skipout = False
    else:
        fig = axs.get_figure()
        skipout = True
    if residuals:
        if axs is None: 
            axs = fig.add_axes((.1,.3,.6,.6))
            plt.setp( axs.get_xticklabels(), visible=False)
            axs_res = fig.add_axes((.1,.1,.6,.2))
        else:
            xlim0 = axs.get_xlim()
            ylim0 = axs.get_ylim()
            rlim0 = axs_res.get_ylim()
            xlim = (min(xlim[0],xlim0[0]),max(xlim[1],xlim0[1]))
            ylim = (min(ylim[0],ylim0[0]),max(ylim[1],ylim0[1]))
            rlim = (min(rlim[0],rlim0[0]),max(rlim[1],rlim0[1]))
        axs_res.set_xlabel('Radius')
        axs_res.set_ylabel('Residuals')
        axs.set_ylabel(method.title())
    else:
        if axs is None:
            axs = plt.subplot(1,1,1)
        axs.set_xlabel('Radius')
        axs.set_ylabel(method.title())
    # Plot profile
    axs.loglog(r[idxplt],ydat[idxplt],c=color,ls='--',label=label+' data')
    axs.loglog(r[idxplt],yfit[idxplt],c=color,ls='-',label=label)
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    # Plot residuals
    if residuals:
        axs_res.semilogx(r[idxplt],(ydat-yfit)[idxplt]/ydat[idxplt],c=color,ls='None',marker='.')
        axs_res.set_xlim(xlim)
        axs_res.set_ylim(rlim)
        # axs_res.yaxis.set_label_coords(labelx,0.5)
        # axs.yaxis.set_label_coords(labelx,0.5)
    # Do things and save/show
    if not skipout:
        # Legend
        leg = axs.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        # Remove bottom tick
        if residuals:
            fig.canvas.draw()
            tics = [item.get_text() for item in axs.get_yticklabels()]
            tics[0] = ''
            axs.set_yticklabels(tics)
        # Save/show
        if plotfile is None:
            plt.show()
        else:
            plt.savefig(plotfile,bbox_inches='tight',bbox_extra_artists=[leg])
            print '    '+plotfile

def calc_nfw(r,m=None,menc=None,method='rhalf',sortby=None,issorted=False,
             hist=False,nbins=100,rmin=None,rmax=None,
             rvir=None,mvir=None,vvir=None,vpeak=None,rvpeak=None,rhalf=None,
             delta=None,rhoc=None,G=None,plotflag=False,plotfile=None,
             label=None,axs=None,axs_res=None,residuals=True,color='k'):
    """Calculates NFW parameters based on particle information. Either m or 
    menc must be provided. If both are provided, menc is used. 

    Args:
        r (np.ndarray): (N,) Particle radii.
        m (Optional[np.ndarray]): (N,) Particle masses. Must be provided if 
            hist is True.
        menc (Optional[np.ndarray]): (N,) Mass enclosed at each radius.
        method (Optional[str]): Method that should be used to determine NFW 
            parameters. (default = 'halfmass') Valid values include:
                'rhalf': half-mass radius relation (See concen_rhalf)
                'vpeak': peak circular velocity (See concen_vpeak)
                'fit'  : leastsq fitting to mass profile
        sortby (Optional[np.ndarray]): (N,) array that masses should be 
            sorted by before determining the mass enclosed. If not provided, 
            r is used to sort.
        issorted (Optional[np.ndarray]): If True, arrays are assumed to be 
            sorted in order to save time. (default = False)
        hist (Optional[np.ndarray]): If True, NFW parameters are determine 
            after creating a mass weighted histogram of particle positions. 
            In this case, m must be provided. (default = False) 
        nbins (Optional[int]): Number of bins that should be used for the 
            histogram. (default = 100)
        rmin (Optional[float]): Minimum radius that should be used in 
            calculations. (default = 0.05*rvir)
        rmax (Optional[float]): Maximum radius that should be used in 
            calculations. (default = rvir)
        rvir (Optional[float]): Virial radius (saves time).
        mvir (Optional[float]): Virial mass (saves time).
        vvir (Optional[float]): Virial velocity (saves time, not much).
        vpeak (Optional[float]): Peak circular velocity (saves time).
        rvpeak (Optional[float]): Radius where vpeak occurs (saves time).
        delta (Optional[float]): Factor determining virial overdensity. 
            (default = `_default_delta`)
        rhoc (Optional[float]): Critical density in the same units as r and 
            menc. (default = `_default_rhoc`)
        G (Optional[float]): Gravitational constant in correct units. 
            (default = 1)

    Returns:
        nfwpar: A dictionary of parameters relation to the determined
            NFW profile. Keys include:
                mvir (float): Virial mass.
                rvir (float): Virial radius.
                vvir (float): Circular velocity at the virial radius.
                rhalf (float): Radius enclosing half the virial mass.
                vpeak (float): Peak circular velocity.
                rvpeak (float): Radius where circular velocity peaks.
                c (float): Concentration.
                ms (float): Scale mass.
                rs (float): Scale radius.
                rho0 (float): Normalizing density. 

    Raises:
        ValueError: If `hist` is True and `m` is not provided.
        ValueError: If provided `method` is unsupported.

    """
    # Set defaults
    if sortby is None and not issorted:
        sortby = r
    if hist and m is None:
        raise ValueError('Using a weighted histogram requires m be provided.')
    # Allow for binning from method
    if method.endswith('_binned'):
        method = copy.deepcopy(method.split('_binned')[0])
        hist = True
    # Virial radius and mass
    if menc is None:
        menc = util.calc_menc(m,sortby=sortby)
    if rvir is None or mvir is None:
        rvir,mvir = calc_virial(r,menc=menc,delta=delta,rhoc=rhoc)
    # Radius bounds
    if rmin is None: rmin = max(0.05*rvir,r.min())
    if rmax is None: rmax = min(rvir,r.max())
    # Get arrays
    if hist:
        rbins = np.logspace(np.log10(rmin),np.log10(rmax),nbins+1)
        mbins,rbins = np.histogram(r,bins=rbins,weights=m)
        rarr = rbins[1:]
        marr = np.cumsum(mbins)+np.sum(m[r<rmin])
        rvir,mvir = calc_virial(rarr,marr,delta=delta,rhoc=rhoc)
    else:
        idx = np.logical_and(r>=rmin,r<=rmax)
        rarr = r[idx]
        marr = menc[idx]
    # Variables that may be used to calculate concentration
    if rhalf is None:
        rhalf = calc_rhalf(rarr,marr,mvir=mvir,delta=delta,rhoc=rhoc)
    if vvir is None:
        vvir = util.calc_vcirc(rvir,mvir,G=G)
    if vpeak is None:
        vcirc = util.calc_vcirc(rarr,marr,G=G)
        if rvpeak is not None:
            idxpeak = np.argmax(np.where(rarr<=rvpeak)[0])
            vpeak = vcirc[idxpeak]
        else:
            idxpeak = np.argmax(vcirc)
            rvpeak = rarr[idxpeak]
            vpeak = vcirc[idxpeak]
    elif rvpeak is None:
        vcirc = util.calc_vcirc(rarr,marr,G=G)
        rvpeak = np.max(rarr[vcirc<=vpeak])
    rs = None ; ms = None
    # Determine concentration
    if method in ['halfmass','rhalf']: 
        c = concen_rhalf(rhalf,rvir,method='leastsq')
    elif method == 'lokas2001': 
        c = concen_rhalf(rhalf,rvir,method='lokas2001')
    elif method in ['vpeak','vmax','prada2012']:
        c = concen_vpeak(vpeak,vvir,method='prada2012')
    elif method in ['fit','fitting','leastsq']:
        ms,rs = fit(rarr,marr,method='mass',ms=mvir/2.,rs=rhalf)
        c = rvir/rs
    else: 
        raise ValueError('Unsupported method: {}'.format(method))
    # Other NFW parameters
    if rs is None: rs = rvir/c
    rho0 = mvir/((4./3.)*np.pi*(rs**3.)*_nfwgc(c))
    if ms is None: ms = rho0*(4./3.)*np.pi*(rs**3.)
    # Plot
    if plotflag:
        if label is None: label=method.title()
        plotfit(rarr,marr,ms=ms,rs=rs,method='mass',axs=axs,axs_res=axs_res,
                residuals=residuals,label=label,color=color)
    # Return dictionary of parameters
    nfwpar = dict(mvir=mvir,rvir=rvir,vvir=vvir,
                  rhalf=rhalf,vpeak=vpeak,rvpeak=rvpeak,
                  c=c,ms=ms,rs=rs,rho0=rho0)
    return nfwpar

def concen_rhalf(rhalf,rvir,method='leastsq'):
    """Determine concentration from the relationship between the half mass and
    virial radii:

        0.5 = [ln(1+s*c)-sc/(1+s*c)]/[ln(1+c)-c/(1+c)], where s=rhalf/rvir

    Args:
        rhalf (float): Radius enclosing half the virial mass.
        rvir (float): Virial radius.
        method (Optional[str]): Method used to solve for concentration. 
            (default = 'leastsq') Supported values include:
                'leastsq'  : leastsq is used to solve the complete relationship.
                'lokas2001': leastsq is used to solve the simplified 
                    relationship from Lokas & Mammon (2001). Eqn. 28 in arxiv 
                    version.

    Returns:
        float: Concentration.

    Raises:
        ValueError: If provided method is not supported.

    """
    from scipy.optimize import leastsq
    if method=='leastsq':
        errfunc = lambda x: 0.5-(_nfwgc(x*rhalf/rvir)/_nfwgc(x))
        c = leastsq(errfunc,x0=3.)[0][0]    
    elif method=='lokas2001':
        errfunc = lambda x: ((rhalf/rvir)-(0.6082-0.1843*np.log10(x)-0.1011*(np.log10(x)**2)+0.03918*(np.log10(x)**3)))
        c = leastsq(errfunc,x0=3.)[0][0]
    else:
        raise ValueError('Unsupported method: {}'.format(method))
    return c

def concen_vpeak(vpeak,vvir,method='prada2012'):
    """Calculate concentration based on relationship between peak and virial
    velocities:

        vpeak/vvir = {0.216*c/[ln(1+c)-c/(1+c)]}**0.5

    Args:
        vpeak (float): Peak circular velocity.
        vvir (float): Circular velocity at virial radius.
        method (Optional[str]): Method used to solve for concentration. 
            (default = 'prada2012') Supported values include:
                'prada2012': leastsq is used to solve the equation from Prada 
                    et al. (2012). Eqn. 9 in arxiv version.

    Returns:
        float: Concentration.

    Raises:
        ValueError: If provided method is not supported.

    """
    from scipy.optimize import leastsq
    if method == 'prada2012':
        errfunc = lambda x: ((vpeak/vvir)-np.sqrt(0.216*x/_nfwgc(x)))
        c = leastsq(errfunc,x0=3.)[0][0]
    else:
        raise ValueError('Unsupported method: {}'.format(method))
    return c
