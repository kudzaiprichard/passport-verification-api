#!/usr/bin/python
"""
tests
=====

This module provides methods for running package tests. These can be broken
into single tests and test series. 

Attributes:

    _outputdir (str): Directory where test output will be saved.
    _installdir (str): Directory containing the distribution installation. 
    _halofmt (int): File format of the test halo snapshots.
    _halodir (str): Directory containing test halo snapshots.
    _copydir (str): Directory containing other realizations of test halo 
        snapshots.
    _list_conc (List[float]): Concentrations of available test halos.
    _list_series (List[str]): Available test series.
    _default_series (str): Test series that ``run_series()`` selects by
        default.
    _default_conc (float): Halo concentration that ``run_test()`` and 
        ``run_series()`` use by default.
    _default_subm (float): Substructure mass that ``run_test()`` and 
        ``run_series()`` use by default for tests with substructure.
    _default_subr (float): Substructure position that ``run_test()`` and 
        ``run_series()`` use by default for tests with substructure.
    _default_subc (float): Substructure concentration that ``run_test()`` and 
        ``run_series()`` use by default for tests with substructure.
    _default_subrho (float): Substructure density that ``run_test()`` and 
        ``run_series()`` use by default for tests with substructure.
    _default_ellip (float): Halo ellipticity that ``run_series()`` uses
        by default for 'triax' series.
    _nfwmeth (List[str]): Methods that are used to measure the concentration 
        for tests.

"""

# Standard packages
import numpy as np
import os,copy

# Sister modules
from . import util
from . import io
from . import voro
from . import config_parser
from . import config

# ------------------------------------------------------------------------------
# DIRECTORIES
# ===========
_outputdir = config['outputdir']
_installdir = os.path.dirname(os.path.realpath(__file__))
_halofmt = config_parser.getint('test-options','snapshot-format')
if config_parser.has_option('test-options','halodir'):
    _halodir = config_parser.get('test-options','halodir').strip()
else:
    _halodir = os.path.join(_installdir,'halos')
if config_parser.has_option('test-options','copydir'):
    _copydir = config_parser.get('test-options','copydir').strip()
else:
    _copydir = _halodir

# ------------------------------------------------------------------------------
# DEFAULT VALUES FOR TESTS
# ========================
_list_conc = map(float,
                 config_parser.get('test-options','avail-conc').split(','))
_list_series = map(str.strip,
                   config_parser.get('test-options','avail-series').split(','))
_default_series = config_parser.get('test-options','default-series').strip()
_default_conc = config_parser.getfloat('test-options','default-conc')
_default_subm = config_parser.getfloat('test-options','default-subm')
_default_subr = config_parser.getfloat('test-options','default-subr')
_default_subc = config_parser.getfloat('test-options','default-subc')
_default_subrho = config_parser.getfloat('test-options','default-subrho')
_default_ellip = config_parser.getfloat('test-options','default-ellip')
_nfwmeth = map(str.strip,
               config_parser.get('test-options','nfw-methods').split(','))

# ------------------------------------------------------------------------------
# METHODS FOR NAMING TESTS
# ========================
def testid(prefix='',c=_default_conc,
           squishy=False,squishz=False,decimate=False,substr=False,
           subm=_default_subm,subr=_default_subr,subc=_default_subc,
           subrho=_default_subrho,version=-1):
    """Creates a standardized string that uniquely identifies a test.
    
    Args:
        prefix (Optional[str]): string to start ID string with (default = '')
        c (Optional[float]): concentration of run (default = `_default_conc`)
        squishy (Optional[float]): factor that y coords should be squished by 
            (default = False)
        squishz (Optional[float]): factor that z coords should be sauished by 
            (default = False)
        decimate (Optional[int]): factor that particle number should be decimated 
            by (default = False)
        substr (Optional[bool]): If true, the test includes substructure 
            (default = False)
        subm (Optional[float]): mass of subhalo relative to parent's Mvir. 
            Only used if `substr` == True. (default = `_default_subm`)
        subr (Optional[float]): radius of subhalo relative to parent's Rvir. 
            Only used if `substr` == True. (default = `_default_subr`)
        subc (Optional[float]): concentration of substructure. 
            Only used if `substr` == True. (default = `_default_subc`)
        subrho (Optional[float]): central density of subhalo relative to parent's central
            density. Only used if `substr` == True. (default = `_default_subrho`)
        version (Optional[int]): test realization. Set to -1 for the fiducial 
            version. (default = -1)

    Returns:
        str: A unique string identifying a test.

    Raises:
        ValueError: If the generated ID string is empty.

    """
    idstr = prefix
    if c is not None:
        idstr+= '_c{}'.format(util.num2str(c))
    # General tags
    if decimate is not None and decimate!=False and decimate > 1: 
        idstr+= '_dec{}'.format(util.num2str(decimate))
    if squishy is not None and squishy!=False and squishy > 0 and squishy!=1:
        idstr+='_sqY{}'.format(util.num2str(float(squishy)))
    if squishz is not None and squishz!=False and squishz > 0 and squishz!=1:
        idstr+='_sqZ{}'.format(util.num2str(float(squishz)))
    # Substructure tags
    if substr:
        idstr+=subid(subm=subm,subr=subr,subc=subc,subrho=subrho)
    # Version number
    if version>=0:
        idstr+='_ver{:03d}'.format(version)
    # Return
    if len(idstr)==0:
        raise ValueError('Empty ID string.')
    elif idstr.startswith('_'):
        return idstr[1:]
    else:
        return idstr

def subid(subm=_default_subm,subr=_default_subr,subc=_default_subc,
          subrho=_default_subrho):
    """Creates a unique ID string for a set of substructure parameters.

    Args:
        subm (Optional[float]): mass of subhalo relative to parent's Mvir. 
            Only used if `substr` == True. (default = _default_subm)
        subr (Optional[float]): radius of subhalo relative to parent's Rvir. 
            Only used if `substr` == True. (default = _default_subr)
        subc (Optional[float]): concentration of substructure. 
            Only used if `substr` == True. (default = _default_subc)
        subrho (Optional[float]): central density of subhalo relative to parent's central
            density. Only used if `substr` == True. (default = _default_subrho)
    """
    idstr=''
    if subm is not None: idstr+='_subm{}'.format(util.num2str(subm))
    if subr is not None: idstr+='_subr{}'.format(util.num2str(subr))
    if subc is not None: idstr+='_subc{}'.format(util.num2str(subc))
    if subrho is not None: idstr+='_subrho{}'.format(util.num2str(subrho))
    return idstr

def _count_copies(c):
    """Determine the number of realizations that exist of a halo with a given 
    concentration.

    Args:
        c (float): Concentration.

    Returns:
        int: The number of realizations that exist.
    """
    import glob
    fcopy = os.path.join(_copydir,'c{}_*.copy'.format(int(c)))
    ncopy = len(glob.glob(fcopy))
    return ncopy

# ------------------------------------------------------------------------------
# METHODS FOR HANDLING TESTS
# ==========================
def param_test(idstr=None,filename=None,overwrite=False,topdir=_outputdir,
               prefix='',snapfile=None,outputdir=None,version=-1,c=_default_conc,
               squishy=False,squishz=False,decimate=False,
               substr=False,subm=_default_subm,subr=_default_subr,
               subrho=_default_subrho,subc=_default_subc):
    """Creates a parameter file for a test. For an example of a valid parameter 
    file including descriptions of each parameter, please see 'example.param'.

    Args:
        idstr (Optional[str]): String that should be used to identify this test. If not
            provided, one is created by the 'testid' method based on the 
            test parameters.
        filename (Optional[str]): Path to file where parameters should be saved. 
            (default = '[topdir]/[idstr]/[idstr].param')
        overwrite (Optional[bool]): If True, any existing parameter file is 
            overwritten. (default = False)
        topdir (Optional[str]): Path to directory where a separate directory 
            should be created for this run based on idstr. (default = `_outputdir`)
        prefix (Optional[str]): String to start generated idstr with. Only used 
            if idstr is not provided. (default = '')

        snapfile (Optional[str]): Path to file where particle snapshot is for 
            this test. If not provided, the appropriate test snapshot is selected.
        outputdir (Optional[str]): Path to directory where run output should be 
            saved. If not provided, the directory containing the generated 
            parameter file is used.

        version (Optional[int]): test realization. Set to -1 for the fiducial 
            version. (default = -1)
        c (Optional[float]): concentration of run (default = `_default_conc`)
        squishy (Optional[float]): factor that y coords should be squished by 
            (default = False)
        squishz (Optional[float]): factor that z coords should be sauished by 
            (default = False)
        decimate (Optional[int]): factor that particle number should be decimated 
            by (default = False)
        substr (Optional[bool]): If true, the test includes substructure 
            (default = False)

        subm (Optional[float]): mass of subhalo relative to parent's Mvir. 
            Only used if `substr` == True. (default = `_default_subm`)
        subr (Optional[float]): radius of subhalo relative to parent's Rvir. 
            Only used if `substr` == True. (default = `_default_subr`)
        subc (Optional[float]): concentration of substructure. 
            Only used if `substr` == True. (default = `_default_subc`)
        subrho (Optional[float]): central density of subhalo relative to parent's central
            density. Only used if `substr` == True. (default = `_default_subrho`)

    Returns:
        dict: Dictionary containing parameters written to the file.

    """
    # Get base id and file names
    baseid = testid(c=c,version=version)
    if version == -1:
        basesnap = os.path.join(_halodir,baseid+'.snap')
    else:
        basesnap = os.path.join(_copydir,baseid+'.copy')
    # Get string ID for this test 
    if idstr is None:
        idstr = testid(c=c,squishy=squishy,squishz=squishz,decimate=decimate,
                       substr=substr,subm=subm,subr=subr,subc=subc,subrho=subrho,
                       version=version,prefix=prefix)
    if substr:
        idstr_sub = subid(subm=subm,subr=subr,subc=subc,subrho=subrho)
    else:
        idstr_sub = ''
    # Get filenames if not provided
    if topdir is None:
        topdir = _outputdir
    if filename is None:
        filename = os.path.join(topdir,idstr,idstr+'.param')
    if outputdir is None:
        outputdir = os.path.dirname(filename)
    if snapfile is None:
        if substr:
            snapfile = os.path.join(outputdir,idstr+'.snap')
        else:
            snapfile = basesnap
    # Create parameter file if needed
    if not os.path.isfile(filename) or overwrite:
        param = {}
        param['FilePrefix'] = baseid
        param['FileSuffix'] = idstr_sub
        param['NumDivide'] = 1
        param['PeriodicBoundariesOn'] = 0
        param['BoxSize'] = 0.
        param['Border'] = 0.1
        param['PositionFileFormat'] = _halofmt
        param['BgTreebiNskip'] = 0
        param['PositionFile'] = snapfile
        param['OutputDir'] = outputdir
        param['DecimateInputBy'] = decimate if decimate>1 else 1
        param['SquishY'] = squishy if squishy>0 else -1
        param['SquishZ'] = squishz if squishz>0 else -1
        param = voro.make_param(filename,overwrite=overwrite,**param)
    # Otherwise, just read it
    else:
        param = voro.read_param(filename)
        # Fix directory
        if param['OutputDir']!=outputdir:
            print 'Fixing miss-matched dir.'
            param['OutputDir'] = outputdir
            param['PositionFile'] = snapfile
            param = voro.make_param(filename,overwrite=True,**param)
    # Return param
    param['idstr'] = idstr
    param['basesnap'] = basesnap
    param['parfile'] = filename
    return param

def load_test(param=None,Mscl=1.,Rscl=1.,**kwargs):
    """Loads & returns particle information for a test based on parameters. 

    Args:
        param (Optional[dict]): Dictionary of parameters defining a test. 
            If param is not provided, additional keywords are passed to 
            ``param_test`` in order to construct a dictionary of test 
            parameters. (default = None)
        Mscl (Optional[float]): Scale factor that particle masses should be 
            scaled by. (default = 1.)
        Rscl (Optional[float]): Scale factor that particle positions (and 
            volumes) should be scaled by. (default = 1.)
        **kwargs: Parameters that should be used to initialize a test
            parameter file if `param` is not provided.

    Returns:
        dict: With keys:

            * **param** (*dict*): The parameter dictionary for this test.
            * **mass** (*np.ndarray*): (N,) array of particle masses.
            * **pos** (*np.ndarray*): (N,3) array of particle positions.
            * **vol** (*np.ndarray*): (N,) array of particle volumes (only 
                included if the voronoi tesselation output exists).

    Raises:
        ValueError: If the ``param['PositionFile]`` is not a valid path.
    """
    # Get parameters
    if param is None:
        kwargs['overwrite'] = False
        param = param_test(**kwargs)
    out = {'param':param}
    # Get file names
    snapfile = param['PositionFile']
    volfile = voro.namefile('vols',param)
    # Load snapshot
    if os.path.isfile(snapfile):
        out['mass'],out['pos'] = voro.read_snapshot(param)
        out['mass']*=Mscl
        out['pos']*=Rscl
    else:
        raise ValueError('Snapshot does not exist.\n    {}'.format(snapfile))
    # Load volumes
    if os.path.isfile(volfile):
        out['vol'] = voro.read_volume(volfile)
        out['vol']*=(Rscl**3.)
    else:
        print 'Could not load volume file.'
        print '    '+volfile
    # Return
    return out

# ------------------------------------------------------------------------------
# CORE METHODS FOR DOING A SERIES OF TESTS
# ========================================
def series_vallist(series,vlist=None,vlim=None,Nv=10):
    """Creates a list of parameter values for a given series.

    Args:
        series (str): String identifying what parameter should be varied. 
            Currently supported values include:
                * 'conc'       : vary concentration of the halo (5,10,50)
                * 'npart'      : vary number of particles (100 to 1,000,000)
                * 'oblate'     : vary oblateness of the halo (0.3 to 0.9)
                * 'prolate'    : vary prolateness of the halo (0.3 to 0.9)
                * 'triax'      : vary triaxiality of the halo (0.1 to 0.9)
                * 'substr_mass': substructure of varying mass (0.01 to 0.2)
                * 'substr_rsep': substructure at varying radii (0.01 to 0.75)
                * 'substr_conc': substructure of varying concentration (5,10,50)
                * 'substr_rho0': substructure of varying density (0.1 to 1)
        vlist (Optional[list]): list of values to vary parameter specified by 
            series over (if not provided, `vlim` and `Nv` are used to generate 
            this list. Default values for `vlim` are different for each series 
            and are in parenthesis next to the description of each series 
            above).
        vlim (Optional[tuple]): (min,max) values to vary parameter over
        Nv (Optional[int]): Number of parameter variations that the series
            should include. (default = 10)

    Returns:
        list: Values that a parameter for a given series should be varied 
            over.

    Raises:
        ValueError: If `series` is not in `_list_series` or `vlist` is
            provided, but not a list.
    """
    # Check that series is in list of existing values
    if series not in _list_series:
        raise ValueError('Series {} is not in list of '.format(series)+
                         'supported values: {}'.format(_list_series))
    # Create list of values
    if isinstance(vlist,(list,np.ndarray)):
        pass
    elif vlist is None:
        if   series in ['conc','substr_conc']: vlist = _list_conc
        else:
            if vlim is None:
                if series == 'oblate': vlim = (0.3,1.0)
                elif series == 'prolate': vlim = (0.3,1.0)
                elif series == 'triax': vlim = (0.1,1.0)
                elif series == 'npart': vlim = (1,10000)
                elif series == 'substr_mass': vlim = (0.01,0.2)
                elif series == 'substr_rsep': vlim = (0.01,0.75)
                elif series == 'substr_rho0': vlim = (0.1,1.0)
            if series in ['substr_mass','substr_rsep','substr_rho0','npart']:
                vlist = np.logspace(np.log10(vlim[0]),np.log10(vlim[1]),Nv)
            else:
                vlist = np.linspace(vlim[0],vlim[1],Nv)
            if series=='npart':
                vlist = vlist.astype(int)
                vlist = np.array([1,5,10,50,100,500,1000,5000,10000])#,50000])
    else:
        raise ValueError('vlist must be a list or array.')
    # Return
    return list(vlist)

def run_series(series,vlist=None,vlim=None,Nv=10,nfwmeth=_nfwmeth,errors=False,
               ownfw=False,owvoro=False,owparam=False,owsnap=False,
               c=_default_conc,decimate=False,squishy=False,squishz=False,
               ellip=_default_ellip,subm=_default_subm,subr=_default_subr,
               subc=_default_subc,subrho=_default_subrho,topdir=None,
               prefix='',plotflag=True,plotfile=None,plotprof=True,**kwargs):
    """Run a series of tests that vary some halo parameter. 

    Args:
        series (str): Identifies what parameter should be varied. Each
            parameter has additional keywords that may be supplied to
            control the test further. These are decribed below.
            (default = 'conc') Supported values are...
                * 'conc'       : Vary concentration of the halo (5,10,50)
                * 'npart'      : Vary number of particles (100 - 1,000,000)
                * 'oblate'     : Vary oblateness of the halo (0.3 - 0.9)
                * 'prolate'    : Vary prolateness of the halo (0.3 - 0.9)
                * 'triax'      : Vary triaxiality of the halo (0.1 - 0.9)
                * 'substr_mass': Substructure of varying mass (0.01 - 0.2)
                * 'substr_rsep': Substructure at varying radii (0.01 - 0.75)
                * 'substr_conc': Substructure of varying concentration (5,10,50)
                * 'substr_rho0': Substructure of varying central density (0.1 - 1)
        vlist (Optional[list]): List of values to vary parameter specified by 
            series over (if not provided, vlim and Nv are used to generate this 
            list. Default values for vlim are different for each series and are 
            in parenthesis next to the description of each series above).
        vlim (Optional[tuple]): (min,max) values to vary parameter over.
        Nv (Optional[int]): Number of variations that should be tested.
        nfwmeth (Optional[list]): List of methods that should be used to find 
            NFW profiles each test value.
        errors (Optional[int]): Number of versions that should be run in order 
            to get errors. If set to True, all versions that can be found are 
            used. If set to False, no errors are included. (default = False)

        c (Optional[float]): Concentration of parent halo (default = 10)
            Used for all series except 'conc' for obvious reasons.
        decimate (Optional[int]): Factor that particle number should be 
            decimated by. (default = False) Invalid for series 'npart'.
        squishy (Optional[float]): Factor that y coords should be squished by 
            (default = False) Invalid if series in ['oblate','prolate','triax'].
        squishz (Optional[float]): Factor that z coords should be sauished by 
            (default = False) Invalid if series in ['oblate','prolate','triax'].
        ellip (Optional[float]): Ellipticity of halo (default = 0.5)
            Only used for 'triax' series to constrain ratio between the
            largest and smallest halo axes.
        subm (Optional[float]): Mass of subhalo in terms of the parent halo's 
            virial mass (default = 0.1) Only used for 'substr_*' series, but 
            not for 'substr_mass'.
        subr (Optional[float]): Distance of subhalo from the center of the 
            parent in terms of the parent halo's virial radius (default = 0.5)
            Only used for 'substr_*' series, but not for 'substr_rsep'.
        subc (Optional[float]): Concentration of subhalo (default = 50) Only 
            used for 'substr_*' series, but not for 'substr_conc'.
        subrho (Optional[float]): Central density of subhalo relative to 
            parent's central density. Only used for 'substr_*' series, but not 
            for 'substr_rho'. (default = 0.5)

        topdir (Optional[str]): Path to directory where a separate directory 
            should be created for this run based on idstr. 
            (default = `_outputdir`)
        prefix (Optional[str]): String to start generated idstr with. 
            (default = '')
        ownfw (Optional[bool]): If true, NFW parameters for individual runs in 
            this series are overwritten. (default = False)
        owvoro (Optional[bool]): If true, any existing vorovol output is 
            overwritten. (default = False)
        owparam (Optional[bool]): If true, any existing parameter files are 
            overwritten. (default = False)
        owsnap (Optional[bool]): If true, any existing snapshot files are 
            overwritten. (default = False)
        plotflag (Optional[bool]): If true, concentrations for this test series 
            are plotted. (default = True)
        plotfile (Optional[str]): File where profile plot should be saved.
        plotprof (Optional[bool]): If true, profiles for each test are plotted .
            (default = True)

        **kwargs: Additional keywords are passed to `avg_run` for each value in 
            the series.

    Returns:
        dict: NFW parameters for each value. Keys are str(value). The 
            items they refer to are dictionaries returned by voro.get_nfw with 
            input methods from nfwmeth. See `voro.get_nfw` for details on their 
            structure.

    Raises:
        ValueError: If `series` is not in `_list_series`.
        ValueError: If `errors` is not an int or bool.

    """
    import pickle
    # Check that series is in list of existing values
    if series not in _list_series:
        raise ValueError('Series {} is not in list of '.format(series)+
                         'supported values: {}'.format(_list_series))
    # List of parameter values
    vlist = series_vallist(series,vlist=vlist,vlim=vlim,Nv=Nv)
    # Force overwrites in case snapshot or parameter files change
    if owparam or owsnap: owvoro = True
    if owvoro: ownfw = True
    # Get series specific keywords
    kwargs.setdefault('verbose',False)
    kws = dict(c=None,decimate=None,squishy=None,squishz=None,
               subm=None,subr=None,subc=None,subrho=None)
    if series != 'npart': kws['decimate'] = decimate
    if series not in ['triax','prolate','oblate']:
        kws.update(squishy=squishy,squishz=squishz)
    if series != 'conc': kws['c'] = c
    if series == 'triax': kws['squishz'] = ellip
    if series.startswith('substr_'):
        kws['substr'] = True
        if series != 'substr_mass': kws['subm'] = subm
        if series != 'substr_rsep': kws['subr'] = subr
        if series != 'substr_conc': kws['subc'] = subc
        if series != 'substr_rho0': kws['subrho'] = subrho
    # Count errors
    if isinstance(errors,bool):
        if errors:
            if series == 'conc':
                nerror = min([_count_copies(ic) for ic in vlist])
            else:
                nerror = _count_copies(kws['c'])
        else:
            nerror = 0
    elif isinstance(errors,int):
        nerror = errors
    else:
        raise ValueError('errors must be a bool or an int specifying number '+
                         'of iterations that should be used.')
    # Create file names
    if topdir is None: 
        topdir = _outputdir
    series_prefix = prefix+'_'+series
    if nerror>0:
        series_prefix+='{:03d}'.format(nerror)
    seriesid = testid(prefix=series_prefix,**kws)
    if plotfile is None:
        plotfile = os.path.join(topdir,'plot',seriesid+'.png')
    # Add additional keywords
    kws.update(**kwargs)
    # Loop over values
    vlist0 = [] ; data = {}
    for v in vlist:
        # Start
        print 80*'-'
        print 'Running {} = {}...'.format(series,v)
        # Set keywords for varied parameter
        if   series == 'conc': 
            kws['c'] = v
        elif series == 'prolate': 
            kws['squishy'] = v
            kws['squishz'] = v
        elif series == 'oblate':
            kws['squishy'] = False
            kws['squishz'] = v
        elif series == 'triax':
            kws['squishy'] = util.triax2squeeze(v,kws['squishz'])
        elif series == 'npart':
            kws['decimate'] = v
        elif series == 'substr_mass':
            kws['subm'] = v
        elif series == 'substr_rsep':
            kws['subr'] = v
        elif series == 'substr_conc':
            kws['subc'] = v
        elif series == 'substr_rho0':
            kws['subrho'] = v
        # Run test for each error
        code,iout = avg_test(owvoro=owvoro,owparam=owparam,owsnap=owsnap,
                             plotflag=plotprof,topdir=topdir,nerror=nerror,
                             ownfw=ownfw,nfwmeth=nfwmeth,prefix=prefix,**kws)
        print 'result ({}) = {}'.format(v,code)
        # Save
        if code==0:
            vlist0.append(v)
            data[str(v)] = iout
    print 80*'-'
    # Plot
    if plotflag:
        if series=='conc':
            clist = vlist0
        else:
            clist = len(vlist0)*[kws['c']]
        plot_series(series,data,vlist0,clist,plotfile=plotfile,nfwmeth=nfwmeth)
    # Return data
    return data

def plot_series(series,data,vlist=None,clist=None,plotfile=None,nfwmeth=_nfwmeth,
                residuals=True,errorbars=True,errorbars_res=False,legend=False,
                verbose=True):
    """Plot performance of each series as a function of value.

    Args:
        series (str): parameter that this series values
        data (dict): dictionary of data that should be plotted or full path to 
            file containing the data that should be plotted.
        vlist (Optional[list]): list of parameter values that series covers
            (optional if data is a file)
        clist (Optional[list]): list of true concentrations for each value
            (optional if data is a file)
        plotfile (Optional[str]): file where this plot should be saved. If not 
            provided the plot is just displayed and not saved.
        nfwmeth (Optional[list]): list of methods that should be plotted
        residuals (Optional[bool]): if True, residuals are also plotted 
            (default = True)
        errorbars (Optional[bool]): if True, errorbars are included
        legend (Optional[bool]): if True, a legend is put on the top. (default = False)
        verbose (Optional[bool]): if True, information about the series is printed. 
            (default = True)

    Raises:
        ValueError: If there arn't enough colors for all of the values.
        ValueError: If vlist or clist are not provided for the data dictionary.

    """
    import matplotlib.pyplot as plt
    import pickle
    import matplotlib as mpl
    mpl.rcParams['axes.linewidth'] = 1
    yloc_res = plt.MaxNLocator(5)
    ptsize = 1
    labelx = None #-0.1
    titles = {'conc':'Concentration','triax':'Triaxiality',
              'prolate':'Prolate Ellipticity','oblate':'Oblate Ellipticity',
              'substr_conc':'Substructure Concentration',
              'substr_mass':'Substructure Mass',
              'substr_rsep':'Substructure Position',
              'substr_rho0':'Substructure Density',
              'npart':'Number of Particles'}
    # colors = ['#1f78b4', # Color blind safe colors - not very distinct
    #           '#33a02c',
    #           '#a6cee3',
    #           '#b2df8a',
    colors = ['b','r','g','m','c']
    limits = {}
    ticks = {}
    if len(nfwmeth)>len(colors):
        raise ValueError('Only have {} colors '.format(len(colors))+
                         'for {} nfwmeths.'.format(len(nfwmeth)))
    if not errorbars:
        errorbars_res = False
    # Load data if its a file
    if isinstance(data,dict):
        if vlist is None or clist is None:
            raise ValueError('Must provide vlist and clist if data is a dictionary.')
    elif isinstance(data,str):
        fd = open(data,'r')
        data = pickle.load(fd)
        fd.close()
    # Ensure arrays
    vlist = np.array(vlist)
    clist = np.array(clist,dtype=float)
    if series=='npart':
        xlist = np.array([data[str(v)]['voronoi']['N'] for v in vlist])
    else:
        xlist = vlist
    # Check that errorbars is possible
    if errorbars:
        for v in vlist:
            for m in nfwmeth:
                if data[str(v)][m].get('Nstd',1) <= 1:
                    errorbars = False
                    break
    # Set up axes
    plt.clf()
    fig = plt.figure(figsize=(9,4))
    if residuals:
        axs1 = fig.add_axes((.1,.3,.6,.6))
        plt.setp( axs1.get_xticklabels(), visible=False)
        axs2 = fig.add_axes((.1,.1,.6,.2))
        axs2.set_xlabel(titles[series])#series.title())
        axs2.yaxis.set_major_locator(yloc_res)
    else:
        axs1 = plt.subplot(1,1,1)
        axs1.set_xlabel(series.title())
    axs1.set_ylabel('Measured c')
    # Determine scales
    if series in ['npart','substr_mass','substr_rho0']:
        if errorbars:
            plotmeth1 = axs1.errorbar
        else:
            plotmeth1 = axs1.semilogx
        if residuals:
            if errorbars_res:
                plotmeth2 = axs2.errorbar
            else:
                plotmeth2 = axs2.semilogx
    else:
        if errorbars:
            plotmeth1 = axs1.errorbar
        else:
            plotmeth1 = axs1.plot
        if residuals:
            if errorbars_res:
                plotmeth2 = axs2.errorbar
            else:
                plotmeth2 = axs2.plot
    # Plot true values
    plotmeth1(xlist,clist,c='k',ls='--',label='True Value')
    if residuals:
        axs2.axhline(0.,c='k',ls='--',label='True Value')
    xlim = (min(xlist),max(xlist))
    ylim = limits.get(series,(min(clist),max(clist)))
    rmax = 0.0
    # Loop over NFW methods
    for i,k in enumerate(nfwmeth):
        # Set plotting keywords
        ikws = dict(c=colors[i],ls='-',label=k.title(),marker='o',markersize=ptsize,zorder=i-20)
        if k.startswith('voro'): 
            ikws['zorder']=1 # Always plot vorovol on top
        if errorbars:
            ikws['yerr'] = np.array([data[str(v)][k]['c_std'] for v in vlist],float)
            Nstd = np.array([data[str(v)][k]['Nstd'] for v in vlist],int)
        # Plot
        ylist = np.array([data[str(v)][k]['c'] for v in vlist],float)
        Nlist = np.array([data[str(v)][k]['Nstd'] for v in vlist],float)
        ipltout = plotmeth1(xlist,ylist,**ikws)
        ylim = (min(ylim[0],min(ylist)),max(ylim[1],max(ylist)))
        # Fix errorbars (bug in matplotlib)
        if errorbars:
            for icap in ipltout[1]:
                icap.zorder = ipltout[0].zorder
            for ierr in ipltout[2]:
                ierr.zorder = ipltout[0].zorder
        # Plot triaxial base line
        if series == 'triax':
            axs1.axhline(ylist[xlist==1.][0],color=colors[i],ls=':',zorder=i-20+len(nfwmeth))
        elif series == 'npart':
            axs1.axhline(ylist[xlist==1000000][0],color=colors[i],ls=':',zorder=i-20+len(nfwmeth))
        # Do residuals
        if residuals:
            if series == 'triax':
                rlist = (ylist-ylist[xlist==1.][0])/(ylist[xlist==1.][0])
            elif series == 'npart':
                rlist = (ylist-ylist[xlist==1000000][0])/(ylist[xlist==1000000][0])
            else:
                rlist = (ylist-clist)/clist
            rlist*=100. # In percents
            rmax = max(rmax,max(np.abs(rlist)))
            if verbose:
                print 'tests.py @ 625:',k
                print '    {:10s} {:5s} {:5s} {:5s} {:7s} {:5s} {:7s} {:5s}'.format('value','ctrue','cmeas','accu','accuR','prec','precR','inside?')
                for cc,ww,xx,rr,ee in zip(clist,xlist,ylist,rlist,ikws['yerr']):
                    print '    {:10.2f} {:5.2f} {:5.2f} {:5.2f} {:7.4f} {:5.2f} {:7.4f} {}'.format(ww,cc,xx,np.abs(cc-xx),rr,ee,100.*ee/cc,np.abs(cc-xx)<ee)
            if errorbars_res:
                # del ikws['yerr']
                ikws['yerr'] = ikws['yerr']/clist
            else:
                if 'yerr' in ikws: del ikws['yerr']
            iresout = plotmeth2(xlist,rlist,**ikws)
            # Fix residual errorbars (bug in matplotlib)
            if errorbars_res:
                for icap in iresout[1]:
                    icap.zorder = iresout[0].zorder
                for ierr in iresout[2]:
                    ierr.zorder = iresout[0].zorder
    # Limits 
    ypad = 0.1*(ylim[1]-ylim[0])
    ylim = (ylim[0]-ypad,ylim[1]+ypad)
    axs1.set_xlim(xlim)
    axs1.set_ylim(ylim)
    if residuals:
        rpad = 0.2*rmax
        axs2.set_xlim(xlim)
        axs2.set_ylim((-rmax-rpad,rmax+rpad))
        if labelx is not None:
            axs1.yaxis.set_label_coords(labelx,0.5)
            axs2.yaxis.set_label_coords(labelx,0.5)
    # Legend
    if legend:
        leg = axs1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                          ncol=3,mode='expand',borderaxespad=0.,
                          prop={'size':13})
    else:
        leg = None
    # Log for error bars
    if series in ['npart','substr_mass','substr_rho0'] and errorbars:
        axs1.set_xscale('log')
        if residuals and errorbars_res:
            axs2.set_xscale('log')
    # Ticks
    if residuals:
        fig.canvas.draw()
        tics = [item.get_text() for item in axs1.get_yticklabels()]
        if (eval(tics[1])-ylim[0])/(ylim[1]-ylim[0])<0.1:
            tics[1] = ''
            axs1.set_yticklabels(tics)
    # Save or show
    if plotfile is None:
        plt.show()
    else:
        extra_artists = []
        if legend:
            extra_artists.append(leg)
        if not os.path.isdir(os.path.dirname(plotfile)):
            os.mkdir(os.path.dirname(plotfile))
        plt.savefig(plotfile,bbox_inches='tight',
                    bbox_extra_artists=extra_artists)
        print '    '+plotfile
    # Return
    return

# ------------------------------------------------------------------------------
# METHODS FOR DOING SINGLE TESTS
# ==============================
def avg_test(nerror=True,verbose=False,**kwargs):
    """Run multiple realiziations of a test in order to get the mean and 
    standard deviation on the parameters determined by different techinques.

    Args:
        nerror (Optional[int]): Number of test runs that should be averaged 
            over. If True, all existing runs are averaged. (default = True)
        verbose (Optional[bool]): If True, information is printed out for each 
            run. (default = False)
        **kwargs: Additional keywords are passed to run_test.

    Returns:
        code (int): Integer code specifying state of the test. Any value other 
            than 0 indicates that there was an error with all of the runs.
        **out** (*dict*): The same as the output dictionary returned by 
            run_test, with the exception that the values returned are averaged 
            over the different realizations. In addtion, each method dictionary
            will include the standard deviations of these values under the
            keys k+'_std' and 'Nstd', the number of realization used to
            compute the average.

    Raises:
        ValueError: If nerror is not a bool or integer.
        ValueError: If filenames are provided as input. Filenames must be
            generated automatically in order to prevent overwrite.

    """
    # List of errors
    errlist = [-1]
    if isinstance(nerror,bool):
        if nerror:
            nerror = _count_copies(kwargs.get('c',_default_conc))
        else:
            nerror = 0
    elif isinstance(nerror,int):
        pass
    else:
        raise ValueError('nerror must be an int specifying number '+
                         'of iterations that should be used or a bool.')
    # Ensure that files are not overwritten
    if nerror > 0:
        errlist+= range(nerror)
        for x in ['idstr','parfile','snapfile','outfile','plotfile']:
            if kwargs.get(x,None) is not None:
                raise ValueError('Specifying {} for multiple runs '.format(x)+
                                 'will result in overwriting data.')
    # Loop over versions
    codes = []
    results = []
    for version in errlist:
        icode,iresult = run_test(version=version,verbose=verbose,
                                 **copy.deepcopy(kwargs))
        if icode == 0:
            codes.append(icode)
            results.append(iresult)
    # Average
    if len(results)==0:
        code = -999
        out = None
    elif len(results)==1:
        code = codes[0]
        out = results[0]
        for m in out.keys():
            out[m]['Nstd'] = 1
    else:
        code = 0
        out = {}
        for m in results[0].keys():
            out[m] = {'Nstd':len(results)}
            for k in results[0][m].keys():
                vals = np.array([results[i][m][k] for i in range(len(results))])
                out[m][k] = np.mean(vals)
                out[m][k+'_std'] = np.std(vals)
                if k == 'c':
                    print str('{:10s}: '+len(vals)*'{:5.2f} ').format(m,*vals)
    # Return
    return code,out

def run_test(idstr=None,topdir=None,prefix='',nfwmeth=_nfwmeth,verbose=True,
             outputdir=None,
             parfile=None,snapfile=None,exefile=None,outfile=None,
             owvoro=False,owparam=False,owsnap=False,version=-1,
             c=_default_conc,squishy=False,squishz=False,decimate=False,
             substr=False,subm=_default_subm,subr=_default_subr,
             subc=_default_subc,subrho=_default_subrho,**kwargs):
    """Run a test.

    This inclues creating the necessary parameter and snapshot files, running
    the tessellation code, and determining the concentration.

    Args:
        idstr (Optional[str]): String that should be used to identify this 
            test. If not provided, one is created by the 'testid' method based 
            on the test parameters.
        topdir (Optional[str]): Path to directory where a separate directory 
            should be created for this run based on idstr. 
            (default = `_outputdir`)
        prefix (Optional[str]): String to start generated idstr with. Only used 
            if idstr is not provided. (default = '')
        verbose (Optional[bool]): If true, infomation about the run is printed. 
            (default = True)

        outputdir (Optional[str]): Path to directory where run output should be 
            saved. If not provided, the directory containing the generated 
            parameter file is used.
        parfile (Optional[str]): Path to file where parameters should be saved
        snapfile (Optional[str]): Path to file where particle snapshot is for 
            this test. If not provided, the appropriate test snapshot is 
            selected.
        exefile (Optional[str]): Path to vorovol executable
        outfile (Optional[str]): Path to file where vorovol runtime output 
            should be piped
        owvoro (Optional[bool]): If true, any existing vorovol output is 
            overwritten. This sets ownfw to True. (default = False)
        owparam (Optional[bool]): If true, any existing parameter file is 
            overwritten (default = False)
        owsnap (Optional[bool]): If true, any existing snapshot file is 
            overwritten. This also sets owvoro to True. (default = False)

        version (Optional[int]): test realization. Set to -1 for the fiducial 
            version. (default = -1)
        c (Optional[float]): concentration of run (default = `_default_conc`)
        squishy (Optional[float]): factor that y coords should be squished by 
            (default = False)
        squishz (Optional[float]): factor that z coords should be sauished by 
            (default = False)
        decimate (Optional[int]): factor that particle number should be 
            decimated by (default = False)
        substr (Optional[bool]): If true, the test includes substructure 
            (default = False)

        subm (Optional[float]): mass of subhalo relative to parent's Mvir. 
            Only used if `substr` == True. (default = `_default_subm`)
        subr (Optional[float]): radius of subhalo relative to parent's Rvir. 
            Only used if `substr` == True. (default = `_default_subr`)
        subc (Optional[float]): concentration of substructure. 
            Only used if `substr` == True. (default = `_default_subc`)
        subrho (Optional[float]): central density of subhalo relative to 
            parent's central density. Only used if `substr` == True. 
            (default = `_default_subrho`)

        **kwargs: Additional keywords are passed to voro.get_nfw. Some keywords 
            are modified by this method. These include:
        nfwmeth (Optional[list]): List of methods that should be used 
            to find NFW profiles each test value.
        nfwfile (Optional[str]): Path to file where NFW fit should be 
            saved. (default = '[outputdir]/[idstr]_nfw.dat')
        ownfw (Optional[bool]): If True, any existing NFW data file is 
            overwritten. If owvoro is True, this is set to True as well 
            in order to make sure the file reflects the most recent 
            tessellation.
        plotfile (Optional[str]): File where profile plot should be 
            saved. (default = '[outputdir]/[idstr]_profile.png')

    Returns:
        code (int): Code describing the success or failure of the test as 
            returned by `voro.run`.
        **out** (*dict*): Results dictionary returned by `voro.get_nfw`.

    Raises:
        ValueError: If the snapfile provided does not exist and there is not 
            a way to create it (Only substructure snapfiles can be created).

    .. todo:: only accept errors for vornoi in run_test

    """
    # Get parameters
    param = param_test(idstr=idstr,prefix=prefix,filename=parfile,snapfile=snapfile,
                       topdir=topdir,outputdir=outputdir,overwrite=owparam,
                       c=c,squishy=squishy,squishz=squishz,
                       decimate=decimate,substr=substr,subm=subm,subr=subr,
                       subc=subc,subrho=subrho,version=version)
    if owsnap: owvoro = True
    if owvoro: kwargs['ownfw'] = True
    # Update directories and file names based on parameters
    idstr = param['idstr']
    basesnap = param['basesnap']
    outputdir = param['OutputDir']
    snapfile = param['PositionFile']
    if outfile is None:
        outfile = os.path.join(outputdir,idstr+'.out')
    kwargs.setdefault('plotfile',os.path.join(outputdir,idstr+'_profile.png'))
    kwargs.setdefault('nfwfile',os.path.join(outputdir,idstr+'_nfw.dat'))
    # Create snapshot if needed
    if not os.path.isfile(snapfile) or owsnap:
        if substr:
            make_substr(snapfile,basesnap,overwrite=owsnap,
                        subm=subm,subr=subr,subc=subc,subrho=subrho,
                        version=version)
        # If its not substructure, it should already exist
        else:
            raise ValueError('Cannot create snapshot: {}'.format(snapfile))
    # Run tessellation
    code = voro.run(param,exefile=exefile,outfile=outfile,
                    overwrite=owvoro,verbose=verbose)
    # Get and return NFW parameters
    if code == 0:
        try:
            out = voro.get_nfw(param,method=nfwmeth,**kwargs)
            return code,out
        except:
            raise
        #   return -999,None
    else:
        return -999,None
    # Return
    return code,out
        
def make_substr(filename,parfile,overwrite=False,version=-1,trialrun=False,
                subm=_default_subm,subr=_default_subr,subc=_default_subc,
                subrho=_default_subrho):
    """Create snapshot containing substructure.

    Args:
        filename (str): File where snapshot with substructure should be saved.
        parfile (str): File containing parent halo that subhalo should be 
            added to.
        overwrite (Optional[bool]): If true, any existing snapshot is replaced.
            (default = False)
        version (Optional[int]): test realization. Set to -1 for the fiducial 
            version. (default = -1)
        trialrun (Optional[bool]): If true, calculations are done and printed, 
            but a snapshot is not created. (default = False)

        subm (Optional[float]): mass of subhalo relative to parent's Mvir. 
            (default = `_default_subm`)
        subr (Optional[float]): radius of subhalo relative to parent's Rvir. 
            (default = `_default_subr`)
        subc (Optional[float]): concentration of substructure. 
            (default = `_default_subc`)
        subrho (Optional[float]): central density of subhalo relative to 
            parent's central density. (default = `_default_subrho`)

    Raises:
        ValueError: If parfile does not exist.
        ValueError: If a file containing the template for the subhalo does 
            not exist.
        ValueError: The particles in the parent and subhalo have different
            masses.
        ValueError: The generated mass array does not have the correct length.
        ValueError: The generated position array does not have the correct 
            shape.
    """
    # Prevent overwrite
    if os.path.isfile(filename) and not overwrite and not trialrun:
        print 'Specified file already exists and overwrite not set.'
        print '    '+filename
        return
    snapdir = os.path.dirname(filename)
    if not os.path.isdir(snapdir):
        os.mkdir(snapdir)
    # Check for necessary files
    if not os.path.isfile(parfile):
        raise ValueError('Parent halo snapshot does not exist: {}'.format(parfile))
    subid = testid(c=subc,version=version)
    if version==-1:
        subfile = os.path.join(_halodir,subid+'.snap')
    else:
        subfile = os.path.join(_copydir,subid+'.copy')
    if not os.path.isfile(subfile):
        raise ValueError('Sub halo snapshot does not exist: {}'.format(subfile))
    # Load base snapshot
    mass1,pos1 = io.read_snapshot(parfile,format=_halofmt)
    N1 = len(mass1)
    r1 = util.pos2rad(pos1)
    rvir1 = np.max(r1)
    mvir1 = np.sum(mass1)
    rho1 = util.calc_rhoenc(mass1,r1,rvir1/100.)
    # Load substr snapshot
    mass2,pos2 = io.read_snapshot(subfile,format=_halofmt)
    N2 = len(mass2)
    r2 = util.pos2rad(pos2)
    rvir2 = np.max(r2)
    mvir2 = np.sum(mass2)
    rho2 = util.calc_rhoenc(mass2,r2,rvir2/100.)
    # Transform mass by downsampling particles
    if mass1[-1]!=mass2[-1]:
        raise ValueError('Halos have different particle masses. m1={} m2={}'.format(mass1[-1],mass2[-1]))
    msub = subm*mvir1
    mfrac = msub/mvir2
    Nsub = int(mfrac*N2)
    Nfrac = float(Nsub)/float(N2)
    # Shuffle order to ensure halo is evenly sampled
    idx2 = np.arange(0,N2)
    np.random.shuffle(idx2)
    mass_sub = mass2[idx2[:Nsub]]
    # Squeeze halo such that the inner density is set to the correct value
    rho_sub = util.calc_rhoenc(mass_sub,r2[idx2[:Nsub]],rvir2/100.)
    # Squeeze halo such that the mean density is set to the correct value
    # irho = 3*M/(4*pi*(iR**3))
    # frho = 3*M/(4*pi*(fR**3))
    # frho/irho = (iR/fR)**3
    # rfrac = (frho/irho)**(-1./3.)
    # frho = subrho*rho1
    rfrac = ((subrho*rho1)/rho_sub)**(-1./3.)
    pos_sub = rfrac*pos2[idx2[:Nsub],:]
    # Print things
    if trialrun:
        r_sub = util.pos2rad(pos_sub)
        rho_sub0 = util.calc_rhoenc(mass_sub,r_sub,rvir2*rfrac/100.)
        print 'rfrac  = {}'.format(rfrac)
        print 'subrho = {} ({})'.format(rho_sub0/rho1,subrho)
        print 'Primary'
        print '  rvir  = {}'.format(rvir1)
        print '  rho   = {}'.format(rho1)
        print 'Secondary'
        print '  rvir  = {}'.format(rvir2)
        print '  rho   = {}'.format(rho2)
        print 'Subhalo'
        print '  rvir  = {} ({})'.format(r_sub.max(),rvir2*rfrac)
        print '  rho   = {} ({})'.format(util.calc_rhoenc(mass_sub,r_sub,rvir2*rfrac/100.),rho_sub)
    # Shift subhalo to correct position in parent
    pos_sub[:,0]+= subr*rvir1
    # Combine arrays
    N = N1+Nsub
    mass = np.concatenate((mass1,mass_sub))
    pos = np.concatenate((pos1,pos_sub))
    if len(mass)!=N:
        raise ValueError('Mass should have {} elements, but it has {}.'.format(N,len(mass)))
    if pos.shape!=(N,3):
        raise ValueError('Position should have shape ({},3) but it has {}'.format(N,pos.shape))
    # Write to file and return
    if not trialrun:
        io.write_snapshot(filename,mass,pos,overwrite=overwrite,format=_halofmt)
    return
