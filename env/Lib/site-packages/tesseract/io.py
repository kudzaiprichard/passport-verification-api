#!/usr/bin/python
"""
io
==

This module provides the classes and methods necessary for reading and 
writing files.

Attributes:
    _snapshot_formats (dict): Maps format codes to the correct snapshot class.

"""

# Standard packages
import numpy as np
import os

# ------------------------------------------------------------------------------
# LOGGING OF SNAPSHOT FORMATS
_snapshot_formats = {}
def register_snapshot_format(code):
    """Decorator to register snapshot format classes. 

    Args:
        code (int): Code that should be associated with the snapshot format
            class being decorated.

    Raises:
        TypeError: If the class being decorated is not a subclass of the
            `io.Snapshot` class.
        Exception: If the provided code is already associated with an existing 
            snapshot format class.

    """
    def wrapper(f):
        if not issubclass(f,Snapshot):
            raise TypeError('Only Snapshot subclasses can be registered.')
        if code in _snapshot_formats:
            raise Exception('Snapshot format code {} was already '.format(code)+
                            'assigned to {}. '.format(_snapshot_formats[code])+
                            'Cannot also be assigned to {}.'.format(f.__name__))
        f._code = code
        _snapshot_formats[code] = f()
        _snapshot_formats[code]._code = code
        _snapshot_formats[code].__name__ = f.__name__
        return f
    return wrapper
def display_snapshot_formats():
    """Prints information on the registered snapshot formats."""
    print 80*'='
    print '{:4s}  {:20s}  {}'.format('code','name','description')
    print 80*'-'
    for c in sorted(_snapshot_formats.keys()):
        f = _snapshot_formats[c]
        print '{:4d}  {:20s}  {}'.format(c,f.__name__,f.__doc__)
    print 80*'='

# ------------------------------------------------------------------------------
# SUPPORTING IO OBJECTS/METHODS
class cStructDict(object):
    """Provides easy read/write of structured data to/from a dictionary."""
    def __init__(self,fields):
        """Initialize object with specified fields.

        Args: 
            fields (List[tuple]): Each tuple contains the name and format of
                field. The order of the tuples specifies the order that the 
                fields should be read from or written to a file.

        """
        self._fields = fields
        self._keys = []
        self._format_list = []
        self._format = ''
        for f in fields:
            self._keys.append(f[0])
            self._format_list.append(f[1])
            self._format+=f[1]
    @property
    def fields(self): 
        """List of field (name,format) tuples."""
        return self._fields
    @property
    def keys(self): 
        """List of field names."""
        return self._keys
    @property
    def format_list(self): 
        """List of field format codes."""
        return self._format_list
    @property
    def format(self): 
        """Format string for the entire data structure."""
        return self._format
    @property
    def size(self):
        """Size of the entire data structure in bytes."""
        if not hasattr(self,'_size'):
            import struct
            self._size = struct.calcsize(self.format)
        return self._size
    def read(self,fd):
        """Reads one instance of data structure from a file.

        Args:
            fd (file): File that data structure should be read from.

        Returns:
            dict: Data structure with fields as key,value pairs.

        """
        return self.unpack(fd.read(self.size))
    def unpack(self,string):
        """Unpacks one instance of data structure from a string.

        Args:
            string (str): String that data structure should be unpacked from.

        Returns:
            dict: Data structure with fields as key,value pairs.

        """
        import struct
        # Get list of values
        values = struct.unpack(self.format,string)
        # Create dictionary
        out = {} ; i = 0
        for f in self.fields:
            f_key = f[0]
            f_fmt = f[1]
            f_len = len(f_fmt)
            # Single value
            if f_len == 1:
                out[f_key] = values[i]
                i+=1
            # Multiple values
            else:
                out[f_key] = values[i:(i+f_len)]
                if len(set(f_fmt))==1:
                    out[f_key] = np.array(out[f_key])
                i+= f_len
        # Return dictionary
        return out

# ------------------------------------------------------------------------------
# GENERAL SNAPSHOT FILE HANDLING
class Snapshot(object):
    """Base class for snapshot formats."""
    @property
    def code(self):
        """Integer code used to reference the snapshot format."""
        if not hasattr(self,'_code'):
            raise AttributeError('This snapshot type does not have a code.')
        return self._code
    def read(self,*args,**kwargs):
        """Dummy read method to return error."""
        raise Exception('This snapshot type does not have a read method.')
    def write(self,*args,**kwargs):
        """Dummy write method to return error."""
        raise Exception('This snapshot type does not have a write method.')
    def load(self,*args,**kwargs):
        """Alias for read."""
        return self.read(*args,**kwargs)
    def parse_voropar(self,param):
        """Maps Vorovol parameters onto input keywords for this class's 
        read/write methods.

        Args:
            param (dict): Vorovol parameters. This snapshot class does not use 
                any of the Vorovol parameters.

        Returns:
            dict: Keyword arguments for this class's read/write methods.

        """
        return {}

def read_snapshot(filename,format=0,**kwargs):
    """Read masses and positions from a snapshot.

    Args:
        filename (str): Full path to the file that should be read.
        format (Optional[int]): Code for the type of snapshot. See output from
            `display_snapshot_formats` for a description of the available 
            codes. (default = 0)
        **kwargs: Additional keywords are passed to the appropriate method for 
            reading.

    Returns:
        mass (np.ndarray): (N,) Particle masses.
        pos (np.ndarray): (N,3) Particle positions.


    Raises:
        ValueError: If the format is not supported.

    """
    if format in _snapshot_formats:
        out = _snapshot_formats[format].read(filename,**kwargs)
    else:
        raise ValueError('No snapshot format registered with code '+
                         '{}.'.format(format))
    return out

def write_snapshot(filename,mass,pos,format=0,**kwargs):
    """Write particle masses and positions to a snapshot.

    Args:
        filename (str): Full path to the file that should be written.
        mass (np.ndarray): (N,) Particle masses.
        pos (np.ndarray): (N,3) Particle positions.
        format (Optional[int]): Integer specifying the type of snapshot. See 
            output from `display_snapshot_formats` for a description of the 
            available codes. (default = 0)
        **kwargs: Additional keywords are passed to the appropriate method for 
            writing.

    Raises:
        ValueError: If the format is not supported.

    """
    if format in _snapshot_formats:
        out = _snapshot_formats[format].write(filename,mass,pos,**kwargs)
    else:
        raise ValueError('No snapshot format registered with code '+
                         '{}.'.format(format))
    return out

def convert_snapshot(filename1,format1,filename2,format2,
                     overwrite=False,**kwargs):
    """Convert one snapshot type into another. The old snapshot is not removed. 

    Args:
        filename1 (str): Full path to the source snapshot that should be read.
        format1 (int): Integer specifying snapshot type of filename1. See
            output from `display_snapshot_formats` for a description of the 
            available codes. 
        filename2 (str): Full path to the destination snapshot that should be 
            created.
        format2 (int): Integer specifying snapshot type of filename2. See
            output from `display_snapshot_formats` for a description of the 
            available codes. 
        overwrite (Optional[bool]): Set to True if existing filename2 should be 
            overwritten. (default = False)
        **kwargs: Additional keywords are passed to the appropriate method for 
            reading.

    Raises:
        Exception: If `filename2` exists and `overwrite` is False.

    """
    # Prevent overwrite
    if os.path.isfile(filename2) and not overwrite:
        raise Exception('Destination snapshot exists and overwrite is not set.')
    # Read in masses and positions from source file
    mass,pos = read_snapshot(filename1,format=format1,**kwargs)
    # Write masses and positions to the destination file
    write_snapshot(filename2,format=format2,overwrite=overwrite)
    # Return
    return

# ------------------------------------------------------------------------------
# F77 UNFORMATED BINARY SNAPSHOT FORMAT
@register_snapshot_format(0)
class UnformattedBinary(Snapshot):
    """Unformatted f77 binary snapshot class."""
    def read(self,*args,**kwargs):
        """See `io.read_unfbi77`."""
        return read_unfbi77(*args,**kwargs)
    def write(self,*args,**kwargs):
        """See `write_unfbi77`."""
        return write_unfbi77(*args,**kwargs)
    def parse_voropar(self,param):
        """Maps Vorovol parameters onto input keywords for this class's 
        read/write methods.

        Args:
            param (dict): Vorovol parameters. This snapshot class performs the
                following parameter to keyword mappings:
                    Unfbi77ArrayType -> dtype (str): Array data type.

        Returns:
            dict: Keyword arguments for this class's read/write methods.

        """
        kwargs = {}
        if 'Unfbi77ArrayType' in param:
            kwargs['dtype'] = param['Unfbi77ArrayType']
        return kwargs

def read_unfbi77(filename,return_npart=False,dtype='float32'):
    """Read unformatted f77 binary snapshot.

    Args:
        filename (str): Full path to file that should be read.
        return_npart (Optional[bool]): If True, only the number of particles 
            is read from the file (default = False).
        dtype (str): Data type of the mass and postion arrays in the snapshot. 
            vorovol currently only reads in float32 arrays. 
            (default = 'float32')

    Returns:
        mass (np.ndarray): (N,) Particle masses.
        pos (np.ndarray): (N,3) Particle positions.

    Raises:
        TypeError: If a type other than `np.float32` is provided.
        IOError: If the 4 bytes before any block does not match the size of the
            block that is read in.

    """
    import struct
    # Check data type
    dtype = np.dtype(dtype)
    if dtype!=np.dtype('float32'):
        raise TypeError("vorvol assumes arrays are 32-bit floats.")
    # Open file
    fd = open(filename,'rb')
    # Read in number of particles
    recl = struct.unpack('i',fd.read(4))[0]
    if recl != np.dtype('int32').itemsize:
        raise IOError('Error reading number of particles from file.')
    nout = struct.unpack('i',fd.read(recl))[0]
    recl = struct.unpack('i',fd.read(4))[0]
    if return_npart: return nout
    # Read in masses
    recl = struct.unpack('i',fd.read(4))[0]
    if recl != (nout*dtype.itemsize):
        raise IOError('Error reading x positions from file.')
    mass = np.fromfile(fd,dtype=dtype,count=nout)
    recl = struct.unpack('i',fd.read(4))[0]
    # Read in x,y,z positions
    pos = np.zeros((nout,3),dtype=np.float32)
    # X
    recl = struct.unpack('i',fd.read(4))[0]
    if recl != (nout*dtype.itemsize):
        raise IOError('Error reading x positions from file.')
    pos[:,0] = np.fromfile(fd,dtype=dtype,count=nout)
    recl = struct.unpack('i',fd.read(4))[0]
    # Y
    recl = struct.unpack('i',fd.read(4))[0]
    if recl != (nout*dtype.itemsize):
        raise IOError('Error reading x positions from file.')
    pos[:,1] = np.fromfile(fd,dtype=dtype,count=nout)
    recl = struct.unpack('i',fd.read(4))[0]
    # Z
    recl = struct.unpack('i',fd.read(4))[0]
    if recl != (nout*dtype.itemsize):
        raise IOError('Error reading x positions from file.')
    pos[:,2] = np.fromfile(fd,dtype=dtype,count=nout)
    recl = struct.unpack('i',fd.read(4))[0]
    # Close and return
    fd.close()
    return mass,pos

def write_unfbi77(filename,mass,pos,overwrite=False):
    """Write unformated f77 binary snapshot.

    Args:
        filename (str): Full path to file that data should be written to.
        mass (np.ndarray): (N,) Particle masses.
        pos (np.ndarray): (N,3) Particle positions.
        overwrite (Optional[bool]): If True and `filename` already exists, it 
            is overwritten. (default = False)

    Raises:
        TypeError: If mass or pos does not have a data type of float32.

    """
    import struct
    # Prevent overwrite
    if os.path.isfile(filename) and not overwrite:
        print 'Specified file already exists and overwrite not set.'
        print '    '+filename
        return
    # Check data type
    dtype = np.dtype('float32')
    if mass.dtype!=dtype or pos.dtype!=dtype:
        raise TypeError("vorvol assumes arrays are 32-bit floats.")
    # Open file
    fd = open(filename,'wb')
    # Write number of particles
    nout = len(mass)
    fd.write(struct.pack('i',np.dtype('int32').itemsize))
    fd.write(struct.pack('i',nout))
    fd.write(struct.pack('i',np.dtype('int32').itemsize))
    # Write masses
    fd.write(struct.pack('i',nout*dtype.itemsize))
    mass.tofile(fd)
    fd.write(struct.pack('i',nout*dtype.itemsize))
    # X Positions
    fd.write(struct.pack('i',nout*dtype.itemsize))
    pos[:,0].tofile(fd)
    fd.write(struct.pack('i',nout*dtype.itemsize))
    # Y Positions
    fd.write(struct.pack('i',nout*dtype.itemsize))
    pos[:,1].tofile(fd)
    fd.write(struct.pack('i',nout*dtype.itemsize))
    # Z Positions
    fd.write(struct.pack('i',nout*dtype.itemsize))
    pos[:,2].tofile(fd)
    fd.write(struct.pack('i',nout*dtype.itemsize))
    # Close and return
    fd.close()
    return

# ------------------------------------------------------------------------------
# GADGET SNAPSHOT FORMAT
@register_snapshot_format(1)
class Gadget2(Snapshot):
    """Gadget2 Type 1 Snapshot"""
    def read(self,*args,**kwargs):
        """Reads data from a Gadget snapshot.
        
        Args:
            format (Optional[int]): Format code specifying what format the 
                Gadget snapshot is in. This code determines what function any 
                additional keywords are passed to. (default = 1) Currently 
                supported values include:
                    1: See `io.read_gadget2_binary1` for details.
            *args: Additional arguments are passed to the appropriate method.
            *kwargs: Additional keyword arguments are passed to the appropriate 
                method.

        Returns:
            mass (np.ndarray): (N,) Particle masses.
            pos (np.ndarray): (N,3) Particle positions.

        Raises:
            ValueError: If specified format is not supported.

        """
        format = kwargs.pop('format',1)
        # Simple binary format
        if format == 1:
            return read_gadget2_binary1(*args,**kwargs)
        # 'Convenient' binary format
        elif format == 2:
            raise ValueError('Gadget snapshot format {} '.format(format)+
                             '(binary2) is not currently supported.')
        # HDF5 format
        elif format == 3:
            raise ValueError('Gadget snapshot format {} '.format(format)+
                             '(hdf5) is not currently supported.')
        # Other
        else:
            raise ValueError('{} is not a valid Gadget-2 '.format(format)+
                             'snapshot format.')
    def parse_voropar(self,param):
        """Maps Vorovol parameters onto input keywords for this class's 
        read/write methods.

        Args:
            param (dict): Vorovol parameters. This snapshot class performs the
                following parameter to keyword mappings:
                    ParticleType -> ptype (int): Particle type.

        Returns:
            dict: Keyword arguments for this class's read/write methods.

        """
        kwargs = {}
        if 'ParticleType' in param:
            kwargs['ptype'] = param['ParticleType']
        return kwargs

def read_gadget2_binary1(filename,ptype=-1,return_npart=False,
                         return_header=False):
    """Read Gadget binary files.

    Args:
        filename (str): Full path to file that should be read.
        ptype (Optional[int,list,tuple,np.ndarray]): Code, or series of codes
            specifying what particle type should be loaded. Supported values 
            include:
               -1: All particles
                0: Gas particles
                1: Dark matter particles
                2: Disk particles
                3: Bulge particles
                4: Star particles
                5: Boundary particles
        return_npart (Optional[bool]): If True, only the number of particles 
            is read from the file (default = False).
        return_header (Optional[bool]): If True, only the file header is read
            from the file (default = False).

    Returns:
        mass (np.ndarray): (N,) Particle masses.
        pos (np.ndarray): (N,3) Particle positions.

    Raises:
        TypeError: If `ptype` is not an integer or series of integers.
        IOError: If file is not read propertly.

    """
    # Set list of particle types
    if isinstance(ptype,int):
        if ptype in range(6):
            typelist = np.array([ptype])
        else:
            typelist = np.arange(6)
    elif isinstance(ptype,(list,tuple)):
        typelist = np.array(ptype)
    elif isinstance(ptype,np.ndarray):
        pass
    else:
        raise TypeError('Unrecognized value for ptype: {}'.format(ptype))
    # Check for multiple files
    if os.path.isfile(filename):
        try:
            idxext = filename0.rindex('.')+1
            ext0 = float(filename[idxext:])
            filebase = filename[:idxext]+'{}'
            filename0 = filebase.format(0)
        except:
            filename0 = filename
    else:
        if os.path.isfile(filename+'.0'):
            filebase = filename+'.{}'
            filename0 = filebase.format(0)
        else:
            filename0 = filename
    # Open first file
    try:
        fd = open(filename0,'rb')
    except IOError:
        raise IOError('Could not open file: {}'.format(filename0))
    # Read in header from first file
    header_struct = GadgetHeaderStruct()
    hint = struct.unpack('i',fd.read(4))[0]
    if hint != header_struct.size:
        raise IOError('Expected block of size {}, '.format(header_struct.size)+
                      'but block header indicates {}.'.format(hint))
    header0 = header_struct.read(fd)
    fint = struct.unpack('i',fd.read(4))[0]
    if fint != header_struct.size:
        raise IOError('Expected block of size {}, '.format(header_struct.size)+
                      'but block footer indicates {}.'.format(fint))
    fd.close()
    if return_header:
        return header0
    # Count particles & return if specified
    nout = 0 ; pc = np.zeros(6) ; mc = np.zeros(6)
    for t in typelist:
        pc[t] = nout
        mc[t] = nout
        nout += header0['npart_total'][t]
    if return_npart:
        return nout
    # Preallocate
    pos = np.zeros((nout,3),dtype=float)
    mass = np.zeros((nout,),dtype=float)
    # Loop over files
    for i in range(header0['num_files']):
        # Open file 
        fd = open(filebase.format(i),'rb')
        # Read header
        if i==0:
            header_i = header0
            fd.seek(4+256+4,0)
        else:
            head_size = header_struct.size
            hint = struct.unpack('i',fd.read(4))[0]
            if hint != head_size:
                raise IOError('Expected block of size {}, '.format(head_size)+
                              'but block header indicates {}.'.format(hint))
            header_i = header_struct.read(fd)
            fint = struct.unpack('i',fd.read(4))[0]
            if fint != head_size:
                raise IOError('Expected block of size {}, '.format(head_size)+
                              'but block footer indicates {}.'.format(fint))
        # Count particles in arrays
        npos_i = 0 ; nmass_i = 0
        for t in typelist:
            npos_i+= header_i['npart'][t]
            if header_i['massarr'][t]==0:
                nmass_i+= header_i['npart'][t]
        # Positions
        pos_size = npos_i*4 # Array of floats
        hint = struct.unpack('i',fd.read(4))[0]
        if hint != pos_size:
            raise IOError('Expected block of size {}, '.format(pos_size)+
                          'but block header indicates {}.'.format(hint))
        for t in range(6):
            ntyp = header_i['npart'][t]
            if ntyp == 0: continue
            if t in typelist:
                pos[pc[t]:(pc[t]+ntyp),:] = \
                    np.fromfile(fd,dtype=np.float32,count=3*ntyp).reshape( \
                    (ntyp,3),order='C')
                pc[t]+= ntyp
            else:
                fd.seek(ntyp*4*3,1)
        fint = struct.unpack('i',fd.read(4))[0]
        if fint != pos_size:
            raise IOError('Expected block of size {}, '.format(pos_size)+
                          'but block footer indicates {}.'.format(fint))
        # Masses
        fd.seek(4 + 3*4*npos_i + 4 +   # velocity block
                4 +   4*npos_i + 4 ,1) # particle ID block
        mass_size = nmass_i*4 # Array of floats
        hint = struct.unpack('i',fd.read(4))[0]
        if hint != mass_size:
            raise IOError('Expected block of size {}, '.format(mass_size)+
                          'but block header indicates {}.'.format(hint))
        for t in range(6):
            ntyp = header_i['npart'][t]
            if ntyp == 0: continue
            if t in typelist:
                if header_i['massarr'][t]!=0:
                    mass[mc[t]:(mc[t]+ntyp)] = header_i['massarr'][t]
                else:
                    mass[mc[t]:(mc[t]+ntyp)] = np.fromfile(fd,dtype=np.float32,
                                                     count=ntyp)
                mc[t]+= ntyp
            else:
                fd.seek(ntyp*4,1)
        fint = struct.unpack('i',fd.read(4))[0]
        if fint != mass_size:
            raise IOError('Expected block of size {}, '.format(mass_size)+
                          'but block footer indicates {}.'.format(fint))
        # Close file
        fd.close()
    # Return
    return mass,pos

class GadgetHeaderStruct(cStructDict):
    """Gadget2 header structure."""
    # Sice should be 256
    def __init__(self):
        fields = [('npart'            ,'IIIIII'),
                  ('massarr'          ,'dddddd'),
                  ('time'             ,'d'),
                  ('redshift'         ,'d'),
                  ('flag_sfr'         ,'i'),
                  ('flag_feedback'    ,'i'),
                  ('npart_total'      ,'IIIIII'),
                  ('flag_cooling'     ,'i'),
                  ('num_files'        ,'i'),
                  ('box_size'         ,'d'),
                  ('Omega0'           ,'d'),
                  ('OmegaLambda'      ,'d'),
                  ('Hubble0'          ,'d'),
                  ('flag_stellarage'  ,'i'),
                  ('flag_metals'      ,'i'),
                  ('NallHW'           ,'IIIIII'),
                  ('flag_entropy_instead_u','i'),
                  ('flag_doubleprecision'  ,'i'),
                  ('flag_ic_info'     ,'i'),
                  ('lpt_scalingfactor','f'),
                  ('padding'          ,'x'*(256-((4*28)+(8*12))))]
        super(GadgetHeaderStruct,self).__init__(fields)

# ------------------------------------------------------------------------------
# BUILDGAL TREEBI FILES
@register_snapshot_format(2)
class BuildgalTreebi(Snapshot):
    """Buildgal TREEBI snapshot."""
    def read(self,*args,**kwargs):
        """See `io.read_bgtreebi`."""
        return read_bgtreebi(*args,**kwargs)
    def write(self,*args,**kwargs):
        """See `io.write_bgtreebi`."""
        return write_bgtreebi(*args,**kwargs)
    def parse_voropar(self,param):
        """Maps Vorovol parameters onto input keywords for this class's 
        read/write methods.

        Args:
            param (dict): Vorovol parameters. This snapshot class performs the
                following parameter to keyword mappings:
                    BgTreebiNskip -> nskip (int): Number of particles to skip.

        Returns:
            dict: Keyword arguments for this class's read/write methods.

        """
        kwargs = {}
        if 'BgTreebiNskip' in param:
            kwargs['nskip'] = param['BgTreebiNskip']
        return kwargs

def read_bgtreebi(filename,nskip=0,return_npart=False):
    """Read Buildgal TREEBI snapshots.

    Args:
        filename (str): Full path to file that should be read.
        nskip (Optional[int]): Number of particles that should be skipped over 
            at the beginning of the file. (default = 0)
        return_npart (Optional[bool]): If True, only the number of particles 
            is read from the file. (default = False)

    Returns:
        mass (np.ndarray): (N,) Particle masses.
        pos (np.ndarray): (N,3) Particle positions.

    """
    fd = open(filename,'r')
    # Read in header
    headline = fd.readline().strip()
    headout = headline.split()
    if headout[0] == '******':
        nout = 1000001
    else:
        nout = int(float(headout[0]))
    dim = int(float(fd.readline().strip()))
    time = float(fd.readline().strip())
    if return_npart: return nout-nskip
    # Allocate
    mass = np.zeros(nout,dtype=np.float32)
    pos = np.zeros((nout,dim),dtype=np.float32)
    # Read masses
    for i in range(nout):
        mass[i] = np.float32(fd.readline().strip())
    # Read positions
    for i in range(nout):
        posstr = fd.readline().strip().split()
        pos[i,0] = np.float32(posstr[0])
        pos[i,1] = np.float32(posstr[1])
        pos[i,2] = np.float32(posstr[2])
    # Close and return
    fd.close()
    return mass[nskip:],pos[nskip:,:]

def write_bgtreebi(filename,mass,pos,overwrite=False):
    """Write Buildgal TREEBI snapshot.

    Args:
        filename (str): Full path to file that data should be written to.
        mass (np.ndarray): (N,) Particle masses.
        pos (np.ndarray): (N,3) Particle positions.
        overwrite (Optional[bool]): If True and `filename` already exists, it 
            is overwritten. (default = False)

    Raises:
        Exception: If the sizes of `mass` and `pos` do not agree.

    """
    # Prevent overwrite
    if os.path.isfile(filename) and not overwrite:
        print 'Specified file already exists and overwrite not set.'
        print '    '+filename
        return
    # Check sizes
    if len(mass)!=pos.shape[0]:
        raise Exception('mass has {} elements, pos has shape {}'.format(len(mass),pos.shape))
    N,dim = pos.shape
    # Open file
    fd = open(filename,'w')
    # Write header
    fd.write(' {:d}     {:d}     {:d}\n'.format(N,0,0))
    fd.write('      {:d}\n'.format(dim))
    fd.write('  {:11.6E}\n'.format(0.))
    # Write masses
    for i in range(N):
        fd.write('  {:11.6E}\n'.format(float(mass[i])))
    # Write positions
    for i in range(N):
        fd.write('  {:11.6E} {:11.6E} {:11.6E}\n'.format(float(pos[i,0]),
                                                         float(pos[i,1]),
                                                         float(pos[i,2])))
    # Close file ane return
    fd.close()
    return

# ------------------------------------------------------------------------------
# BGC2 HALO CATALOGUE
@register_snapshot_format(3)
class Bgc2HaloCatalogue(Snapshot):
    """BGC2 Halo Catalogue"""
    def read(self,*args,**kwargs):
        """See `io.read_bgc2halo`."""
        return read_bgc2halo(*args,**kwargs)
    def parse_voropar(self,param):
        """Maps Vorovol parameters onto input keywords for this class's 
        read/write methods.

        Args:
            param (dict): Vorovol parameters. This snapshot class performs the
                following parameter to keyword mappings:
                    Bgc2HaloId -> haloid (int): Halo ID.

        Returns:
            dict: Keyword arguments for this class's read/write methods.

        """
        kwargs = {}
        if 'Bgc2HaloId' in param:
            kwargs['haloid'] = param['Bgc2HaloId']
        return kwargs


def read_bgc2halo(filename0,haloid=-1,return_npart=False,return_header=False):
    """Read BGC2 halo catalogue.

    Args:
        filename (str): Full path to file that should be read.
        haloid (Optional[int]): ID of halo that should be loaded. If -1, all 
            particles are loaded. (default = -1)
        return_npart (Optional[bool]): If True, only the number of particles 
            is read from the file (default = False).
        return_header (Optional[bool]): If True, only the file header is read
            from the file (default = False).

    Returns:
        mass (np.ndarray): (N,) Particle masses.
        pos (np.ndarray): (N,3) Particle positions.

    Raises:
        IOError: If file cannot be read correctly.

    .. todo:: preallocate groups

    """
    import re,struct
    # Get base file name
    idxext = filename0.rindex('.')-5
    if re.search('\.[0-9][0-9][0-9][0-9]\.bgc2',filename0[idxext:]) is None:
        filename = filename0
    else:
        filename = filename0[:idxext]
    # Open first file
    try:
        fd = open(filename,'rb')
    except IOError:
        try:
            fd = open(filename+'.{:04d}.bgc2'.format(0),'rb')
        except IOError:
            raise IOError('Could not open file: {}'.format(filename0))
    # Read in header from first file
    header_struct = Bgc2HeaderStruct()
    hint = struct.unpack('i',fd.read(4))[0]
    if hint != header_struct.size:
        raise IOError('Expected block of size {}, '.format(header_struct.size)+
                      'but block header indicates {}.'.format(hint))
    header = header_struct.read(fd)
    fint = struct.unpack('i',fd.read(4))[0]
    if fint != header_struct.size:
        raise IOError('Expected block of size {}, '.format(header_struct.size)+
                      'but block footer indicates {}.'.format(fint))
    fd.close()
    if return_header:
        return header
    # Select correct formats
    group_struct = Bgc2GroupStruct(header['format_group_data'])
    part_struct = Bgc2PartStruct(header['format_part_data'])
    # Loop over files, reading in group info
    nfiles = header['num_files']
    nout = 0 ; ng = 0 ; gtot = 0
    groups = [] 
    #    groups = header['ngroups_total']*[0]
    for i in range(nfiles):
        # Open the right file
        if nfiles > 1:
            fname = '{}.{:04d}.bgc2'.format(filename,i)
        else:
            fname = filename
        fd = open(fname,'rb')
        # Read header
        hint = struct.unpack('i',fd.read(4))[0]
        if hint != header_struct.size:
            raise IOError('Expected block of size '+
                          '{}, '.format(header_struct.size)+
                          'but block header indicates {}.'.format(hint))
        iheader = header_struct.read(fd)
        fint = struct.unpack('i',fd.read(4))[0]
        if fint != header_struct.size:
            raise IOError('Expected block of size '+
                          '{}, '.format(header_struct.size)+
                          'but block footer indicates {}.'.format(fint))
        # Read groups
        hint = struct.unpack('i',fd.read(4))[0]
        if hint != iheader['ngroups']*group_struct.size:
            raise IOError('Expected block of size '+
                          '{}, '.format(iheader['ngroups']*group_struct.size)+
                          'but block header indicates {}.'.format(hint))
        for g in range(iheader['ngroups']):
            groups.append(group_struct.read(fd))
            #groups[gtot] = group_struct.read(fd)
            if haloid < 0: haloid = groups[gtot]['id']
            if groups[gtot]['id'] == haloid:
                nout+= groups[gtot]['npart']
                ng+=1
            gtot+=1
        fint = struct.unpack('i',fd.read(4))[0]
        if fint != iheader['ngroups']*group_struct.size:
            raise IOError('Expected block of size '+
                          '{}, '.format(iheader['ngroups']*group_struct.size)+
                          'but block footer indicates {}.'.format(fint))
        # Close the file
        fd.close()
    # Return number of particles
    if return_npart:
        return nout
    # Allocate
    mass = np.ones(nout,dtype=np.float32)*header['part_mass']
    pos = np.zeros((nout,3),dtype=np.float32)
    # Loop over files, collecting particle info
    gtot = 0 ; ntot = 0
    for i in range(nfiles):
        # Open the right file
        if nfiles > 1:
            fname = '{}.{:04d}.bgc2'.format(filename,i)
        else:
            fname = filename
        fd = open(fname,'rb')
        # Read header
        hint = struct.unpack('i',fd.read(4))[0]
        if hint != header_struct.size:
            raise IOError('Expected block of size '+
                          '{}, '.format(header_struct.size)+
                          'but block header indicates {}.'.format(hint))
        iheader = header_struct.read(fd)
        fint = struct.unpack('i',fd.read(4))[0]
        if fint != header_struct.size:
            raise IOError('Expected block of size '+
                          '{}, '.format(header_struct.size)+
                          'but block footer indicates {}.'.format(fint))
        # Skip group data
        fd.seek(4+iheader['ngroups']*group_struct.size+4,1)
        # Loop over groups
        for g in range(iheader['ngroups']):
            # Group header
            hint = struct.unpack('i',fd.read(4))[0]
            if hint != groups[gtot]['npart']*part_struct.size:
                raise IOError('Expected block of size '+
                              '{}, '.format(groups[gtot]['npart']*part_struct.size)+
                              'but block header indicates {}.'.format(hint))
            # Read only particles in this halo
            if groups[gtot]['id'] == haloid:
                for n in range(groups[gtot]['npart']):
                    ipart = part_struct.read(fd)
                    pos[ntot,:] = ipart['pos']
                    ntot+=1
            else:
                fd.seek(groups[gtot]['npart']*part_struct.size,1)
            # Group footer
            fint = struct.unpack('i',fd.read(4))[0]
            if fint != groups[gtot]['npart']*part_struct.size:
                raise IOError('Expected block of size '+
                              '{}, '.format(groups[gtot]['npart']*part_struct.size)+
                              'but block footer indicates {}.'.format(fint))
            # Count groups
            gtot+=1
        # Close this file
        fd.close()
    return mass,pos

class Bgc2HeaderStruct(cStructDict):
    """BGC2 header structure"""
    # Size should be 1024
    def __init__(self):
        fields = [('magic'            ,'Q'), # A magic number to identify this as a BGC file.
                  ('version'          ,'q'), # File version number.
                  ('num_files'        ,'q'), # number of files output is distributed into
                  ('file_id'          ,'q'), # this file's ID (number) if multiple files are output
                  ('snapshot'         ,'q'), # Snapshot ID
                  ('format_group_data','q'), # output group data format identifier
                  ('format_part_data' ,'q'), # output particle data format identifier
                  ('group_type'       ,'q'), # FOF, SO, etc
                  ('ngroups'          ,'q'), # number of groups stored LOCALLY in this file
                  ('ngroups_total'    ,'q'), # number of groups stored GLOBALLY over all output BGC files
                  ('npart'            ,'q'), # number of particles bound in groups LOCALLY in this file
                  ('npart_total'      ,'q'), # number of particles bound in groups GLOBALLY over all output BGC files
                  ('npart_orig'       ,'q'), # number of particles from original simulation input
                  ('max_npart'        ,'q'), # maximum number of particles in one group LOCALLY in this file
                  ('max_npart_total'  ,'q'), # maximum number of particles in one group GLOBALLY over all output BGC files
                  ('min_group_part'   ,'q'), # minimum number of particles in a group
                  ('valid_part_ids'   ,'q'), # valid particle IDs mean they match input snapshot
                  ('linkinglength'    ,'d'), # for FOF halos, what linking length is used
                  ('overdensity'      ,'d'), # mostly SO: overdensity with respect to mean
                  ('time'             ,'d'), # time of the input snapshot 
                  ('redshift'         ,'d'), # redshift of the input snapshot
                  ('box_size'         ,'d'), # input BoxSize
                  ('box_min'          ,'ddd'), # alternative to center, Gadget assumes (0,0,0)
                  ('bounds'           ,'dddddd'), # Spatial bounds of the halo centers contained in this file
                  ('part_mass'        ,'d'), # mass of particles if only one
                  ('Omega0'           ,'d'), # Matter density at z=0
                  ('OmegaLambda'      ,'d'), # Dark energy density at z=0
                  ('Hubble0'          ,'d'), # NOT ALWAYS SET
                  ('GravConst'        ,'d'), # NOT ALWAYS SET
                  ('padding'          ,'x'*(1024-(36*8)))]
        super(Bgc2HeaderStruct,self).__init__(fields)

class Bgc2GroupStruct(cStructDict):
    """BGC2 group structure"""
    def __init__(self,struct_type):
        fields = None
        # GDATA_FORMAT_ID
        if struct_type == 10:
            pass
        # GDATA_FORMAT_RM
        elif struct_type == 20:
            pass
        # GDATA_FORMAT_RMPV
        elif struct_type == 30:
            pass
        # GDATA_FORMAT_RMPVMAX
        elif struct_type == 40:
            fields = [('id','q'),
                      ('parent_id','q'),
                      ('npart','Q'),
                      ('npart_self','Q'),
                      ('radius','f'),
                      ('mass','f'),
                      ('pos','fff'),
                      ('vel','fff'),
                      ('vmax','f'),
                      ('rvmax','f')]
        # Error
        if fields is None:
            raise ValueError('Unsupported group data format: {}.'.format(struct_type))
        super(Bgc2GroupStruct,self).__init__(fields)

class Bgc2PartStruct(cStructDict):
    """BGC2 particle structure"""
    def __init__(self,struct_type):
        fields = None
        # PDATA_FORMAT_NULL = 0,
        if struct_type == 0:
            pass
        # PDATA_FORMAT_ID = 10,
        elif struct_type == 10:
            pass
        # PDATA_FORMAT_IDBE = 15,
        elif struct_type == 15:
            pass
        # PDATA_FORMAT_POS = 20,
        elif struct_type == 20:
            pass
        # PDATA_FORMAT_POSBE = 25,
        elif struct_type == 25:
            pass
        # PDATA_FORMAT_PV = 30,
        elif struct_type == 30:
            fields = [('part_id','q'),
                      ('pos','fff'),
                      ('vel','fff')]
        # PDATA_FORMAT_PVBE = 35,
        elif struct_type == 35:
            pass
        # PDATA_FORMAT_PVM = 40,
        elif struct_type == 40:
            pass
        # PDATA_FORMAT_PVMBE = 45,
        elif struct_type == 45:
            pass
        # PDATA_FORMAT_GPVM = 50
        elif struct_type == 50:
            pass
        # Error
        if fields is None:
            raise ValueError('Unsupported particle data format: {}.'.format(struct_type))
        super(Bgc2PartStruct,self).__init__(fields)

# ------------------------------------------------------------------------------
# TIPSY FORMAT
@register_snapshot_format(4)
class Tipsy(Snapshot):
    """Tipsy Snapshot"""
    def read(self,*args,**kwargs):
        """See `io.read_tipsy`."""
        return read_tipsy(*args,**kwargs)
    def parse_voropar(self,param):
        """Maps Vorovol parameters onto input keywords for this class's 
        read/write methods.

        Args:
            param (dict): Vorovol parameters. This snapshot class performs the
                following parameter to keyword mappings:
                    ParticleType -> ptype (int): Particle type.

        Returns:
            dict: Keyword arguments for this class's read/write methods.

        """
        kwargs = {}
        if 'ParticleType' in param:
            kwargs['ptype'] = param['ParticleType']
        return kwargs

def read_tipsy(filename,ptype=-1,return_npart=False,return_header=False,
               dtype_pos='f',dtype_vel='f'):
    """Read a Tipsy snapshot.

    Args:
        filename (str): Full path to file that should be read.
        ptype (Optional[int]): Code specifying what particle type should be 
            loaded. Supported values include:
               -1: All particles
                0: Gas particles
                1: Dark matter particles
                2: Star particles
        return_npart (Optional[bool]): If True, only the number of particles 
            is read from the file (default = False).
        return_header (Optional[bool]): If True, only the file header is read
            from the file (default = False).
        dtype_pos (Optional[str]): Format code for particle positions. 
            (default = 'f')
        dtype_vel (Optional[str]): Format code for particle velocities.
            (default = 'f')

    Returns:
        mass (np.ndarray): (N,) Particle masses.
        pos (np.ndarray): (N,3) Particle positions.

    Raises:
        TypeError: If `ptype` is not an integer or series of integers.
        IOError: If particle numbers in header do not match.

    """
    pstructs = [TipsyPartStructGas(dtype_pos=dtype_pos,dtype_vel=dtype_vel),
                TipsyPartStructDM(dtype_pos=dtype_pos,dtype_vel=dtype_vel),
                TipsyPartStructStar(dtype_pos=dtype_pos,dtype_vel=dtype_vel)]
    # Set list of particle types
    if isinstance(ptype,int):
        if ptype in range(3):
            typelist = np.array([ptype])
        else:
            typelist = np.arange(3)
    elif isinstance(ptype,(list,tuple)):
        typelist = np.array(ptype)
    elif isinstance(ptype,np.ndarray):
        pass
    else:
        raise TypeError('Unrecognized value for ptype: {}'.format(ptype))
    # Open file
    fd = open(filename,'rb')
    # Read header
    header_struct = TipsyHeaderStruct()
    header = header_struct.read(fd)
    if header['ntot']!=sum(header['npart']):
        print 'ntot = ',header['ntot']
        print 'npart = ',header['npart']
        raise IOError('Error reading file: {}.'.format(filename))
    if return_header:
        fd.close()
        return header
    # Random 4 bytes
    mystery = struct.unpack('i',fd.read(4))
    print 'Tipsy Mystery Number = {}'.format(mystery)
    # Count particles & return if specified
    nout = 0 ; pc = np.zeros(3) ; mc = np.zeros(3)
    for t in typelist:
        pc[t] = nout
        nout += header['npart'][t]
    if return_npart:
        return nout
    # Preallocate
    pos = np.zeros((nout,3),dtype=float)
    mass = np.zeros((nout,),dtype=float)
    # Loop over particle families
    for t in range(3):
        if t in typelist:
            for i in range(header['npart'][t]):
                ipart = pstructs[t].read(fd)
                pos[pc[t],:] = ipart['pos']
                mass[pc[t]] = ipart['mass']
                pc[t]+=1
        else:
            fd.seek(pstructs[t].size*header['npart'][t],1)
    # Close file & return
    fd.close()
    return mass,pos

class TipsyHeaderStruct(cStructDict):
    """Tipsy header structure"""
    # Size should be 28
    def __init__(self):
        fields = [('time'            ,'d'), 
                  ('ntot'            ,'i'),
                  ('ndim'            ,'i'),
                  ('npart'           ,'iii')]
        super(TipsyHeaderStruct,self).__init__(fields)
class TipsyPartStructGas(cStructDict):
    """Tipsy gas particle structure"""
    def __init__(self,dtype_pos='f',dtype_vel='f'):
        fields = [('mass'  ,'f'        ),
                  ('pos'   ,3*dtype_pos),
                  ('vel'   ,3*dtype_vel),
                  ('rho'   ,'f'        ),
                  ('eps'   ,'f'        ),
                  ('metals','f'        ),
                  ('phi'   ,'f'        )]
        super(TipsyPartStructGas,self).__init__(fields)
class TipsyPartStructDM(cStructDict):
    """Tipsy dark matter particle structure"""
    def __init__(self,dtype_pos='f',dtype_vel='f'):
        fields = [('mass'  ,'f'        ),
                  ('pos'   ,3*dtype_pos),
                  ('vel'   ,3*dtype_vel),
                  ('eps'   ,'f'        ),
                  ('phi'   ,'f'        )]
        super(TipsyPartStructDM,self).__init__(fields)
class TipsyPartStructStar(cStructDict):
    """Tipsy star particle structure"""
    def __init__(self,dtype_pos='f',dtype_vel='f'):
        fields = [('mass'  ,'f'        ),
                  ('pos'   ,3*dtype_pos),
                  ('vel'   ,3*dtype_vel),
                  ('metals','f'        ),
                  ('tform' ,'f'        ),
                  ('eps'   ,'f'        ),
                  ('phi'   ,'f'        )]
        super(TipsyPartStructStar,self).__init__(fields)

# ------------------------------------------------------------------------------
# SUPPORTING FUNCTIONS
def _read_c_struct(fd,fields):
    """Read c structure from file as a dictionary.

    Args:
        fd (file): Open file object at the place where read should begin.
        fields (List[tuple]): Each element is a length 2 tuple containing the 
            field name and data type string (See `struct` documentation for 
            info on specifying different types). Fields with multiple 
            type characters are returned as arrays if the types are all the
            same and lists if the types are not.

    Returns: 
        dict: Containing struct info read from the file.

    Example:
        Assuming the file is structured in the correct format,::

            >>> fd = open('somefile.dat','rb')
            >>> fields = [('N','i'), ('pos','fff'), ('etc','if')]
            >>> x = _read_binary_dict(fd, fields)
 
        the above will return a dictionary with 3 fields::

            >>> x.keys()
            ['N','pos','etc']

        including a single integer,::

            >>> type(x['N'])
            int

        an array of three floats,::

            >>> type(x['pos'])
            numpy.ndarray
            >>> x['pos'].dtype
            dtype('float32')
            >>> len(x['pos'])
            3

        and a list containing one integer and one floats.::

            >>> type(x['etc'])
            list
            >>> len(x['etc'])
            2
            >>> type(x['etc'][0])
            int
            >>> type(x['etc'][1])
            float

    """
    import struct
    # Get format string & size
    format = '' ; keys = []
    for f in fields:
        keys.append(f[0])
        format+=f[1]
    size = struct.calcsize(format)
    # Get list of values
    values = struct.unpack(format,fd.read(size))
    # Create dictionary
    out = {} ; i = 0
    for f in fields:
        f_key = f[0]
        f_fmt = f[1]
        f_len = len(f_fmt)
        # Single value
        if f_len == 0:
            out[f_key] = values[i]
            i+=1
        # Multiple values
        else:
            out[f_key] = values[i:(i+f_len)]
            if len(set(f_fmt))==1:
                out[f_key] = np.array(out[f_key])
            i+= f_len
    # Return dictionary
    return out
        

