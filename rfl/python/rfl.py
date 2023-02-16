"""
Response Function Library
Primary Author: Paul O'Brien, The Aerospace Corporation paul.obrien@aero.org

A note on units:
    The standard allows for L_UNIT (length unit) and E_UNIT (energy unit)
    But internally, all parameters are represented using cm and MeV
    So, units will be converted on initialization/load

glossary:
    E - energy, MeV
    alpha - local pitch angle, degrees
    beta - gyrophase angle, degrees
    theta - sensor polar angle, degrees
    phi - sensor azimuth angle, degrees
    alpha0 - local pitch angle of particle incident at theta=0, degrees
    beta0 - gyrophase angle of particle incident at theta=0, degrees
    phib - sensor azimuth angle of magnetic field vector, degrees
    Egrid - energy grid, MeV
    alphagrid - alpha grid, degrees, 1-d numpy array or scalar
    betagrid - beta grid, degrees, 1-d numpy array or scalar
    thetagrid - theta grid, degrees, 1-d numpy array or scalar
    phigrid - phi grid, degrees, 1-d numpy array or scalar
    tgrid - time grid, seconds, 1-d numpy array or scalar
    hE - energy response, units of MeV (CROSSCALIB applied)
    hA* - angular response, units of cm^2(-s) (CROSSCALIB *not* applied)
    h* - energy-angle response, units of MeV-cm^2-sr(-s) (CROSSCALIB applied)
    
    
The ChannelResponse and its internal classes EnergyResponse and AngleResponse
have implicit factory constructors via the FactoryConstructorMixin. 
That means passing the constructor a data tree that actually defines a subclass, 
then the appropriate subclass is returned. see FactoryConstructorMixin
for information on how to define the is_mine method in any subclasses
defined from ChannelResponse, EnergyResponse, and AngleResponse

"""

# NOTE: in this code where the docstring is *INHERIT* that means
# the string *INHERIT* will be replaced by the docstring for the
# parent class's corresponding method

import numpy as np
from scipy.interpolate import interp1d,RegularGridInterpolator

# utilities

# degrees versions of common trig ops
sind = lambda x: np.sin(np.radians(x))
cosd = lambda x: np.cos(np.radians(x))
tand = lambda x: np.tan(np.radians(x))
asind = lambda x: np.degrees(np.arcsin(x))
acosd = lambda x: np.degrees(np.arccos(x))
atand = lambda x: np.degrees(np.arctan(x))
atan2d = lambda x,y: np.degrees(np.arctan2(y,x))


def inherit_docstrings(cls,parent = None):
    """
    inherit_docstrings(cls,parent = None)
    propagate docstring down the class tree
    replacing string "*INHERIT*" with docstring
    from corresponding item in parent class
    """
    if parent:
        for name in dir(cls):
            if name.startswith('__'): continue # don't mess with special items
            m = getattr(cls,name)
            if hasattr(parent,'__doc__') and hasattr(m,'__doc__') and (getattr(m,'__doc__') is not None):
                if type(m).__name__ == 'method':
                    # set the method's function's docstring
                    m.__func__.__doc__ = m.__func__.__doc__.replace('*INHERIT*',parent.__doc__)
                else:
                    # set the object's docstring
                    m.__doc__ = m.__doc__.replace('*INHERIT*',parent.__doc__)
    for c in cls.__subclasses__():
        inherit_docstrings(c,parent=cls)

def get_list_neighbors(lst,q):
    """
    i = get_list_neighbors(lst,q)
    find bounding neighbors of q in list lst
    if q is scalar:
        returns i = (2,), where lst[i[0]] <= q < lst[i[1]]
    else:
        returns i = (len(q),2), where lst[i[:,0]] <= q < lst[i[:,1]]
    i == -1 implies q out of range of lst
    """
    lst = np.array(lst).ravel()
    N = len(lst)
    is_scalar = np.isscalar(q)
    q = np.atleast_1d(q).ravel()
    i = np.full((len(q),2),-1)
    iq = np.isfinite(q)
    i[iq,0] = interp1d(lst,np.arange(N,dtype=int),'nearest',bounds_error=False)(q[iq])
    f = (i[:,0]>=0)
    if any(f):
        i[f,0] = i[f,0] - (q[f]<lst[i[f,0]]) # force to greatest list <= q
        i[f,1] = i[f,0]+1 # next neighbor to right
    # limit checks
    i[i>=N] = -1
    i[i<0] = -1
    if is_scalar:
        i = i[0,:]
    return i

def make_deltas(grid,int_method='trapz',**kwargs):
    """
    d = make_deltas(grid,int_method='trapz',...)
    d are the default weights for a given 1-d grid
    d is (N,), unless N<=1, in which case d is scalar
    for empty grid or None, returns 1
    for singleton grid, returns grid (assumes it's actually dt)
    keyword args besides int_method ignored
    """
    if grid is None:
        return 1
    grid = np.atleast_1d(grid).ravel()
    if len(grid) == 0:
        return 1
    elif len(grid) == 1:
        return grid[0] # grid is actually dt
    
    d = np.full(grid.shape,np.nan)
    if int_method == 'trapz':
        d[0] = (grid[1]-grid[0])/2
        d[1:-1] = (grid[2:]-grid[:-2])/2
        d[-1] = (grid[-1]-grid[-2])/2
    else:
        raise ValueError('Unknown int_method "%s"' % int_method)
    
    return d

def broadcast_grids(*args,mesh=False):
    """
    X1,X2,... = broadcast_grids(x1,x2,...,mesh=False)
    prepare 2 or more grids for broadcasting
    if mesh==False, then the grids are output with enough singleton deimensions
      added before and after the data dimension so as to be mutually broadcastable
      e.g., X1 is (N1,1), X2 is (1,N2) in 2-d case
      e.g., X1 is (N1,1,1), X2 is (1,N2,1), and X3 is (1,1,N2) in 3-d case
    if mesh==True, then the grids are output as shape (len(x1),len(x2),len(x3))
    """
    if mesh:
        return np.meshgrid([x.ravel() for x in args],indexing='ij')
    else:
        NX = len(args)
        X = [None]*NX
        for i,x in enumerate(args):
            s = [1]*NX
            s[i] = np.array(x).size
            X[i] = np.reshape(x,s)
        return(tuple(X))

def validate_grid(grid):
    """
    returns true if grid is value:
        1-d numpy array and strictly increasing
    
    """
    if grid.shape != (grid.size,):
        return False
    if any(np.diff(grid)<0):
        return False
    return True

def default_thetagrid(**kwargs):
    """
    betagrid = thetagrid(**kwargs)
    returns default theta grid, as modified by keyword args
    """
    # TODO: implement some keyword controls
    return np.linspace(0.0,180.0,181) # 1 degree spacing, closed endpoints

def default_phigrid(**kwargs):
    """
    phigrid = default_phigrid(**kwargs)
    returns default phi grid, as modified by keyword args
    """
    # TODO: implement some keyword controls
    return np.linspace(0.0,360.0,361) # 1 degree spacing, closed endpoints

def default_alphagrid(**kwargs):
    """
    alphagrid = alphagrid(**kwargs)
    returns default alpha grid, as modified by keyword args
    """
    # TODO: implement some keyword controls
    return np.linspace(0.0,180.0,181) # 1 degree spacing, closed endpoints

def default_betagrid(**kwargs):
    """
    betagrid = betagrid(**kwargs)
    returns default beta grid, as modified by keyword args
    """
    # TODO: implement some keyword controls
    return np.linspace(0.0,360.0,361) # 1 degree spacing, closed endpoints

def alphabeta2thetaphi(alpha,beta,alpha0,beta0,phib):
    """
    (theta,phi) = alphabeta2thetaphi(alpha,beta,alpha0,beta0,phib)
    Convert pitch angle and gyrophase to instrument angles. All angles
    degrees. When theta=0,180 phi=0.
    alpha, beta, should all be mutually broadcastable
    alpha0, beta0, phib should all be mutually broadcastable
    output will be (*broadcast(alpha,beta).shape,*broadcast(alpha0,beta0,phib).shape)
    """
    
    try:
        _ = np.broadcast(alpha,beta)
    except ValueError:
        raise ValueError('alpha, and beta must be mutually broadcastable')
    try:
        _ = np.broadcast(alpha0,beta0,phib)
    except ValueError:
        raise ValueError('alpha0, beta0, and phib must be mutually broadcastable')

        
    # define rotation matrix
    # columns are coefficients of c,d,b terms of inst basis vectors
    # rows are coefficients of s1,s2,s0 terms of mag basis vectors
    # note: python inner lists are columns, in matlab row of code is row of matrix, so this code is transposed from matlab
    R = np.array([
        [-sind(beta0)*sind(phib)-cosd(alpha0)*cosd(beta0)*cosd(phib), cosd(beta0)*sind(phib)- cosd(alpha0)*cosd(phib)*sind(beta0),cosd(phib)*sind(alpha0)]
        [cosd(phib)*sind(beta0)-cosd(alpha0)*cosd(beta0)*sind(phib), -cosd(beta0)*cosd(phib)-cosd(alpha0)*sind(beta0)*sind(phib),sind(alpha0)*sind(phib)]
        [cosd(beta0)*sind(alpha0), -cosd(beta0)*cosd(phib)-cosd(alpha0)*sind(beta0)*sind(phib), sind(alpha0)*sind(beta0),cosd(alpha0)]
        ])


    x = np.array([sind(alpha)*cosd(beta),sind(alpha)*sind(beta),cosd(alpha)])
    # y = R*x ( sum over 2nd dim of R and 1st dim of x )
    y = np.tensordot(x,R,axies=((0,),(1,))) # do it in this order to get alpha,beta shape before alpha0,beta0,phib shape
    theta = acosd(y[3])
    phi = atan2d(y[2],y[1])
    phi[(theta==0.0)|(theta==180.0)] = 0
    return theta,phi

def thetaphi2alphabeta(theta,phi,alpha0,beta0,phib):
    """
    (alpha,beta) = thetaphi2alphabeta(theta,phi,alpha0,beta0,phib)
    convert instrument angles to pitch angle and gyrophase. All angles
    degrees. When alpha=0,180 beta=0.
    theta,phi, should all be mutually broadcastable
    alpha0, beta0, phib should all be mutually broadcastable
    output will be (*broadcast(theta,phi).shape,*broadcast(alpha0,beta0,phib).shape)
    """
    # the alpha,beta -> theta,phi transform is identical under
    # alpha <-> theta
    # beta <-> phi
    # beta0 <-> phib (makes R <-> R')
    # so just call alphabeta2thetaphi with phib,beta0 swapped

    return alphabeta2thetaphi(theta,phi,alpha0,phib,beta0) # phib <-> beta0


# classes


class RFLError(Exception):
    """
    Unknown RFL error
    accepts optional message argument
    """
    def __init__(self,message = 'Unknown RFL Error'):
        super().__init__(message)

class ArgSizeError(ValueError):
    """
    Incorrect/Inconsistent Argument Size
    accepts optional message argument
    """
    def __init__(self,message = 'Incorrect/Inconsistent Argument Size'):
        super().__init__(message)

class KeywordError(RFLError):
    """
    Missing/Invalid Keyword
    accepts optional message argument
    """
    def __init__(self,message = 'Missing/Invalid Keyword'):
        super().__init__(message)

# RFL helper functions
def keyword_check(*keywords,**kwargs):
    """
    bool = keyword_check(*keywords,**kwargs)
    check one or more keywords in kwargs
    if present and str, return True
    otherwise raise appropriate exception
    """
    for key in keywords:
        if key not in kwargs: raise KeywordError('Keyword %s required, absent' % key)
        if not isinstance(kwargs[key],str): raise KeywordError('%s must be a string(str)' % key)
    return True

def keyword_check_bool(*keywords,**kwargs):
        """
        bool = keyword_check_bool(*keywords,**kwargs)
        check one or more keywords in kwargs
        if present and str, return True
        otherwise returns False
        """
        for key in keywords:
            if key not in kwargs: return False
            if not isinstance(kwargs[key],str): return False
        return True
    
def keyword_check_numeric(*keywords,**kwargs):
    """
    bool = keyword_check_numeric(*keywords,**kwargs)
    check one or more keywords in kwargs
    if present and np.numeric, return True
    otherwise raise appropriate exception
    """
    for key in keywords:
        if key not in kwargs: raise KeywordError('Keyword %s required, absent' % key)
        val = kwargs[key]
        while isinstance(val,(list,np.ndarray)):
            val = val[0]
        if not np.issubdtype(type(val),np.number): raise KeywordError('%s must be numeric' % key)
    return True

def keyword_check_numeric_bool(*keywords,**kwargs):
    """
    bool = keyword_check_numeric_bool(*keywords,**kwargs)
    check one or more keywords in kwargs
    if present and np.numeric, return True
    otherwise return False
    """
    for key in keywords:
        if key not in kwargs: return False
        val = kwargs[key]
        while isinstance(val,(list,np.ndarray)):
            val = val[0]
        if not np.issubdtype(type(val),np.number): return False
    return True

def squeeze(val):
    """
    val = squeeze(val)
    squeeze a value by removing trailing singleton dimensions
    and if a single value, make it a scalar
    """
    val = np.squeeze(val)
    if val.size==1:
        val = val.item() # reduce singleton to scalar (replaces .asscalar)
    return val

def interp_weights_1d(xgrid,xhat,extrap_left=False,extrap_right=False,period=None):
    """
        v = interp_weights_1d(xgrid,xhat,extrap_left=False,extrap_right=False,period=None)
        xgrid is list of grid x values (Nx,)
        xhat is list of query x values (N,) or scalar
        extrap_left - 
            0,False - zero for xhat < xgrid[0]
            True - linearly extrapolate for xhat < xgrid[0]
            'fixed' xhat < xgrid[0] use value at xgrid[0]
        extrap_right - 
            0,False - zero for xhat > xgrid[-1]
            True - linearly extrapolate for xhat > xgrid[-1]
            'fixed' xhat > xgrid[-1] use value at xgrid[-1]
        period - specify period over which x values wrap around
        v is the value of the interpolation weights from the xgrid to points xhat
        v is (N,Nx)
        if xhat is a scalar, then v is (Nx,)
    """
    xgrid = np.array(xgrid).ravel()
    isscalar = np.isscalar(xhat)
    xhat = np.atleast_1d(xhat).ravel()
    if period:
        xhat = np.remainder(xhat,period)
    Nx = xgrid.size
    if Nx < 2:
        raise ValueError('xgrid must have at least 2 entries')
    N = xhat.size
    v = np.zeros((N,Nx))
    I = get_list_neighbors(xgrid,xhat) # (N,2)
    # prepare to broadacst
    
    for (k,(ileft,iright)) in enumerate(I):
        if (ileft>=0) & (ileft<Nx-1): #xgrid[i] <= x < xgrid[i+1]
            v[k,ileft] = (xgrid[ileft+1]-xhat[k])/(xgrid[ileft+1]-xgrid[ileft])
        if (iright>0) & (iright<Nx): # xgrid[i-1] <= x < xgrid[i]
            v[k,iright] = (xhat[k]-xgrid[iright-1])/(xgrid[iright]-xgrid[iright-1])
            

    if extrap_left == 'fixed':
        f = xhat < xgrid[0]
        v[f,0] = 1.0
    elif extrap_left:
        dx = (xgrid[1]-xgrid[0])
        f = (xhat < xgrid[0]) & (xhat > xgrid[0]-dx)
        v[f,0] = 1.0-(xgrid[0]-xhat[f])/dx # extrapolate left
        

    if period:
        dx = xgrid[0]+period-xgrid[-1]
        f = xhat<xgrid[0]
        v[f,-1] = (xgrid[0]-xhat[f])/dx
        v[f,0] = 1-v[f,-1]
        f = xhat>xgrid[-1]
        v[f,0] = (xhat[f]-xgrid[-1])/dx
        v[f,-1] = 1-v[f,0]
    elif extrap_right == 'fixed':
        f = xhat > xgrid[-1]
        v[f,-1] = 1.0
    elif extrap_right:
        dx = xgrid[-1]-xgrid[-2]
        f = (xhat >= xgrid[-1]) & (xhat < xgrid[-1]+dx)
        v[f,-1] = 1.0-(xhat[f]-xgrid[-1])/dx
    if isscalar:
        v = v[0,:] # (1,Nx) -> (Nx,)
    return v

def interp_weights_3d(xgrid,xhat,ygrid,yhat,zgrid=None,zhat=None,xopts={},yopts={},zopts={}):
    """
        v = interp_weights_2d(xgrid,xhat,ygrid,yhat,zgrid=None,zhat=None,xopts={},yopts={},zopts={})
        xgrid is list of grid x values (Nx,)
        xhat is list of query x values (N,) or scalar
        ygrid is list of grid y values (Ny,)
        yhat is list of query y values (N,) or scalar
        zgrid is list of grid y values (Nz,)
        zhat is list of query z values (N,) or scalar
        xopts is a dict of options for interp_weigths_1d for x
        yopts is a dict of options for interp_weigths_1d for y
        zopts is a dict of options for interp_weigths_1d for z
        v is the value of the interpolation weights from the 
          xgrid x ygrid x zgrid to points xhat,yhat,zhat
        v is (N,Nx,Ny,Nz)
        if xhat, yhat and zhat are scalar, then v will be (Nx,Ny,Nz)
        if zgrid and zhat are None, produce 2-d interpolation weights 
          without Nz axis on output
    """
    if (zgrid is None) != (zhat is None):
        raise ValueError('if either zgrid or zhat is None, both must be')

    isscalar = np.isscalar(xhat) and np.isscalar(yhat)
    xhat = np.atleast_1d(xhat).ravel()
    yhat = np.atleast_1d(yhat).ravel()
    N = max(len(xhat),len(yhat),len(zhat))
    if zhat is not None:
        isscalar = isscalar and np.isscalar(zhat)
        zhat = np.atleast_1d(zhat).ravel()
        N = max(N,len(zhat))

    if len(xhat) not in [1,N]:
        raise ValueError('xhat must be broadcastable to same size as yhat,zhat')
    wx = interp_weights_1d(xgrid,xhat,**xopts) # len(xhat) x Nx
    if len(yhat) not in [1,N]:
        raise ValueError('yhat must be broadcastable to same size as xhat,zhat')
    wy = interp_weights_1d(ygrid,yhat,**yopts) # len(yhat) x Ny
    if zhat is not None:
        if len(zhat) not in [1,N]:
            N = max(N,len(zhat))
            raise ValueError('zhat must be broadcastable to same size as xhat,yhat')
        wz = interp_weights_1d(zgrid,zhat,**zopts) # len(zhat) x Nz

    if zhat is not None: # make broadcastable to (N,Nx,Ny,Nz)
        wy = np.expand_dims(np.expand_dims(wy,1),3) 
        wz = np.expand_dims(np.expand_dims(wz,1),1)
    else: # make broadcastable to (N,Nx,Ny)
        wy = np.expand_dims(wy,1)
        wz = 1 # dummy

    w = wx*wy*wz
    if isscalar: # remove first dim
        w = np.squeeze(w,axis=0)
    return w

def interp_weights_2d(xgrid,xhat,ygrid,yhat,xopts={},yopts={}):
    """
        v = interp_weights_2d(xgrid,xhat,ygrid,yhat,xopts={},yopts={})
        xgrid is list of grid x values (Nx,)
        xhat is list of query x values (N,) or scalar
        ygrid is list of grid y values (Ny,)
        yhat is list of query y values (N,) or scalar
        xopts is a dict of options for interp_weigths_1d for x
        yopts is a dict of options for interp_weigths_1d for y
        v is the value of the interpolation weights from the 
          xgrid x ygrid to points xhat,yhat
        v is (N,Nx,Ny)
        if xhat and yhat are scalar, then v will be (Nx,Ny)
    """
    return interp_weights_3d(xgrid,xhat,ygrid,yhat,xopts,yopts)

class FactoryConstructorMixin(object):
    """
    class FactoryConstructorMixin
    a mixin to support letting an object constructor
    traverse its subclass tree to instantiate the correct
    subclass based on the initialization arguments

    Class trees that inherit from FactoryConstructorMixin
    should implement is_mine in every subclass to "claim" the
    initialiation data and insantiate that subclass rather than
    the parent constructor that was called.    

    @classmethod
    def is_mine(cls,*args,**kwargs):
        # search args,kwargs dict/tree to determine if this subclass should handle
        # the data, and if so, return True
    
    Any class that inherites directly from FactoryConstructorMixin
    should implement is_mine that raises an exception rather than returning False
    
    An subclass that does not define is_mine will replace its parent in the
    factory hierarchy (use this approach to modify the behavior of a subclass).
    
    Any subclass that defines its own is_mine will extend its parent into a
    distinct, new subclass in the factory hierarchy (use this approach if
    parent and child are alternatives that should co-exist)
    
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        """
            bool = is_mine(cls,*args,**kwargs)
            returns true if args, kwargs describe an object of this class
        """
        raise Exception("Class '%s' failed to implement is_mine classmethod or didn't claim case" % cls.__name__)
    def __new__(cls,*args,**kwargs):
        # check children first
        if args or kwargs:
            for c in cls.__subclasses__():
                res = c.__new__(c,*args,**kwargs)
                if res: return res
        if cls.is_mine(*args,**kwargs):
            return super().__new__(cls) # python will pass kwargs to init
        return None            
    def __init__(self,*args,**kwargs):
        super().__init__() # avoid trying to call object initializer with arguments
 

class EnergyResponse(FactoryConstructorMixin):
    """
    EnergyResponse
    abstract class representing energy responses
    initialization accepts CROSSCALIB and EPS (default to 1 if absent)
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        """*INHERIT*"""
        keyword_check('E_TYPE',**kwargs)
        raise KeywordError('Supplied keywords to not define a recognized EnergyResponse')
    
    def __init__(self,**kwargs):
        for arg in ['CROSSCALIB','EPS']:
            if keyword_check_numeric_bool(arg):
                setattr(self,arg,squeeze(kwargs[arg]))
            else:
                print('%s not found or not numeric in channel description, assuming unity' % arg)
                setattr(self,arg,1.0) # default to unity
        
    def RE(self,E):
        """
        R = RE(E)
        Channel response at energy E (w/o CROSSCALIB)
        output is dimensionless
        """
        raise NotImplementedError # abstract base class
    def hE(self,Egrid,**kwargs):
        """
        hE = .hE(Egrid,...)
        return energy response weights on Egrid
        Egrid and hE are 1-d numpy arrays of the sampe size
        CROSSCALIB applied
        """
        raise NotImplementedError  # abstract base class
    def hE0(self,Egrid=None,**kwargs):
        """
        hE = .hE(Egrid=None,...)
        return energy response integrated over Egrid
        if Egrid is not provided, attempts to provide its own
        hE is a scalar, MeV
        CROSSCALIB applied
        """
        if Egrid is None:
            if self._E0 is None:
                raise ValueError('Egrid must be supplied for %s.hE0' % self.__class__.__name__)
            return self._E0
        else:
            hE = self.hE(Egrid,**kwargs)
            return hE.sum() # integrate over energy

class ER_Diff(EnergyResponse):
    """
    ER_Diff
    class representing differential energy channel response
    initialization requires EPS, E0, DE
    EPS is optional (defaults to 1)
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        """*INHERIT*"""
        return ('E_TYPE' in kwargs) and (kwargs['E_TYPE'] == 'DIFF')

    def __init__(self,**kwargs):
        super().__init__(**kwargs) # handles CROSSCALIB and EPS
        keyword_check_numeric('E0','DE',**kwargs)
        self.E0 = squeeze(kwargs['E0']) # TODO - convert to MeV
        self.DE = squeeze(kwargs['DE']) # TODO - convert to MeV
        self._hE0 = self.DE*self.EPS/self.CROSSCALIB # hE for flat spectrum

    def RE(self,E):
        """*INHERIT*"""
        return (E==self.E0)*self.EPS # DE gets ignored, which is bad, but unaviodable

    def hE(self,Egrid,**kwargs):
        """*INHERIT*"""
        hE = interp_weights_1d(Egrid,self.E0)*self.dE*self.EPS
        return hE / self.CROSSCALIB
        
class ER_Int(EnergyResponse):
    """
    ER_Int
    class representing integral energy channel response
    initialization requires EPS, E0
    propertis:
        hE0 - integrated energy response for flat spectrum, MeV
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        """*INHERIT*"""
        return ('E_TYPE' in kwargs) and (kwargs['E_TYPE'] == 'INT')

    def __init__(self,**kwargs):
        super().__init__(**kwargs) # handles CROSSCALIB and EPS
        keyword_check_numeric('E0',**kwargs)            
        self.E0 = squeeze(kwargs['E0']) # TODO - convert to MeV
        self._hE0 = None # not defined  for integral channel

    def RE(self,E):
        """*INHERIT*"""
        return (E>=self.E0)*self.EPS

    def hE(self,Egrid,**kwargs):
        """*INHERIT*"""
        Egrid = np.array(Egrid)
        NE = len(Egrid)
        hE = np.zeros(Egrid.shape) 
        I = get_list_neighbors(Egrid,self.E0)
        # Egrid[I[0]] <= inst_info.E0 < Egrid[I[1]]
        # or I[1] = -1 or both I=-1
        
        dE = make_deltas(Egrid,**kwargs)
        hE = dE # default, eventually only kept for Egrid>Egrid[I[1]]
        hE[0:(I[0]+1)] = 0 # zero out below E0
        
        # left side
        i = I[0]
        if (i>=0) and (i<NE-1): # E[i] <= E0 < E[i+1]
            hE[i] = (Egrid[i+1]-self.E0)**2/2/(Egrid[i+1]-Egrid[i])
        
        # right side
        i = I[1]
        if (i>=1) and (i<NE): # E[i-1] < E0 < E[i]
            hE[i] = (self.E0-Egrid[i-1])**2/2/(Egrid[i]-Egrid[i-1])
            
        return hE *self.EPS / self.CROSSCALIB
    
class ER_Wide(EnergyResponse):
    """
    ER_Wide
    class representing wide differential energy channel response
    initialization requires EPS, E0, E1
    EPS is optional (defaults to 1)
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        """*INHERIT*"""
        return ('E_TYPE' in kwargs) and (kwargs['E_TYPE'] == 'WIDE')

    def __init__(self,**kwargs):
        super().__init__(**kwargs) # handles CROSSCALIB and EPS
        keyword_check_numeric('E0','E1',**kwargs)            
        self.E0 = squeeze(kwargs['E0']) # TODO - convert to MeV
        self.E1 = squeeze(kwargs['E1']) # TODO - convert to MeV
        if self.E1<=self.E0:
            raise ValueError('E1 must be greater than E0 for wide channel')
        self._hE0 = (self.E1-self.E0)*self.EPS/self.CROSSCALIB  # for flat spectrum
        # build two integral channels to difference them
        tmp = {**kwargs,'E_TYPE':'INT'} # integral channel w/ same E0
        self.low = ER_Int(**tmp)
        tmp['E0'] = self.E1 # integral channel above this one
        self.high = ER_Int(**tmp)

    def RE(self,E):
        """*INHERIT*"""
        return ((E>=self.E0) & (E<=self.E1))*self.EPS

    def hE(self,Egrid,**kwargs):
        """*INHERIT*"""
        # treat as difference between two integral channels
        return self.low.hE(Egrid)-self.high.hE(Egrid)
        
class ER_Table(EnergyResponse):
    """
    ER_Table
    class representing tabular energy channel response
    initialization requires EPS, E_GRID
    EPS is optional (defaults to 1)
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        """*INHERIT*"""
        return ('E_TYPE' in kwargs) and (kwargs['E_TYPE'] == 'TBL')

    def __init__(self,**kwargs):
        super().__init__(**kwargs) # handles CROSSCALIB and EPS
        keyword_check_numeric('E_GRID',**kwargs)
        self.E_GRID = squeeze(kwargs['E_GRID']) # TODO convert to MeV
        if not validate_grid(self.E_GRID):
            raise ValueError('E_GRID is not a valid grid: 1-d, unique')

        if self.EPS.shape != self.E_GRID.shape:
            raise ArgSizeError('E_GRID and EPS are not the same shape')
            
        self._hE0 = np.sum(self.E_GRID*self.EPS*make_deltas(self.E_GRID))/self.CROSSCALIB  # for flat spectrum
        self._RE = None

    def RE(self,E):
        """*INHERIT*"""
        if self._RE is None:
            # extrapolate beyond grid with nearest edge (first/last) value
            self._RE = interp1d(self.E_GRID,self.EPS,'linear',bounds_error=False,fill_value=[self.EPS[0],self.EPS[-1]],assum_sorted=True)
        return self._RE(E)

    def hE(self,Egrid,**kwargs):
        """*INHERIT*"""
        
        dE = make_deltas(Egrid,**kwargs)
        hE = self.RE(Egrid)*dE/self.CROSSCALIB
        return hE

class AngleResponse(FactoryConstructorMixin):
    """
    AngleResponse
    class representing angular responses
    all responses are forward-only and
    add backward response via self.backward when bidirectional
    initialization accepts BIDIRECTIONAL (False if missing)
    
    properties:
        bidirectional - bool, is the angular response bidirectional?
        backward - backward angular response (null class when not bidirectional)
        hA0 - goemetric factor for forward and (if bidirectional) backward hemispheres
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        """*INHERIT*"""     
        if not kwargs:
            return True # null response
        else:
            raise KeywordError('The data provided did not define a recognized AngleResponse')
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.bidirectional = ('BIDIRECTIONAL' in kwargs) and (kwargs['BIDIRECTIONAL'] in [True,'TRUE'])
    def A(self,theta,phi):
        """
        A = .A(theta,phi)
        return angular response (effective area) at specified theta,phi cm^2
        should add forward and backward response together
        theta and phi must broadcast together
        A is mutual broadcast shape of theta,phi
        see module glossary for input meanings
        """
        raise NotImplementedError # abstract base class
    def hAthetaphi(self,thetagrid,phigrid,**kwargs):
        """
        hAthetaphi = .hAthetaphi(thetagrid,phigrid,...)
        return angular response weights on theta x phi grid, cm^2 sr
        hAthetaphi is a scalar or shape (len(thetagrid),len(phigrid))
        see module glossary for input meanings
        """
        dcostheta = np.abs(make_deltas(-cosd(thetagrid),**kwargs))
        dphi = make_deltas(np.radians(phigrid),**kwargs)
        thetaI,phiI = broadcast_grids(thetagrid,phigrid)
        dcosthetaI,dphiI = broadcast_grids(dcostheta,dphi)
        # includes backward response if present
        return self.A(thetaI,phiI)*dcosthetaI*dphiI
    def hAtheta(self,thetagrid,phigrid=None,**kwargs):
        """
        hAtheta = .hAtheta(thetagrid,phigrid,...)
        return angular response weights on theta grid, cm^2 sr
        if phigrid is None, integrated over 0-360
        hAtheta is a scalar or 1-d numpy array of same length as thetagrid
        see module glossary for input meanings
        """
        if phigrid is None:
            phigrid = default_phigrid(**kwargs)
        h = self.hAthetaphi(thetagrid,phigrid,**kwargs)
        return h.sum(axis=1) # integrate over beta
    def hAalphabeta(self,alpha0,beta0,phib,alphagrid,betagrid,tgrid=None,**kwargs):
        """
        hAalphabeta = .hAalphabeta(alpha0,beta0,phib,alphagrid,betagrid,tgrid=None,...)
        return angular response weights on alpha x beta grid
        see module glossary for input meanings
        if tgrid is supplied, result is integrated over time grid
            expects alpha0,beta0,phi0 to depend on tgird
        output units are cm^2-sr or cm^2-sr-s if tgrid supplied
        hAalphabeta is scalar or shape (len(alphagrid),len(betagrid))
        """
        dcosa = make_deltas(-cosd(alphagrid),**kwargs) # (Na,)
        db = np.radians(make_deltas(betagrid,**kwargs)) # (Nb,)
        dt = make_deltas(tgrid,**kwargs); #(Nt,) or scalar
        dcosa,db = broadcast_grids(dcosa,db) # (Na,Nb)
        theta,phi = alphabeta2thetaphi(alphagrid,betagrid,alpha0,beta0,phib) # (Na,Nb,Nt) or (NA,Nb)
        A = self.A(E,theta,phi) # (Na,Nb,*Nt)
        if np.isscalar(dt):
            h = A*dt
        else:
            h = np.tensordot(A,dt,axes=((2,),(0,))) # dot product along Nt dimension
        # now h is (Na,Nb)
        h = h*dcosa*db/self.CROSSCALIB
        return h
    def hAalpha(self,alpha0,beta0,phib,alphagrid,betagrid=None,tgrid=None,**kwargs):
        """
        hAalpha = .hAalpha(alpha0,beta0,phib,alphagrid,betagrid=None,tgrid=None,...)
        return angular response weights on alpha grid (integrated over beta)
        see module glossary for input meanings
        if betagrid is None, integrated over 0-360
        if tgrid is supplied, result is integrated over time grid
            expects alpha0,beta0,phi0 to depend on tgird
        output units are cm^2-sr or cm^2-sr-s if tgrid supplied
        hAalpha is a scalar or 1-d numpy array of same length as alphagrid
        """
        if betagrid is None:
            betagrid = default_betagrid(**kwargs)
        h = self.hAalphabeta(alpha0,beta0,phib,alphagrid,betagrid,tgrid,**kwargs)
        return h.sum(axis=1) # integrate over beta
    @property
    def hA0(self):
        """
        return the nominal geometric factor in cm^2 sr
        including bidirectional effects
        """
        hA0 = self.G
        if self.bidirectional:
            hA0 += self.backward.hA0
        return hA0

class AR_csym(AngleResponse):
    """
    AR_csym 
    virtual base class representing phi-symmetric angular responses
    forward response only
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        """*INHERIT*"""
        return False # abstract class, never the right answer
                
    def hAtheta(self,thetagrid,phigrid=None,**kwargs):
        """*INHERIT*"""
        dcostheta = np.abs(make_deltas(-cosd(thetagrid),**kwargs))
        if phigrid is None:
            dphi = 2*np.phi
        else:
            dphi = make_deltas(np.radians(phigrid),**kwargs).sum()
        phigrid = 0.0 # dummy value to save calculation
        # includes backward response if present
        return self.A(thetagrid,phigrid)*dcostheta*dphi

class AR_Omni(AR_csym):
    """
    AR_Omni omnidirectional response
    class representing omnidirectional angular responses
    initialization requires G (None OK for derived classes)
    Note: an omni is a sphere (same apparent area from all directions). 
    For flat, single-element detectors, use DISK and SLAB
    properties:
        area - detector area, cm^2
        G - geometric factor for forward hemisphere particles (theta<=90) cm^2 sr
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        """*INHERIT*"""
        return keyword_check_bool('RESP_TYPE',**kwargs) and \
            (kwargs['RESP_TYPE'] == '[E]')
    def __init__(self,**kwargs):
        super().__init__(**kwargs) # handles BIDIRECTIONAL
        G = None
        if ('G' in kwargs) and kwargs['G'] is not None:
            # derived class may provide G as none
            keyword_check_numeric('G',**kwargs)
            # setup for half omni: full G over 2pi
            G = squeeze(kwargs['G']) # TODO: convert to cm^2 sr       
            self.G = G
            self.area = self.G/2/np.pi # equivalent sphere's cross-section 
        if self.bidirectional:
            # backwards is also omni
            backw = kwargs.copy()
            backw['BIDIRECTIONAL'] = 'FALSE'
            # divide forward part in half:
            if G is not None:
                self.G /= 2
                self.aera /= 2
                backw['G'] = self.G # already divided by 2
                backw['L_UNIT'] = 'cm'
            self.backward = self.__class__(**backw)
    def A(self,theta,phi):
        """*INHERIT*"""
        A = self.area*(theta<=90)
        if self.bidirectiona:
            A += self.backward.area*(theta>90)
        return  A

class AR_SingleElement(AR_Omni):
    """
    AR_SingleElement
    abstract class representing single-element detectors (i.e., w/ cos projection effect)
    derived class must define 
      self.area, the detector area
      self.G, the geometric factor
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        "*INHERIT*"
        return False # abstract class, never the right answer
    def __init__(self,**kwargs):
        # supply fake G to super, area and G will be set in derived class
        super().__init__(G=None,**kwargs) # handles BIDIRECTIONAL
    def A(self,theta,phi):
        """*INHERIT*"""
        return super().A(theta,phi)*np.abs(cosd(theta))
    @property
    def G(self):
        """
        Geometric factor for forward going particles
        computed as pi*area
        """
        return self.area*np.pi
class AR_Disk(AR_SingleElement):
    """
    AR_Disk
    class representing disk angular responses
    initialization requires R1, R2, D
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        return keyword_check_bool('TH_TYPE',**kwargs) and \
            (kwargs['TH_TYPE'] == 'DISK')
    def __init__(self,**kwargs):
        super().__init__(**kwargs)  # handles BIDIRECTIONAL
        keyword_check_numeric('R1',**kwargs)
        self.R1 = squeeze(kwargs['R1']) # TODO: convert to cm
        self.area = np.pi*self.R1**2

class AR_Slab(AR_SingleElement):
    """
    AR_Slab
    class representing rectangular slab angular responses
    initialization requires W1, H1
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        return keyword_check_bool('TH_TYPE',**kwargs) and \
            (kwargs['TH_TYPE'] == 'SLAB')
    def __init__(self,**kwargs):
        super().__init__(**kwargs) # handles BIDIRECTIONAL
        keyword_check_numeric('W1','H1',**kwargs)
        self.W1 = squeeze(kwargs['W1']) # TODO: convert to cm
        self.H1 = squeeze(kwargs['H1']) # TODO: convert to cm
        self.area = self.W1*self.H1

class AR_Tele_Cyl(AR_csym):
    """
    AR_Tele_Cyl
    class representing cylindrically-symmetric telescope angular responses
    Using Sullivan's 1971 paper as updated in ~2010
    initialization requires R1, R2, D
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        return keyword_check_bool('TH_TYPE',**kwargs) and \
            (kwargs['TH_TYPE'] == 'CYL_TELE')
    def __init__(self,**kwargs):
        super().__init__(**kwargs) # handles BIDIRECTIONAL
        keyword_check_numeric('R1','R2','D',**kwargs)
        self.R1 = squeeze(kwargs['R1']) # TODO: convert to cm
        self.R2 = squeeze(kwargs['R1']) # TODO: convert to cm
        self.D = squeeze(kwargs['D']) # TODO: convert to cm
        self.Rs = min(self.R1,self.R2)
        self.thetac = atand(np.abs(self.R1-self.R2)/self.D)
        self.thetam = atand((self.R1+self.R2)/self.D)
        self.area = np.pi*self.Rs**2 # area at normal incidence
        tmp = self.R1**2+self.R2**2+self.D**2
        self.G = np.pi**2/2*(tmp-np.sqrt(tmp**2-4*self.R1**2*self.R2**2)) # eqn 8
        if self.bidirectional:
            backw = kwargs.copy()
            backw['BIDIRECTIONAL'] = 'FALSE'
            # backward telescope reverses R2 and R1
            backw['R1'] = self.R1
            backw['R2'] = self.R2
            backw['L_UNIT'] = 'cm'                
            self.backward = self.__class__(**backw)
    def A(self,theta,phi):
        """*INHERIT*"""
        # uses Sullivan's updated paper, equation 10
        A = np.zeros(theta.shape)
        thetarads = np.radians(theta)
        f = theta <= self.thetac
        if any(f):
            A[f] = np.pi*self.Rs**2*np.cos(thetarads[f])
        f = (theta>self.thetac) & (theta < self.thetam)
        if any(f):
            tantheta = self.tan(thetarads[f])
            Psi1 = np.arccos((self.R1**2+self.D**2*tantheta**2-self.R2**2)/(2*self.D*self.R1*tantheta)) # rads
            Psi2 = np.arccos((self.R2**2+self.D**2*tantheta**2-self.R1**2)/(2*self.D*self.R2*tantheta)) # rads
            A[f] = np.cos(thetarads[f])/2*(self.R1**2*(2*Psi1-np.sin(2*Psi1))+self.R2**2*(2*Psi2-np.sin(2*Psi2)))
        if self.bidirectional:
            A += self.backward.A(180.0-theta,phi)
        return np.broadcast_to(A,np.broadcast(theta,phi).shape)

class AR_Table_sym(AngleResponse):
    """
    AR_Table_sym
    class representing phi-symmetric tabular angular responses
    initialization requires TH_GRID, A
    Ignores bidirectional: For bidirectional tabluar response, 
        provide A on TH_GRID that spans 0 to 180 degrees
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        return keyword_check_bool('TH_TYPE',**kwargs) and \
            (kwargs['TH_TYPE'] == 'TBL')
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        keyword_check_numeric('TH_GRID','A',**kwargs)
        self.TH_GRID = squeeze(kwargs['TH_GRID'])
        if not validate_grid(self.TH_GRID):
            raise ValueError('TH_GRID is not a valid grid: 1-d, unique')
        self._A = squeeze(kwargs['A']) # TODO: convert to cm^2
        self._Ainterpolator = None
        if self.TH_GRID.shape != self.A.shape:
            raise ArgSizeError('TH_GRID and A are not the same shape')
    def A(self,theta,phi):
        """*INHERIT*"""
        if self._Ainterpolator is None:
            # extrapolate beyond grid with nearest edge with zero
            self._Ainterpolator = interp1d(self.TH_GRID,self.A,'linear',bounds_error=False,fill_value=0.0,assum_sorted=True)
        A = self._Ainterpolator(theta)
        return np.broadcast_to(A,np.broadcast(theta,phi).shape)


class AR_Pinhole(AngleResponse):
    """
    AR_Pinhole
    class representing pinhole angular responses (delta function at THETA=0)
    initialization requires G
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):            
        return keyword_check_bool('TH_TYPE',**kwargs) and \
            (kwargs['TH_TYPE'] == 'PINHOLE')
    def __init__(self,**kwargs):
        super().__init__(**kwargs) # handles BIDIRECTIONAL
        keyword_check_numeric('G',**kwargs)
        self.G = squeeze(kwargs['G']) # TODO: convert to cm^2 sr  
        self.area = 0
        if self.bidirectional:
            self.G /= 2 # divide G between forward and backward
            backw = kwargs.copy()
            backw['BIDIRECTIONAL'] = 'FALSE'                
            backw['G'] = self.G # already divided in half
            backw['L_UNIT'] = 'cm'
            self.backward = self.__class__(**backw)
    def A(self,theta,phi):
        """*INHERIT*"""
        A = self.G*(theta==0)
        if self.bidirectional:
            A += self.G*(theta==180)
        return np.broadcast_to(A,np.broadcast(theta,phi).shape)
    def hAthetaphi(self,thetagrid,phigrid,**kwargs):
        """*INHERIT*"""
        # compute interpolation weight in -cos(theta)
        # so that integral over dcos(theta) gives unity
        w = interp_weights_1d(-cosd(thetagrid),-1.0) # -cos(0)=-1
        w = np.broadcast_to(w,np.broadcast(thetagrid,phigrid).shape)
        tmp = w/w.sum()*self.G # ensure tmp sums to G
        if self.bidirectional: # special case of actually calling hAtheta
            tmp += self.backward.hAtheta(thetagrid,phigrid,**kwargs)
        return tmp
    def hAtheta(self,thetagrid,phigrid=None,**kwargs):
        """*INHERIT*"""
        # compute interpolation weight in -cos(theta)
        # so that integral over dcos(theta) gives unity
        w = interp_weights_1d(-cosd(thetagrid),-1.0) # -cos(0)=-1
        tmp = w/w.sum()*self.G # ensure tmp sums to G
        if self.bidirectional: # special case of actually calling hAtheta
            tmp += self.backward.hAtheta(180-thetagrid,phigrid,**kwargs)
        return tmp
    def hAalphabeta(self,alpha0,beta0,phib,alphagrid,betagrid,tgrid=None,**kwargs):
        """*INHERIT*"""
        dt = np.atleast_1d(make_deltas(tgrid,**kwargs))
        alpha0 = np.broadast_to(alpha0,dt.shape)
        beta0 = np.broadast_to(beta0,dt.shape)
        # compute interpolation weight in -cos(alpha)
        # so that integral over dcos(alpha) gives unity
        dha = interp_weights_1d(-cosd(alphagrid),-cosd(alpha0)) # len(tgrid) x len(alphagrid)
        dhb = interp_weights_1d(betagrid,beta0) # len(tgrid) x len(betagrid)
        # reshape to allow broadcast
        dha.shape = (len(dt),len(alphagrid),1)
        dhb.shape = (len(dt),1,len(betagrid))
        dt.reshape = (len(dt),1,1)
        h = (dt*dha*dhb).sum(axis=0) # integrate over time
        h = h*self.G*sum(dt)/sum(h) # force to integrate to Gdt
        if self.bidirectional:
            h += self.backward.hAalpha(alpha0,beta0,phib,alphagrid,betagrid,tgrid,**kwargs)
        return h
    def hAalpha(self,alpha0,beta0,phib,alphagrid,betagrid=None,tgrid=None,**kwargs):
        """*INHERIT*"""
        dt = np.atleast_1d(make_deltas(tgrid,**kwargs))
        alpha0 = np.broadast_to(alpha0,dt.shape)
        # compute interpolation weight in -cos(alpha)
        # so that integral over dcos(alpha) gives unity
        dh = interp_weights_1d(-cosd(alphagrid),-cosd(alpha0)) # len(tgrid) x len(alphagrid)
        h = np.dot(dt,dh) # sum over time
        h = h*self.G*sum(dt)/sum(h) # force to integrate to Gdt
        if self.bidirectional:
            h += self.backward.hAalpha(alpha0,beta0,phib,alphagrid,betagrid,tgrid,**kwargs)
        return h

class AR_Tele_Rect(AngleResponse):
    """
    AR_Tele_Rect
    class representing telescope with rectangular elements
    initialization requires W1, H1, W2, H2, D
    Using Sullivan's 1971 paper as updated in ~2010
    mainly based on equation (13)
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        return keyword_check_bool('TP_TYPE',**kwargs) and \
            (kwargs['TP_TYPE'] == 'RECT_TELE')
    def __init__(self,**kwargs):
        super().__init__(**kwargs) # handles BIDIRECTIONAL
        keyword_check_numeric('W1','H1','W2','H2',**kwargs)
        self.W1 = squeeze(kwargs['W1']) # TODO: convert to cm
        self.H1 = squeeze(kwargs['H1']) # TODO: convert to cm
        self.W2 = squeeze(kwargs['W2']) # TODO: convert to cm
        self.H2 = squeeze(kwargs['H2']) # TODO: convert to cm
        self.D = squeeze(kwargs['D']) # TODO: convert to cm
        # initialize some Sullivan variables, radians
        alpha = (self.H1+self.H2)/2
        beta = (self.W1+self.W2)/2
        gamma = (self.H1-self.H2)/2
        delta = (self.W1-self.W2)/2
        # Sullivan's 11:
        self.G = self.D**2*np.log(
                ((self.D**2+alpha**2+delta**2)/(self.D**2+alpha**2+beta**2))  
                *(self.D**2+gamma**2+beta**2)/(self.D**2+gamma**2+delta**2)) \
                +2*alpha*np.sqrt(self.D**2+beta**2) \
                *np.arctan(alpha/np.sqrt(self.D**2+beta**2)) \
                +2*beta*np.sqrt(self.D**2+alpha**2) \
                *np.arctan(beta/np.sqrt(self.D**2+alpha**2)) \
                -2*alpha*np.sqrt(self.D**2+delta**2) \
                *np.arctan(alpha/np.sqrt(self.D**2+delta**2)) \
                -2*beta*np.sqrt(self.D**2+gamma**2) \
                *np.arctan(beta/np.sqrt(self.D**2+gamma**2)) \
                -2*gamma*np.sqrt(self.D**2+beta**2) \
                *np.arctan(gamma/np.sqrt(self.D**2+beta**2)) \
                -2*delta*np.sqrt(self.D**2+alpha**2) \
                *np.arctan(delta/np.sqrt(self.D**2+alpha**2)) \
                +2*gamma*np.sqrt(self.D**2+delta**2) \
                *np.arctan(gamma/np.sqrt(self.D**2+delta**2)) \
                +2*delta*np.sqrt(self.D**2+gamma**2) \
                *np.arctan(delta/np.sqrt(self.D**2+gamma**2))

        if self.bidirectional:
            backw = kwargs.copy()
            backw['BIDIRECTIONAL'] = 'FALSE'                
            backw['W1'] = self.W2
            backw['H1'] = self.H2
            backw['W2'] = self.W1
            backw['H2'] = self.H1
            backw['D'] = self.D
            backw['L_UNIT'] = 'cm'
            self.backward = self.__class__(**backw)
    def _X(self,zeta,a1,a2):
        """
        Sullivan's X(zeta,a1,a2)
        """
        return min(zeta+a1/a2,a2/2) - max(zeta-a1/a2,-a2/2)
    def A(self,theta,phi):
        """*INHERIT*"""
        A = np.zeros(np.broadcast(theta,phi).shape)
        thetarads = np.radians(theta)
        phirads = np.radians(phi)
        tantheta = np.tan(thetarads)
        zeta = self.D*tantheta*np.cos(phirads)
        eta = self.D*tantheta*np.sin(phirads)
        X = self._X(zeta,self.W1,self.W2)
        Y = self._X(eta,self.H1,self.H2)
        costheta = np.cos(thetarads)
        f = (costheta>0) & (X>0) &  (Y>0) # apply Heaviside implicitly and resolve tand(90)=inf
        if any(f):
            A[f] = costheta[f]*X[f]*Y[f] # eqn 13 w/o Heaviside functions            
        
        if self.bidirectional:
            A += self.backward.A(180-theta,phi)
        return A

class AR_Table_asym(AngleResponse):
    """
    AR_Table_asym
    class representing phi-asymmetric tabular angular responses
    initialization requires TH_GRID, PH_GRID, A
    Ignores bidirectional: For bidirectional tabluar response, 
        provide A on TH_GRID that spans 0 to 180 degrees
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        return keyword_check_bool('TP_TYPE',**kwargs) and \
            (kwargs['TP_TYPE'] == 'TBL')
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        keyword_check_numeric('TH_GRID','PH_GRID','A',**kwargs)
        for arg in ['TH_GRID','PH_GRID']:
            setattr(self,arg,squeeze(kwargs[arg]))
            if not validate_grid(getattr(self,arg)):
                raise ValueError('%s is not a valid grid: 1-d, unique' % arg)
        self._A = squeeze(kwargs['A']) # TODO: convert to cm^2
        self._Ainterpolator = None
        if self.A_GRID.shape != (len(self.TH_GRID),len(self.PH_GRID)):
            raise ArgSizeError('A does not have shape (TH_GRID,PH_GRID)')
        if (self.PH_GRID[0] !=0) or (self.PH_GRID[-1] != 360):
            raise ValueError('PH_GRID should span [0,360]')
    def A(self,theta,phi):
        """*INHERIT*"""
        if self._Ainterpolator is None:
            # extrapolate beyond grid with nearest edge with zero
            self._Ainterpolator = RegularGridInterpolator((self.TH_GRID,self.PH_GRID),self._A,'linear',bounds_error=False,fill_value=0.0)
        return self._Ainterpolator(theta,phi)
            

class ChannelResponse(FactoryConstructorMixin):
    """
    Virtual base class for channel responses
    Constructor is a factory that will return 
    appropriate subclass from response dictionary
    
    Methods:
        Working in (E,alpha,beta,time) (results integrated over time)
        hEiso - make weights for numerical integral over E (isotropic)
        hEalpha - make weights for double numerical integral over E,alpha (gyrotropic)
        hEalphabeta - make weights for triple numerical integral over E,alpha,beta
        halpha - make weights for numerical integral over alpha (flat spectrum, gyrotropic)
        halphabeta - make weights for double numerical integral over alpha,beta (flat spectrum)
        Working in (E,theta,phi)
        R - return response function itself in units of effective area, cm^2
        hE - make weights for numerical integral over E, (isotropic)
        hEtheta - make weights for double numerical integral over E,theta (ignores flux phi dependence)
        hEthetaphi - make weights for triple numerical integral over E,theta,phi
        htheta - make weights for numerical integral over theta (flat spectrum, ignores flux phi dependence)
        hthetaphi - make weights for double numerical integral over theta,phi (flat spectrum)
    """
        
    @classmethod
    def is_mine(cls,*args,**kwargs):
        """*INHERIT*"""
        keyword_check('RESP_TYPE',**kwargs)
        raise KeywordError('Supplied keywords to not define a recognized ChannelResponse')
    def copy(self):
        """
        channel.copy() return a copy of this object
        """
        return ChannelResponse(**self.to_dict())
    def to_dict(self):
        """
        channel.to_dict() return a nested dict of this object
        """
        raise NotImplementedError # TODO
    def R(self,E,theta,phi):
        """
        R = .R(E,theta,phi)
        return 3-d response (effective area * efficiency) at 
         specified E,theta,phi cm^2
        E, theta and phi must broadcast together
        R is mutual broadcast shape of E,theta,phi
        see module glossary for input meanings
        """
        raise NotImplementedError # abstract base class
    def hEthetaphi(self,Egrid,thetagrid,phigrid,**kwargs):
        """
        hEthetaphi = .hEthetaphi(Egrid,thetagrid,phigrid,...)
        return angular response weights on E x theta x phi grid, cm^2 sr MeV
        hEthetaphi is a scalar or shape (len(Egrid),len(thetagrid),len(phigrid))
        if Egrid is None, object will attempt to supply its own
        see module glossary for input meanings
        """
        if Egrid is None:
            if hasattr(self,'E_GRID'):
                Egrid = self.E_GRID
            else:
                raise ValueError('Egrid input cannot be None for object without its own E_GRID')
        dE = make_deltas(Egrid,**kwargs)
        dcostheta = np.abs(make_deltas(-cosd(thetagrid),**kwargs))
        dphi = make_deltas(np.radians(phigrid),**kwargs)
        EI,thetaI,phiI = broadcast_grids(Egrid,thetagrid,phigrid)
        dEI,dcosthetaI,dphiI = broadcast_grids(dE,dcostheta,dphi)
        return self.R(EI,thetaI,phiI)*dEI*dcosthetaI*dphiI
    def hEtheta(self,Egrid,thetagrid,phigrid=None,**kwargs):
        """
        hEtheta = .hEtheta(Egrid,thetagrid,phigrid=None,...)
        return channel response weights on theta grid, cm^2 sr MeV
        if phigrid is none, integrated over 0-360
        hEtheta is a scalar or shape (len(Egrid),len(thetagrid))
        if Egrid is None, object will attempt to supply its own
        see module glossary for input meanings
        """
        if phigrid is None:
            if hasattr(self,'PH_GRID'):
                phigrid = self.PH_GRID
            else:
                phigrid = default_phigrid(**kwargs)
        h = self.hEthetaphi(Egrid,thetagrid,phigrid,**kwargs)                
        return h.sum(axis=1) # integrate over phi
    def hE(self,Egrid,thetagrid=None,phigrid=None,**kwargs):
        """
        hEtheta = .hEtheta(Egrid,thetagrid=None,phigrid=None,...)
        return channel response weights on E grid, cm^2 sr MeV
        if thetagrid is none, integrated over 0-180
        if phigrid is none, integrated over 0-360
        hEtheta is a scalar or 1-d numpy array of same length as Egrid
        if Egrid is None, object will attempt to supply its own
        see module glossary for input meanings
        """
        if thetagrid is None:
            if hasattr(self,'TH_GRID'):
                thetagrid = self.TH_GRID
            else:
                thetagrid = default_thetagrid(**kwargs)
        h = self.hEtheta(Egrid,thetagrid,phigrid,**kwargs)                
        return h.sum(axis=1) # integrate over theta
    def hEalphabeta(self,alpha0,beta0,phib,Egrid,alphagrid,betagrid,tgrid=None,**kwargs):
        """
        hEalphabeta = .hEalphabeta(alpha0,beta0,phib,Egrid,alphagrid,betagrid,tgrid=None,...)
        return channel response weights on alpha x beta grid
        if tgrid is supplied, result is integrated over time grid
            expects alpha0,beta0,phi0 to depend on tgird
        output units are cm^2-sr-MeV or cm^2-sr-MeV-s if tgrid supplied
        if Egrid is None, object will attempt to supply its own
        hEalphabeta is scalar or shape (len(Egrid),len(alphagrid),len(betagrid))
        see module glossary for input meanings
        """
        dE = make_deltas(Egrid,**kwargs) # (NE,)
        dcosa = make_deltas(-cosd(alphagrid),**kwargs) # (Na,)
        db = np.radians(make_deltas(betagrid,**kwargs)) # (Nb,)
        dt = make_deltas(tgrid,**kwargs); #(Nt,) or scalar
        dE,dcosa,db = broadcast_grids(dE,dcosa,db) # (NE,Na,Nb)
        theta,phi = alphabeta2thetaphi(alphagrid,betagrid,alpha0,beta0,phib) # (Na,Nb,Nt) or (NA,Nb)
        E = np.reshape(Egrid,(len(Egrid),*np.ones(theta.ndim))) # (NE,1,1,*1)
        theta = np.expand_dims(theta,axis=0) # (1,Na,Nb,*Nt)
        phi = np.expand_dims(phi,axis=0) # (1,Na,Nb,*Nt)
        R = self.R(E,theta,phi) # (NE,Na,Nb,*Nt)
        if np.isscalar(dt):
            h = R*dt
        else:
            h = np.tensordot(h,dt,axes=((3,),(0,))) # dot product along Nt dimension
        # now h is (NE,Na,Nb)
        h = h*dE*dcosa*db/self.CROSSCALIB
        return h
    def hEalpha(self,alpha0,beta0,phib,Egrid,alphagrid,betagrid=None,tgrid=None,**kwargs):
        """
        hEalpha = .hEalpha(alpha0,beta0,phib,Egrid,alphagrid,betagrid=None,tgrid=None,...)
        return channel response weights on alpha grid (integrated over beta)
        if Egrid is None, object will attempt to supply its own
        if betagrid is None, supplies numeric beta grid when needed
        if tgrid is supplied, result is integrated over time grid
            expects alpha0,beta0,phi0 to depend on tgird
        output units are cm^2-sr-MeV or cm^2-sr-MeV-s if tgrid supplied
        hEalpha is a scalar or 1-d numpy array of same length as alphagrid
        see module glossary for input meanings
        """
        if betagrid is None:
            betagrid = default_betagrid(**kwargs)
        h = self.hEalphabeta(alpha0,beta0,phib,Egrid,alphagrid,betagrid,tgrid,**kwargs)
        return h.sum(axis=2) # integrate over beta
    def hEiso(self,alpha0,beta0,phib,Egrid,alphagrid,betagrid=None,tgrid=None,**kwargs):
        """
        hEiso = .hEiso(alpha0,beta0,phib,Egrid,alphagrid=None,betagrid=None,tgrid=None,...)
        return channel response weights on alpha grid (integrated over beta)
        if Egrid is None, object will attempt to supply its own
        if betagrid is None, supplies numeric beta grid when needed
        if tgrid is supplied, result is integrated over time grid
            expects alpha0,beta0,phi0 to depend on tgird
        output units are cm^2-sr-MeV or cm^2-sr-MeV-s if tgrid supplied
        hEiso is a scalar or 1-d numpy array of same length as alphagrid
        see module glossary for input meanings
        """
        if alphagrid is None:
            alphagrid = default_alphagrid(**kwargs)
        h = self.hEalpha(alpha0,beta0,phib,Egrid,alphagrid,betagrid,tgrid,**kwargs)
        return h.sum(axis=1) # integrate over alpha
    def hthetaphi(self,Egrid,thetagrid,phigrid,**kwargs):
        """
        hthetaphi = .hthetaphi(Egrid,thetagrid,phigrid,...)
        return angular response weights on theta x phi grid, cm^2 sr MeV
        integrated over E for a flat spectrum
        hthetaphi is a scalar or shape (len(thetagrid),len(phigrid))
        if Egrid is None, object will attempt to supply its own
        see module glossary for input meanings
        """
        h = self.hEthetaphi(Egrid,thetagrid,phigrid,**kwargs)
        return h.sum(axis=0) # integrate over E
    def htheta(self,Egrid,thetagrid,phigrid=None,**kwargs):
        """
        htheta = .htheta(Egrid,thetagrid,phigrid,...)
        return channel response weights on theta grid, cm^2 sr MeV
        integrated over E for a flat spectrum
        integrated over phigrid (or 0-360 if None)
        htheta is a scalar or 1-d numpy array of same length as thetagrid
        if Egrid is None, object will attempt to supply its own
        see module glossary for input meanings
        """
        h = self.hEtheta(Egrid,thetagrid,phigrid,**kwargs)
        return h.sum(axis=0)  # integrate over E
    def halphabeta(self,alpha0,beta0,phib,Egrid,alphagrid,betagrid,tgrid=None,**kwargs):
        """
        halphabeta = .halphabeta(alpha0,beta0,phib,Egrid,alphagrid,betagrid,tgrid=None,...)
        return channel response weights on alpha x beta grid
        integrated over E for a flat spectrum
        if tgrid is supplied, result is integrated over time grid
            expects alpha0,beta0,phi0 to depend on tgird
        output units are cm^2-sr-MeV or cm^2-sr-MeV-s if tgrid supplied
        if Egrid is None, object will attempt to supply its own
        halphabeta is scalar or shape (len(alphagrid),len(betagrid))
        see module glossary for input meanings
        """
        h = self.hEalphabeta(alpha0,beta0,phib,Egrid,alphagrid,betagrid,tgrid,**kwargs)
        return h.sum(axis=0) # integrate over E
    def halpha(self,alpha0,beta0,phib,Egrid,alphagrid,betagrid=None,tgrid=None,**kwargs):
        """
        hAalpha = .hAalpha(alpha0,beta0,phib,Egrid,alphagrid,betagrid=None,tgrid=None,...)
        return channel response weights on alpha grid (integrated over E,beta)
        if Egrid is None, object will attempt to supply its own
        if betagrid is None, supplies numeric beta grid when needed
        if tgrid is supplied, result is integrated over time grid
            expects alpha0,beta0,phi0 to depend on tgird
        output units are cm^2-sr-MeV or cm^2-sr-MeV-s if tgrid supplied
        hAalph ais a scalar or 1-d numpy array of same length as alphagrid
        see module glossary for input meanings
        """
        h = self.halphabeta(alpha0,beta0,phib,Egrid,alphagrid,betagrid,tgrid,**kwargs)
        return h.sum(axis=1) # integrate over beta

class CR_Esep(ChannelResponse):
    """
    CR_Esep
    abstract base class for all energy-separable responses
    """
    
    @classmethod
    def is_mine(cls,*args,**kwargs):
        return False # abstract class is never the right answer
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.er = EnergyResponse(**kwargs)
        self.ar = AngleResponse(**kwargs)
    def R(self,E,theta,phi):
        """*INHERIT*"""
        return self.er.RE(E)*self.ar.A(theta,phi)
        raise NotImplementedError # abstract base class
    def _merge_hEhA(hE,hA):
        """
        hEA = merge hE size (NE,) with hA size (N1,) or (N1,N2)
        into hEhA of size (NE,N1) or size (NE,N1,N2)
        hEA = hE*hA
        """
        hE = np.reshape(hE,(hE.size,np.ones(hA.ndim,type=int)))
        hA = np.reshape(hE,(1,*hA.shape))
        return hE*hA
    def hEthetaphi(self,Egrid,thetagrid,phigrid,**kwargs):
        """*INHERIT*"""
        if Egrid is None:
            if hasattr(self,'E_GRID'):
                Egrid = self.E_GRID
            else:
                raise ValueError('Egrid input cannot be None for object without its own E_GRID')
        hE = self.er.hE(Egrid,**kwargs)
        hA = self.ar.hAthetaphi(thetagrid,phigrid,**kwargs)
        return self._merge_hEhA(hE,hA)
    def hEtheta(self,Egrid,thetagrid,phigrid=None,**kwargs):
        """*INHERIT*"""
        if phigrid is None:
            if hasattr(self,'PH_GRID'):
                phigrid = self.PH_GRID
            else:
                phigrid = default_phigrid(**kwargs)
        hE = self.er.hE(Egrid,**kwargs)
        hA = self.ar.hAtheta(thetagrid,phigrid,**kwargs)
        return self._merge_hEhA(hE,hA)
    def hE(self,Egrid,thetagrid=None,phigrid=None,**kwargs):
        """*INHERIT*"""
        if thetagrid is None:
            if hasattr(self,'TH_GRID'):
                thetagrid = self.TH_GRID
            else:
                thetagrid = default_thetagrid(**kwargs)
        hE = self.er.hE(Egrid,**kwargs)
        return hE*self.hA0
    def hEalphabeta(self,alpha0,beta0,phib,Egrid,alphagrid,betagrid,tgrid=None,**kwargs):
        """*INHERIT*"""
        hE = self.er.hE(Egrid,**kwargs)
        hA = self.ar.hAalphabeta(alpha0,beta0,phib,alphagrid,betagrid,tgrid,**kwargs)
        return self._merge_hEhA(hE,hA)
    def hEalpha(self,alpha0,beta0,phib,Egrid,alphagrid,betagrid=None,tgrid=None,**kwargs):
        """*INHERIT*"""
        hE = self.er.hE(Egrid,**kwargs)
        hA = self.ar.hAalpha(alpha0,beta0,phib,alphagrid,betagrid,tgrid,**kwargs)
        return self._merge_hEhA(hE,hA)
    def hEiso(self,alpha0,beta0,phib,Egrid,alphagrid,betagrid=None,tgrid=None,**kwargs):
        """*INHERIT*"""
        hE = self.er.hE(Egrid,**kwargs)
        hA = self.ar.hA0*sum(make_deltas(tgrid,**kwargs))
        return hE*hA
    def hthetaphi(self,Egrid,thetagrid,phigrid,**kwargs):
        """*INHERIT*"""
        hE = self.er.hE0(Egrid,**kwargs)
        hA = self.ar.hAthetaphi(thetagrid,phigrid,**kwargs)
        return hE*hA
    def htheta(self,Egrid,thetagrid,phigrid=None,**kwargs):
        """*INHERIT*"""
        hE = self.er.hE0(Egrid,**kwargs)
        hA = self.ar.hAtheta(thetagrid,phigrid,**kwargs)
        return hE*hA
    def halphabeta(self,alpha0,beta0,phib,Egrid,alphagrid,betagrid,tgrid=None,**kwargs):
        """*INHERIT*"""
        hE = self.er.hE0(Egrid,**kwargs)
        hA = self.ar.hAalphabeta(alpha0,beta0,phib,alphagrid,betagrid,tgrid,**kwargs)
        return hE*hA
    def halpha(self,alpha0,beta0,phib,Egrid,alphagrid,betagrid=None,tgrid=None,**kwargs):
        """*INHERIT*"""
        hE = self.er.hE0(Egrid,**kwargs)
        hA = self.ar.hAalpha(alpha0,beta0,phib,alphagrid,betagrid,tgrid,**kwargs)
        return hE*hA

class CR_Omni(CR_Esep):
    """
    CR_Omni
    omnidirectional channel response
    """
    
    @classmethod
    def is_mine(cls,*args,**kwargs):
        return ('RESP_TYPE' in kwargs) and (kwargs['RESP_TYPE'] == '[E]')
    # no further customization required
    
class CR_Esep_sym(CR_Esep):
    """
    CR_Esep_sym
    energy-separable, phi-symmetric channel response
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        return ('RESP_TYPE' in kwargs) and (kwargs['RESP_TYPE'] == '[E][TH]')
    # no further customization required
    
class CR_Esep_asym(CR_Esep):
    """
    CR_Esep_asym
    energy-separable, phi-asymmetric channel response
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        return ('RESP_TYPE' in kwargs) and (kwargs['RESP_TYPE'] == '[E][TH,PH]')
    # no further customization required

class CR_Table_sym(ChannelResponse):
    """
    CR_Table_sym
    energy-inseparable, phi-symmetric channel response
    response defined by 2-D table
    initializer requires ET_TYPE, E_GRID, TH_GRID, and R
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        keyword_check('RESP_TYPE',**kwargs)
        return ('RESP_TYPE' in kwargs) and \
                (kwargs['RESP_TYPE'] == '[E,TH]') \
                ('ET_TYPE' in kwargs) and \
                (kwargs['ET_TYPE'] == 'TBL')
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        keyword_check_numeric('E_GRID','TH_GRID','R',**kwargs)
        for arg in ['E_GRID','TH_GRID']:
            setattr(self,arg,squeeze(kwargs[arg]))
            if not validate_grid(getattr(self,arg)):
                raise ValueError('%s is not a valid grid: 1-d, unique' % arg)
        # TODO convert E_GRID to MeV
        self._R = squeeze(kwargs['R']) # TODO: convert to cm^2
        if self._R.shape != (self.E_GRID.size,self.TH_GRID.size):
            raise ArgSizeError('R is not shape (E_GRID x TH_GRID')
        self._Rinterpolator = None
    def R(self,E,theta,phi):
        """*INHERIT*"""
        if self._Rinterpolator is None:
            # extrapolate beyond grid with nearest edge with zero
            self._Rinterpolator = RegularGridInterpolator((self.E_GRID,self.TH_GRID),self._R,'linear',bounds_error=False,fill_value=0.0)
        return self._Rinterpolator(E,theta)
    
class CR_Table_asym(ChannelResponse):
    """
    CR_Table_asym
    energy-inseparable, phi-asymmetric channel response
    response defined by 3-D table
    initializer requires ETP_TYPE, E_GRID, TH_GRID, PH_GRID, and R
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        return ('RESP_TYPE' in kwargs) and \
                (kwargs['RESP_TYPE'] == '[E,TH,PH]') \
                ('ETP_TYPE' in kwargs) and \
                (kwargs['ETP_TYPE'] == 'TBL')
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        keyword_check_numeric('E_GRID','TH_GRID','PH_GRID','R',**kwargs)
        for arg in ['E_GRID','TH_GRID','PH_GRID']:
            setattr(self,arg,squeeze(kwargs[arg]))
            if not validate_grid(getattr(self,arg)):
                raise ValueError('%s is not a valid grid: 1-d, unique' % arg)
        # TODO convert E_GRID to MeV
        self._R = squeeze(kwargs['R']) # TODO: convert to cm^2
        if self._R.shape != (self.E_GRID.size,self.TH_GRID.size,self.PH_GRID.size):
            raise ArgSizeError('R is not shape (E_GRID x TH_GRID x PH_GRID')
        self._Rinterpolator = None
        if (self.PH_GRID[0] !=0) or (self.PH_GRID[-1] != 360):
            raise ValueError('PH_GRID should span [0,360]')
    def R(self,E,theta,phi):
        """*INHERIT*"""
        if self._Rinterpolator is None:
            # extrapolate beyond grid with nearest edge with zero
            self._Rinterpolator = RegularGridInterpolator((self.E_GRID,self.TH_GRID,self.PH_GRID),self._R,'linear',bounds_error=False,fill_value=0.0)
        return self._Rinterpolator(E,theta,phi)

# inherit dosctrings for everythin descended from FactoryConstructorMixin
# includes ChannelResponse and its internal classes EnergyResponse and AngleResponse
inherit_docstrings(FactoryConstructorMixin)

if __name__ == '__main__':

    # factory constructor demo
    class A(FactoryConstructorMixin):  # abstract base class
        @classmethod
        def _is_mine(cls,*args,**kwargs):
            # abstract base class never right answer
            raise Exception("Provided initialization data did not correspond to any implemented subclass of '%s'" % cls.__name__)
        def __init__(self,*args,**kwargs):
            if 'type' in kwargs:
                self.type = kwargs['type']
            else:
                self.type = None
        def __str__(self):
            return 'Object of class %s type=%s' % (self.__class__.__name__,self.type)
            
    class B(A):
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return ('type' in kwargs) and (kwargs['type'] == 'B')
    class C(A):
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return ('type' in kwargs) and (kwargs['type'] == 'C')
    class D(C): # extends C, C and D co-exist
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return ('type' in kwargs) and (kwargs['type'] == 'D')
    
    class E(C):
        pass # replaces C in hierarchy
    
    try:    
        print('A()',A())
    except Exception as e:
        print(e)
    try:    
        print('A',A(type='A'))
    except Exception as e:
        print(e)
    print('B',A(type='B'))
    print('C',A(type='C'))
    print('D',A(type='D'))
    try:    
        print('E',A(type='E'))
    except Exception as e:
        print(e)
    
    # very simple test of RFL class instantiation
    print('Omni int',ChannelResponse(RESP_TYPE='[E]',E_TYPE='INT',E0=2))
    print('Tele_sym wide',ChannelResponse(RESP_TYPE='[E][TH]',E_TYPE='WIDE',E0=2,E1=3,TH_TYPE='CYL_TELE',R1=1,R2=2,D=3.0))
    print('Tele_asym diff',ChannelResponse(RESP_TYPE='[E][TH,PH]',E_TYPE='DIFF',E0=2,DE=1,TP_TYPE='RECT_TELE',W1=1.0,H1=2.0,W2=0.5,H2=1.0,D=3.0))
