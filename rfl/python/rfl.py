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
    alphagrid - alpha grid, degrees
    betagrid - beta grid, degrees
    thetagrid - theta grid, degrees
    phigrid - phi grid, degrees
    tgrid - time grid, seconds
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

"""
TODO:
    apply CROSSCALIB for hE and h* (but not hA*)
"""

# NOTE: in this code where the docstring is *INHERIT* that means
# the string *INHERIT* will be replaced by the docstring for the
# parent class's corresponding method

import numpy as np
from scipy.interpolate import interp1d

# utilities
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
        returns i = (2,), where lst[i[0]] <= q <= lst[i[1]]
    else:
        returns i = (len(q),2), where lst[i[:,0]] <= q <= lst[i[:,1]]
    i == -1 implies q out of range of lst
    """
    lst = np.array(lst).ravel()
    N = len(lst)
    is_scalar = np.isscalar(q)
    q = np.atleast_1d(q).ravel()
    i = np.full((len(q),2),np.nan)
    iq = np.isfinite(q)
    i[iq,0] = interp1d(lst,np.arange(N),'nearest',bounds_error=False)(q[iq])
    f = np.isfinite(i[:,0])
    if any(f):
        i[f,0] = i[f,0] - (q[f]<lst[i[f,0]]) # force to greatest list <= q
        i[f,1] = i[f,0]+1 # next neighbor to right
    # limit checks
    i[i>=N] = -1
    i[i<0] = -1
    i[np.logical_not(np.isfinite(i))] = -1
    if is_scalar:
        i = i[0,:]
    return i

def make_deltas(grid,int_method='trapz',**kwargs):
    """
    d = make_deltas(grid,int_method='trapz',...)
    d are the default weights for a given 1-d grid
    d is (N,), unless N<=1, in which case d is scalar
    for empty grid, returns 1
    for singleton grid, returns grid (assumes it's actually dt)
    keyword args besides int_method ignored
    """
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

def squeeze(val):
    """
    val = squeeze(val)
    squeeze a value by removing trailing singleton dimensions
    and if a single value, make it a scalar
    """
    val = np.squeeze(val)
    if val.size==1:
        val = np.asscalar(val) # reduce singleton to scalar
    return val
    

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
 

class ChannelResponse(FactoryConstructorMixin):
    """
    Virtual base class for channel responses
    Constructor is a factory that will return 
    appropriate subclass from response dictionary
    """
        
    @classmethod
    def is_mine(cls,*args,**kwargs):
        """*INHERIT*"""
        keyword_check('RESP_TYPE',**kwargs)
        raise RFLException('Supplied keywords to not define a recognized ChannelResponse')

    def __init__(self,**kwargs):
        print('Class not implemented yet: ' + self.__class__.__name__)
    
    def copy(self):
        """
        channel.copy() return a copy of this object
        """
        return ChannelResponse(**self.to_dict())

    def to_dict(self):
        """
        channel.to_dict() return a nested dict of this object
        """
        raise NotImplementedError

    class EnergyResponse(FactoryConstructorMixin):
        """
        EnergyResponse
        class representing energy responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            """*INHERIT*"""
            keyword_check('E_TYPE',**kwargs)
            raise RFLException('Supplied keywords to not define a recognized EnergyResponse')
        
        def __init__(self,**kwargs):
            if 'CROSSCALIB' not in kwargs:
                print('CROSSCALIB not found in channel description, assuming unity')
                self.CROSSCALIB = 1.0 # default
            else:
                self.CROSSCALIB = kwargs['CROSSCALIB']
            
        def RE(self,E):
            """
            R = RE(E)
            Channel response at energy E (w/o CROSSCALIB)
            output is dimensionless
            """
            raise NotImplementedError
        def hE(self,Egrid,**kwargs):
            """
            hE = .hE(Egrid,...)
            return energy response weights on Egrid
            CROSSCALIB applied
            """
            raise NotImplementedError

    class ER_Diff(EnergyResponse):
        """
        ER_Diff
        class representing differential energy channel response
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            """*INHERIT*"""
            keyword_check('E_TYPE',**kwargs)
            return (kwargs['E_TYPE'] == 'DIFF')

        def __init__(self,**kwargs):
            super().__init__(**kwargs)
            keyword_check_numeric('EPS','E0','DE',**kwargs)            
            self.EPS = squeeze(kwargs['EPS'])
            self.E0 = squeeze(kwargs['E0']) # TODO - convert to MeV
            self.DE = squeeze(kwargs['DE']) # TODO - convert to MeV
            self.hE0 = self.DE*self.EPS/self.CROSSCALIB # hE for flat spectrum

        def RE(self,E):
            """*INHERIT*"""
            return (E==self.E0)*self.EPS # DE gets ignored, which is bad, but unaviodable

        def hE(self,Egrid,**kwargs):
            """*INHERIT*"""
            Egrid = np.array(Egrid)
            NE = len(Egrid)
            hE = np.zeros(Egrid.shape) 
            I = get_list_neighbors(Egrid,self.E0)
            # Egrid(I(1)) <= inst_info.E0 < Egrid(I(2))
            # or I(2) = -1 or both are -1
            
            # left side
            i = I[0]
            if (i>=0) and (i<NE-1): # E[i] <= E0 < E[i+1]
                hE[i] = self.dE*self.EPS*(Egrid[i+1]-self.E0)/(Egrid[i+1]-Egrid[i])
            
            # right side
            i = I[1]
            if (i>=1) and (i<NE): # E[i-1] < E0 < E[i]
                hE[i] = self.dE*self.EPS*(self.E0-Egrid[i-1])/(Egrid[i]-Egrid[i-1])
                
            return hE / self.CROSSCALIB

            
    class ER_Int(EnergyResponse):
        """
        ER_Int
        class representing integral energy channel response
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            """*INHERIT*"""
            keyword_check('E_TYPE',**kwargs)
            return (kwargs['E_TYPE'] == 'INT')

        def __init__(self,**kwargs):
            super().__init__(**kwargs)
            keyword_check_numeric('EPS','E0',**kwargs)            
            self.EPS = squeeze(kwargs['EPS'])
            self.E0 = squeeze(kwargs['E0']) # TODO - convert to MeV
            self.hE0 = self.DE*self.EPS/self.CROSSCALIB  # should be inf for flat spectrum, but assume energy bandwidth==1

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
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            """*INHERIT*"""
            keyword_check('E_TYPE',**kwargs)
            return (kwargs['E_TYPE'] == 'WIDE')

        def __init__(self,**kwargs):
            super().__init__(**kwargs)
            keyword_check_numeric('EPS','E0','E1',**kwargs)            
            self.EPS = squeeze(kwargs['EPS'])
            self.E0 = squeeze(kwargs['E0']) # TODO - convert to MeV
            self.E1 = squeeze(kwargs['E1']) # TODO - convert to MeV
            self.hE0 = (self.E1-self.E0)*self.EPS/self.CROSSCALIB  # for flat spectrum

        def RE(self,E):
            """*INHERIT*"""
            return ((E>=self.E0) & (E<=self.E1))*self.EPS

        def hE(self,Egrid,**kwargs):
            """*INHERIT*"""
            
            # treat as difference between two integral channels
            low = ChannelResponse.DE_Int(**kwargs)
            high = {**kwargs,'E0':self.E1}
            return low.hE(Egrid)-high.hE(Egrid)
            
    class ER_Table(EnergyResponse):
        """
        ER_Table
        class representing tabular energy channel response
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            """*INHERIT*"""
            keyword_check('E_TYPE',**kwargs)
            return (kwargs['E_TYPE'] == 'TBL')

        def __init__(self,**kwargs):
            super().__init__(**kwargs)
            keyword_check_numeric('EPS','E_GRID',**kwargs)
            self.EPS = squeeze(kwargs['EPS'])
            self.E_GRID = squeeze(kwargs['E_GRID'])
            self.hE0 = np.sum(self.E_GRID*self.EPS*make_deltas(self.E_GRID))/self.CROSSCALIB  # for flat spectrum
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
        when initialized with no arguments, returns a null response
        
        properties:
            bidirectional - bool, is the angular response bidirectional?
            backward - backward angular response (null class when not bidirectional)
            area - detector area, cm^2
            G - geometric factor for forward hemisphere particles (theta<=90) cm^2 sr
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
            if kwargs: # backward present only when called with initialization data
                self.backward = ChannelResponse.AngleResponse() # backward response is null by default
            self.bidirectional = ('BIDIRECTIONAL' in kwargs) and (kwargs['BIDIRECTIONAL'] == 'TRUE')
            self.area = 0.0  # used by null
            self.G = 0.0 # used by null
        def hAalphabeta(self,alpha0,beta0,phib,alphagrid,betagrid,tgrid=None,**kwargs):
            """
            hAalphabeta = .hAalphabeta(alpha0,beta0,phib,alphagrid,betagrid,tgrid=None,...)
            return angular response weights on alpha x beta grid
            see module glossary for input meanings
            if tgrid is supplied, result is integrated over time grid
                expects alpha0,beta0,phi0 to depend on tgird
            output units are cm^2 or cm^2-s if tgrid supplied
            """
            return 0.0 # used by null
        def hAalpha(self,alpha0,beta0,phib,alphagrid,betagrid=None,tgrid=None,**kwargs):
            """
            hAalpha = .hAalpha(alpha0,beta0,phib,alphagrid,betagrid=None,tgrid=None,...)
            return angular response weights on alpha grid (integrated over beta)
            see module glossary for input meanings
            if betagrid is None, supplies numeric beta grid when needed
            if tgrid is supplied, result is integrated over time grid
                expects alpha0,beta0,phi0 to depend on tgird
            output units are cm^2 or cm^2-s if tgrid supplied
            """
            return 0.0 # used by null
        def A(self,theta,phi):
            """
            A = .A(theta,phi)
            return angular response (effective area) at specified theta,phi cm^2
            should add forward and backward response together
            see module glossary for input meanings
            """
            return 0.0 # used by null
        def hAthetaphi(self,thetagrid,phigrid,**kwargs):
            """
            hAthetaphi = .hAthetaphi(thetagrid,phigrid,...)
            return angular response weights on theta x phi grid, cm^2
            see module glossary for input meanings
            """
            return 0.0 # used by null
        def hAtheta(self,thetagrid,phigrid=None,**kwargs):
            """
            hAtheta = .hAtheta(thetagrid,phigrid,...)
            return angular response weights on theta grid, cm^2
            integrated over phigrid (or 0-2pi if None)
            see module glossary for input meanings
            """
            return 0.0 # used by null
        @property
        def hA0(self):
            """
            return the nominal geometric factor in cm^2 sr
            including bidirectional effects
            """
            return self.G + self.backward.G

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
        def hAthetaphi(self,thetagrid,phigrid,**kwargs):
            """*INHERIT*"""
            dcostheta = make_deltas(-np.cos(np.radians(thetagrid)),**kwargs)
            dphi = make_deltas(np.radians(phigrid),**kwargs)
            return self.A(thetagrid,phigrid)*dcostheta*dphi
        def hAtheta(self,thetagrid,phigrid=None,**kwargs):
            """*INHERIT*"""
            dcostheta = make_deltas(-np.cos(np.radians(thetagrid)),**kwargs)
            if phigrid is None:
                dphi = 2*np.phi
                phigrid = 0.0 # dummy value
            else:
                dphi = make_deltas(np.radians(phigrid),**kwargs).sum()
            return self.A(thetagrid,phigrid)*dcostheta*dphi

    class AR_Omni(AR_csym):
        """
        AR_Omni omnidirectional response
        class representing omnidirectional angular responses
        Note: an omni is a sphere (same apparent area from all directions). 
        For flat, single-element detectors, use DISK and SLAB
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            """*INHERIT*"""
            return keyword_check_bool('RESP_TYPE',**kwargs) and \
                (kwargs['RESP_TYPE'] == '[E]')
        def __init__(self,**kwargs):
            super().__init__(**kwargs)
            if ('G' in kwargs) and kwargs['G'] is not None:
                # derived class may provide G as none
                keyword_check_numeric('G',**kwargs)
                self.G = squeeze(kwargs['G'])            
            # assum a half omni
            self.area = self.G/2/np.pi # equivalent sphere's cross-section 
            if self.bidirectional:
                # backwards is also omni
                backw = kwargs.copy()
                backw['BIDIRECTIONAL'] = 'FALSE'
                self.backward = self.__class__(**backw)
        def A(self,theta,phi):
            """*INHERIT*"""
            return self.area*(theta<=90) + self.backward.area*(theta>90)

    class AR_SingleElement(AR_Omni):
        """
        AR_SingleElement
        class representing single-element detectors (i.e., w/ cos projection effect)
        derived class must define 
          self.area, the detector area
          self.G, the geometric factor
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            "*INHERIT*"
            return False # abstract class, never the right answer
        def __init__(self,**kwargs):
            super().__init__(G=None,**kwargs) # supply fake G, area and G will be set in derived class
        def A(self,theta,phi):
            """*INHERIT*"""
            return super().A(theta,phi)*np.abs(np.cos(np.radians(theta,phi)))
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
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return keyword_check_bool('TH_TYPE',**kwargs) and \
                (kwargs['TH_TYPE'] == 'DISK')
        def __init__(self,**kwargs):
            super().__init__(**kwargs)
            keyword_check_numeric('W1','H1',**kwargs)
            self.W1 = squeeze(kwargs['W1'])
            self.H1 = squeeze(kwargs['H1'])
            self.area = self.W1*self.H1

    class AR_Slab(AR_SingleElement):
        """
        AR_Slab
        class representing rectangular slab angular responses
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return keyword_check_bool('TH_TYPE',**kwargs) and \
                (kwargs['TH_TYPE'] == 'SLAB')
        def __init__(self,**kwargs):
            super().__init__(**kwargs)
            keyword_check_numeric('R1',**kwargs)
            self.R1 = squeeze(kwargs['R1'])
            self.area = np.pi*self.R1**2

    class AR_Tele_sym(AR_csym):
        """
        AR_Tele_sym
        class representing phi-symmetric telescope angular responses
        Using Sullivan's 1971 paper as updated in ~2010
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return keyword_check_bool('TH_TYPE',**kwargs) and \
                (kwargs['TH_TYPE'] == 'CYL_TELE')
        def __init__(self,**kwargs):
            super().__init__(**kwargs)
            keyword_check_numeric('R1','R2','D',**kwargs)
            self.R1 = squeeze(kwargs['R1'])
            self.R2 = squeeze(kwargs['R1'])
            self.D = squeeze(kwargs['D'])
            self.Rs = min(self.R1,self.R2)
            self.thetac = np.degrees(np.atan(np.abs(self.R1-self.R2)/self.D))
            self.thetam = np.degrees(np.atand((self.R1+self.R2)/self.D))
            self.area = np.pi*self.Rs**2 # area at normal incidence
            tmp = self.R1**2+self.R2**2+self.D**2
            self.G = np.pi**2/2*(tmp-np.sqrt(tmp**2-4*self.R1**2*self.R2**2)) # eqn 8
            if self.bidirectional:
                backw = kwargs.copy()
                backw['BIDIRECTIONAL'] = 'FALSE'
                backw['R1'] = self.R1
                backw['R2'] = self.R2
                self.backward = ChannelResponse.AR_Tele_sym(**backw)
        def A(self,theta,phi):
            """*INHERIT*"""
            # uses Sullivan's updated paper, equation 10
            A = np.zeros(np.broadcast(theta,phi).shape)
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
            return A

    class AR_Table_sym(AngleResponse):
        """
        AR_Table_sym
        class representing phi-symmetric tabular angular responses
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return keyword_check_bool('TH_TYPE',**kwargs) and \
                (kwargs['TH_TYPE'] == 'TBL')

    class AR_Pinhole(AngleResponse):
        """
        AR_Pinhole
        class representing pinhole angular responses
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):            
            return keyword_check_bool('TH_TYPE',**kwargs) and \
                (kwargs['TH_TYPE'] == 'PINHOLE')

    class AR_Tele_asym(AngleResponse):
        """
        AR_Tele_asym
        class representing phi-asymmetric telescope angular responses
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return keyword_check_bool('TP_TYPE',**kwargs) and \
                (kwargs['TP_TYPE'] == 'RECT_TELE')

    class AR_Table_asym(AngleResponse):
        """
        AR_Table_asym
        class representing phi-asymmetric tabular angular responses
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return keyword_check_bool('TP_TYPE',**kwargs) and \
                (kwargs['TP_TYPE'] == 'TBL')


class CR_Omni(ChannelResponse):
    """
    CR_Omni
    omnidirectional channel response
    """
    
    @classmethod
    def is_mine(cls,*args,**kwargs):
        keyword_check('RESP_TYPE',**kwargs)
        return (kwargs['RESP_TYPE'] == '[E]')
    
class CR_Esep_sym(CR_Omni):
    """
    CR_Esep_sym
    energy-separable, phi-symmetric channel response
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        keyword_check('RESP_TYPE',**kwargs)
        return (kwargs['RESP_TYPE'] == '[E][TH]')
    
class CR_Esep_asym(CR_Esep_sym):
    """
    CR_Esep_asym
    energy-separable, phi-asymmetric channel response
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        keyword_check('RESP_TYPE',**kwargs)
        return (kwargs['RESP_TYPE'] == '[E][TH,PH]')

class CR_insep_sym(CR_Omni):
    """
    CR_insep_sym
    energy-inseparable, phi-symmetric channel response
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        keyword_check('RESP_TYPE',**kwargs)
        return (kwargs['RESP_TYPE'] == '[E,TH]')
    
class CR_insep_asym(CR_insep_sym):
    """
    CR_insep_asym
    energy-inseparable, phi-asymmetric channel response
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        keyword_check('RESP_TYPE',**kwargs)
        return (kwargs['RESP_TYPE'] == '[E,TH,PH]')

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
    print('Omni',ChannelResponse(RESP_TYPE='[E]'))
    print('CR_insep_sym',ChannelResponse(RESP_TYPE='[E,TH]'))
    print('CR_Esep_asym',ChannelResponse(RESP_TYPE='[E][TH,PH]'))
    print('AR_Tele_asym',ChannelResponse.AngleResponse(TP_TYPE='RECT_TELE'))
