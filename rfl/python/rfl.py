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
    use __subclasses__ in create factory so user-defined sub-classes work
    do this in the __new__ constructor
    use a class method is_mine to check if this class is appropriate
    when sbuclassing, if is_mine is not overloaded, subclass will always
    take precedence over its parent class
"""

import numpy as np
from scipy.interpolate import interp1d



# utilities
def get_list_neighbors(lst,q):
    """
    i = rfl_get_list_neighbors(lst,q)
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
        raise RFL_ValueEx('Unknown int_method "%s"' % int_method)
    
    return d


# classes

class RFL_Exception(Exception):
    """Unknown exception (base class of other RFL Exceptions)"""
    def __init__(self,message = 'Unknown RFL Exception'):
        super().__init__(message)

class RFL_ArgSizeEx(RFL_Exception):
    """Incorrect/Inconsistent Argument Size"""
    def __init__(self,message = 'Incorrect/Inconsistent Argument Size'):
        super().__init__(message)

class RFL_KeywordEx(RFL_Exception):
    """Missing/Invalid Keyword"""
    def __init__(self,message = 'Missing/Invalid Keyword'):
        super().__init__(message)

class RFL_ValueEx(RFL_Exception):
    """Unexpected/Invalid Value"""
    def __init__(self,message = 'Unexpected/Invalid Value'):
        super().__init__(message)

def keyword_check(*keywords,**kwargs):
    """
    bool = keyword_check(*keywords,**kwargs)
    check one or more keywords in kwargs
    if present and str, return True
    otherwise raise appropriate exception
    """
    for key in keywords:
        if key not in kwargs: raise RFL_KeywordEx('Keyword %s required, absent' % key)
        if not isinstance(kwargs[key],str): raise RFL_KeywordEx('%s must be a string(str)' % key)
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

class FactoryConstructorMixin(object):
    """
    class FactoryConstructorMixin
    a mixin to support letting an object constructor
    traverse its subclass tree to instantiate the correct
    subclass based on the initialization arguments

    Class trees that inherit from FactoryConstructorMixin
    should implement is_mine in every subclass to "claim" the
    initialiation data and insantiate that subclass rather than
    the parent constructor thta was called.    

    @classmethod
    def is_mine(cls,*args,**kwargs):
        # search kwargs dict/tree to determine if this subclass should handle
        # the data, and if so, return True
    
    Any class that inherites directly from FactoryConstructorMixin
    should implement is_mine that raises an exception rather than returning False
    
    An subclass that does not define is_mine will replace its parent in the
    factory hierarchy (use this approach to modify the behavior of a subclass).
    
    Any subclass that defines its own is_mine will extend its parent into a
    distinct, new subclass in the factory hierarchy (use this approach if
    parent and chiled are alternatives that should co-exist)
    
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
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
        keyword_check('RESP_TYPE',**kwargs)
        raise RFL_Exception('Supplied keywords to not define a recognized ChannelResponse')

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
            keyword_check('E_TYPE',**kwargs)
            raise RFL_Exception('Supplied keywords to not define a recognized EnergyResponse')
        
        def __init__(self,**kwargs):
            if 'CROSSCALIB' not in kwargs:
                print('CROSSCALIB not found in channel description, assuming unity')
                self.CROSSCALIB = 1.0 # default
            else:
                self.CROSSCALIB = kwargs['CROSSCALIB']
            
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
            keyword_check('E_TYPE',**kwargs)
            return (kwargs['E_TYPE'].upper() == 'DIFF')

        def __init__(self,**kwargs):
            super().__init__(**kwargs)
            self.EPS = kwargs['EPS']
            self.E0 = kwargs['E0'] # TODO - convert to MeV
            self.DE = kwargs['DE'] # TODO - convert to MeV
            self.hE0 = self.DE*self.EPS/self.CROSSCALIB # hE for flat spectrum

        def RE(self,E):
            """
            R = RE(E)
            Channel response at energy E (w/o CROSSCALIB)
            output is dimensionless
            """
            return (E==self.E0)*self.EPS # DE gets ignored, which is bad, but unaviodable

        def hE(self,Egrid,**kwargs):
            """
            hE = .hE(Egrid,...)
            return energy response weights on Egrid
            CROSSCALIB applied
            """
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
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            keyword_check('E_TYPE',**kwargs)
            return (kwargs['E_TYPE'].upper() == 'INT')

        def __init__(self,**kwargs):
            super().__init__(**kwargs)
            self.EPS = kwargs['EPS']
            self.E0 = kwargs['E0'] # TODO - convert to MeV
            self.hE0 = self.DE*self.EPS/self.CROSSCALIB  # should be inf for flat spectrum, but assume energy bandwidth==1

        def RE(self,E):
            """
            R = RE(E)
            Channel response at energy E (w/o CROSSCALIB)
            output is dimensionless
            """
            return (E>=self.E0)*self.EPS

        def hE(self,Egrid,**kwargs):
            """
            hE = .hE(Egrid,...)
            return energy response weights on Egrid
            CROSSCALIB applied
            """
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
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            keyword_check('E_TYPE',**kwargs)
            return (kwargs['E_TYPE'].upper() == 'WIDE')

        def __init__(self,**kwargs):
            super().__init__(**kwargs)
            self.EPS = kwargs['EPS']
            self.E0 = kwargs['E0'] # TODO - convert to MeV
            self.E1 = kwargs['E1'] # TODO - convert to MeV
            self.hE0 = (self.E1-self.E0)*self.EPS/self.CROSSCALIB  # for flat spectrum

        def RE(self,E):
            """
            R = RE(E)
            Channel response at energy E (w/o CROSSCALIB)
            output is dimensionless
            """
            return ((E>=self.E0) & (E<=self.E1))*self.EPS

        def hE(self,Egrid,**kwargs):
            """
            hE = .hE(Egrid,...)
            return energy response weights on Egrid
            CROSSCALIB applied
            """
            
            # treat as difference between two integral channels
            low = ChannelResponse.DE_Int(**kwargs)
            high = {**kwargs,'E0':self.E1}
            return low.hE(Egrid)-high.hE(Egrid)
            
    class ER_Table(EnergyResponse):
        """
        ER_Table
        class representing tabular energy channel response
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            keyword_check('E_TYPE',**kwargs)
            return (kwargs['E_TYPE'].upper() == 'TBL')

        def __init__(self,**kwargs):
            super().__init__(**kwargs)
            self.EPS = kwargs['EPS']
            self.E_GRID = kwargs['E_GRID']
            self.hE0 = np.sum(self.E_GRID*self.EPS*make_deltas(self.E_GRID))/self.CROSSCALIB  # for flat spectrum
            self._RE = None

        def RE(self,E):
            """
            R = RE(E)
            Channel response at energy E (w/o CROSSCALIB)
            output is dimensionless
            """
            if self._RE is None:
                # extrapolate beyond grid with nearest edge (first/last) value
                self._RE = interp1d(self.E_GRID,self.EPS,'linear',bounds_error=False,fill_value=[self.EPS[0],self.EPS[-1]],assum_sorted=True)
            return self._RE(E)

        def hE(self,Egrid,**kwargs):
            """
            hE = .hE(Egrid,...)
            return energy response weights on Egrid
            CROSSCALIB applied
            """
            
            dE = make_deltas(Egrid,**kwargs)
            hE = self.RE(Egrid)*dE/self.CROSSCALIB
            return hE

    class AngleResponse(FactoryConstructorMixin):
        """
        AngleResponse
        class representing angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        # omni, TBL,CYL_TELE,RECT_TELE,SLAB,DISK,PIN_HOLE
        @classmethod
        def is_mine(cls,*args,**kwargs):
            raise RFL_Exception('Supplied keywords to not define a recognized AngleResponse')
        def hAalphabeta(self,alpha0,beta0,phib,alphagrid,betagrid,tgrid=None,**kwargs):
            """
            hAalphabeta = .hAalphabeta(alpha0,beta0,phib,alphagrid,betagrid,tgrid=None,...)
            return angular response weights on alpha x beta grid
            see module glossar for input meanings
            if tgrid is supplied, result is integrated over time grid
                expects alpha0,beta0,phi0 to depend on tgird
            output units are cm^2 or cm^2-s if tgrid supplied
            """
            raise NotImplementedError
        def hAalpha(self,alpha0,beta0,phib,alphagrid,betagrid=None,tgrid=None,**kwargs):
            """
            hAalpha = .hAalpha(alpha0,beta0,phib,alphagrid,betagrid=None,tgrid=None,...)
            return angular response weights on alpha grid (integrated over beta)
            see module glossar for input meanings
            if betagrid is None, supplies numeric beta grid when needed
            if tgrid is supplied, result is integrated over time grid
                expects alpha0,beta0,phi0 to depend on tgird
            output units are cm^2 or cm^2-s if tgrid supplied
            """
            raise NotImplementedError
        def hthetaphi(self,thetagrid,phigrid,**kwargs):
            """
            hAalphabeta = .hAalphabeta(thetagrid,phigrid,...)
            return angular response weights on theta x phi grid, cm^2
            see module glossar for input meanings
            """
            raise NotImplementedError
        def hAtheta(self,thetagrid,phigrid=None,**kwargs):
            """
            halphabeta = .halphabeta(thetagrid,phigrid,...)
            return angular response weights on theta grid, cm^2
            see module glossar for input meanings
            """
            raise NotImplementedError
        @property
        def G0(self):
            """
            return the nominal geometric factor in cm^2 sr
            """
            raise NotImplementedError

    class AR_Omni(AngleResponse):
        """
        AR_Omni
        class representing omnidirectional angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return keyword_check_bool('RESP_TYPE',**kwargs) and \
                (kwargs['RESP_TYPE'].upper() == '[E]')

    class AR_Pinhole(AngleResponse):
        """
        AR_Pinhole
        class representing pinhole angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):            
            return keyword_check_bool('TH_TYPE',**kwargs) and \
                (kwargs['TH_TYPE'].upper() == 'PINHOLE')

    class AR_Disk(AngleResponse):
        """
        AR_Disk
        class representing disk angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return keyword_check_bool('TH_TYPE',**kwargs) and \
                (kwargs['TH_TYPE'].upper() == 'DISK')

    class AR_Slab(AngleResponse):
        """
        AR_Slab
        class representing rectangular slab angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return keyword_check_bool('TH_TYPE',**kwargs) and \
                (kwargs['TH_TYPE'].upper() == 'SLAB')

    class AR_Tele_sym(AngleResponse):
        """
        AR_Tele_sym
        class representing phi-symmetric telescope angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return keyword_check_bool('TH_TYPE',**kwargs) and \
                (kwargs['TH_TYPE'].upper() == 'CYL_TELE')

    class AR_Table_sym(AngleResponse):
        """
        AR_Table_sym
        class representing phi-symmetric tabular angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return keyword_check_bool('TH_TYPE',**kwargs) and \
                (kwargs['TH_TYPE'].upper() == 'TBL')

    class AR_Tele_asym(AngleResponse):
        """
        AR_Tele_asym
        class representing phi-asymmetric telescope angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return keyword_check_bool('TP_TYPE',**kwargs) and \
                (kwargs['TP_TYPE'].upper() == 'RECT_TELE')

    class AR_Table_asym(AngleResponse):
        """
        AR_Table_asym
        class representing phi-asymmetric tabular angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def is_mine(cls,*args,**kwargs):
            return keyword_check_bool('TP_TYPE',**kwargs) and \
                (kwargs['TP_TYPE'].upper() == 'TBL')


class CR_Omni(ChannelResponse):
    """
    CR_Omni
    omnidirectional channel response
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    
    @classmethod
    def is_mine(cls,*args,**kwargs):
        keyword_check('RESP_TYPE',**kwargs)
        return (kwargs['RESP_TYPE'].upper() == '[E]')
    
class CR_Esep_sym(CR_Omni):
    """
    CR_Esep_sym
    energy-separable, phi-symmetric channel response
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        keyword_check('RESP_TYPE',**kwargs)
        return (kwargs['RESP_TYPE'].upper() == '[E][TH]')
    
class CR_Esep_asym(CR_Esep_sym):
    """
    CR_Esep_asym
    energy-separable, phi-asymmetric channel response
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        keyword_check('RESP_TYPE',**kwargs)
        return (kwargs['RESP_TYPE'].upper() == '[E][TH,PH]')

class CR_insep_sym(CR_Omni):
    """
    CR_insep_sym
    energy-inseparable, phi-symmetric channel response
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        keyword_check('RESP_TYPE',**kwargs)
        return (kwargs['RESP_TYPE'].upper() == '[E,TH]')
    
class CR_insep_asym(CR_insep_sym):
    """
    CR_insep_asym
    energy-inseparable, phi-asymmetric channel response
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def is_mine(cls,*args,**kwargs):
        keyword_check('RESP_TYPE',**kwargs)
        return (kwargs['RESP_TYPE'].upper() == '[E,TH,PH]')

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
