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


# class factory demo

class A(object):
    @classmethod
    def is_mine(cls,**kwargs):
        return None # virtual base class never right answer
    def __new__(cls,**kwargs):
        # check children first
        for c in cls.__subclasses__():
            res = c.__new__(c,**kwargs)
            if res: return res
        if cls.is_mine(**kwargs):
            return cls
        return None            
        
class B(A):
    @classmethod
    def is_mine(cls,**kwargs):
        return ('type' in kwargs) and (kwargs['type'] == 'B')
class C(A):
    @classmethod
    def is_mine(cls,**kwargs):
        return ('type' in kwargs) and (kwargs['type'] == 'C')
class D(C):
    @classmethod
    def is_mine(cls,**kwargs):
        return ('type' in kwargs) and (kwargs['type'] == 'D')

print('A()',A())
print('A',A(type='A'))
print('B',A(type='B'))
print('C',A(type='C'))
print('D',A(type='D'))
raise Exception('stop')

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

class ChannelResponse(object):
    """
    Virtual base class for channel responses
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def create(cls,**kwargs):
        """
        channel = ChannelResponse.create(**kwargs)
        returns an initialized response object of the specified type
        """
        if 'RESP_TYPE' not in kwargs: raise RFL_KeywordEx('RESP_TYPE required, absent')
        if not isinstance(kwargs['RESP_TYPE'],str): raise RFL_KeywordEx('RESP_TYPE must be a string(str)')
        for c in cls.__subclasses__():
            obj = c.create(**kwargs)
            if obj: return obj
        raise RFL_Exception('Supplied keywords to not define a recognized ChannelResponse')

    def __init__(self,**kwargs):
        print('Not Implemented yet: ' + self.__class__.__name__)
    
    def copy(self):
        """
        channel.copy() return a copy of this object
        """
        ChannelResponse.create(**self.to_dict())

    def to_dict(self):
        """
        channel.to_dict() return a nested dict of this object
        """
        raise NotImplementedError

    class EnergyResponse(object):
        """
        EnergyResponse
        class representing energy responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def create(cls,**kwargs):
            if 'E_TYPE' not in kwargs: raise RFL_KeywordEx('E_TYPE required, absent')
            if not isinstance(kwargs['E_TYPE'],str): raise RFL_KeywordEx('E_TYPE must be a string(str)')
            # TODO: rewrite this to use create on subclasses
            E_TYPE = kwargs['E_TYPE'].upper()
            if E_TYPE == 'TBL':
                return ChannelResponse.ER_Table.create(**kwargs)
            elif E_TYPE == 'INT':
                return ChannelResponse.ER_Int.create(**kwargs)
            elif E_TYPE == 'WIDE':
                return ChannelResponse.ER_Wide.create(**kwargs)
            elif E_TYPE == 'DIFF':
                return ChannelResponse.ER_Diff.create(**kwargs)
            else:
                raise RFL_ValueEx('E_TYPE had unexpected value %s' % E_TYPE)
        
        def __init__(self,**kwargs):
            if 'CROSSCALIB' not in kwargs:
                raise Warning('CROSSCALIB not found in channel description, assuming unity')
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
        def create(cls,**kwargs):
            """
            Create an ER_Diff object
            """
            return cls(**kwargs)

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
        def create(cls,**kwargs):
            """
            Create an ER_Int object
            """
            return cls(**kwargs)

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
        def create(cls,**kwargs):
            """
            Create an ER_Wide object
            """
            return cls(**kwargs)

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
        def create(cls,**kwargs):
            """
            Create an ER_Table object
            """
            return cls(**kwargs)

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

    class AngleResponse(object):
        """
        AngleResponse
        class representing angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        # omni, TBL,CYL_TELE,RECT_TELE,SLAB,DISK,PIN_HOLE
        @classmethod
        def create(cls,**kwargs):
            # TODO: write this to use create on subclasses
            raise NotImplementedError
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
        def create(cls,**kwargs):
            raise NotImplementedError

    class AR_Pinhole(AngleResponse):
        """
        AR_Pinhole
        class representing pinhole angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def create(cls,**kwargs):
            raise NotImplementedError

    class AR_Disk(AngleResponse):
        """
        AR_Disk
        class representing disk angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def create(cls,**kwargs):
            raise NotImplementedError

    class AR_Slab(AngleResponse):
        """
        AR_Slab
        class representing rectangular slab angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def create(cls,**kwargs):
            raise NotImplementedError

    class AR_Tele_sym(AngleResponse):
        """
        AR_Tele_sym
        class representing phi-symmetric telescope angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def create(cls,**kwargs):
            raise NotImplementedError

    class AR_Table_sym(AngleResponse):
        """
        AR_Table_sym
        class representing phi-symmetric tabular angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def create(cls,**kwargs):
            raise NotImplementedError

    class AR_Tele_asym(AngleResponse):
        """
        AR_Tele_asym
        class representing phi-asymmetric telescope angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def create(cls,**kwargs):
            raise NotImplementedError

    class AR_Table_asym(AngleResponse):
        """
        AR_Table_asym
        class representing phi-asymmetric tabular angular responses
        use factory class method .create() to create
        appropriate subclass from response dictionary
        """
        @classmethod
        def create(cls,**kwargs):
            raise NotImplementedError


class CR_Omni(ChannelResponse):
    """
    CR_Omni
    omnidirectional channel response
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def create(cls,**kwargs):
        if 'RESP_TYPE' not in kwargs: raise RFL_KeywordEx('RESP_TYPE required, absent')
        if not isinstance(kwargs['RESP_TYPE'],str): raise RFL_KeywordEx('RESP_TYPE must be a string(str)')
        if kwargs['RESP_TYPE'].upper() == '[E]':
            return cls(**kwargs)
        for c in cls.__subclasses__():
            obj = c.create(**kwargs)
            if obj: return obj
        return None
    
class CR_Esep_sym(CR_Omni):
    """
    CR_Esep_sym
    energy-separable, phi-symmetric channel response
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def create(cls,**kwargs):
        if 'RESP_TYPE' not in kwargs: raise RFL_KeywordEx('RESP_TYPE required, absent')
        if not isinstance(kwargs['RESP_TYPE'],str): raise RFL_KeywordEx('RESP_TYPE must be a string(str)')
        if '[E]' not in kwargs['RESP_TYPE'].upper(): return None # not Esep
        for c in cls.__subclasses__():
            obj = c.create(**kwargs)
            if obj: return obj
        if kwargs['RESP_TYPE'].upper() == '[E][TH]':
            return cls(**kwargs)
        else:
            return None
    
class CR_Esep_asym(CR_Esep_sym):
    """
    CR_Esep_asym
    energy-separable, phi-asymmetric channel response
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def create(cls,**kwargs):
        if 'RESP_TYPE' not in kwargs: raise RFL_KeywordEx('RESP_TYPE required, absent')
        if not isinstance(kwargs['RESP_TYPE'],str): raise RFL_KeywordEx('RESP_TYPE must be a string(str)')
        for c in cls.__subclasses__():
            obj = c.create(**kwargs)
            if obj: return obj
        if kwargs['RESP_TYPE'].upper() == '[E][TH,PH]':
            return cls(**kwargs)
        else:
            return None
    
class CR_insep_sym(CR_Omni):
    """
    CR_insep_sym
    energy-inseparable, phi-symmetric channel response
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def create(cls,**kwargs):
        if 'RESP_TYPE' not in kwargs: raise RFL_KeywordEx('RESP_TYPE required, absent')
        if not isinstance(kwargs['RESP_TYPE'],str): raise RFL_KeywordEx('RESP_TYPE must be a string(str)')
        for c in cls.__subclasses__():
            obj = c.create(**kwargs)
            if obj: return obj
        if kwargs['RESP_TYPE'].upper() == '[E,TH]':
            return cls(**kwargs)
        else:
            return None
    
class CR_insep_asym(CR_insep_sym):
    """
    CR_insep_asym
    energy-inseparable, phi-asymmetric channel response
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def create(cls,**kwargs):
        if 'RESP_TYPE' not in kwargs: raise RFL_KeywordEx('RESP_TYPE required, absent')
        for c in cls.__subclasses__():
            obj = c.create(**kwargs)
            if obj: return obj
        if kwargs['RESP_TYPE'].upper() == '[E,TH,PH]':
            return cls(**kwargs)
        else:
            return None

if __name__ == '__main__':
    # very simple test
    print(ChannelResponse.create(RESP_TYPE='[E,TH]'))
