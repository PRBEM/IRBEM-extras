"""
rpp.py - right parallel piped calculations

Public Globals:
    VERBOSE - report new cache computations

Public Functions:
    weibull_diff_cross_section - differential Weibull cross section function
    weibull_integral_cross_section - cumulative integral Weibull cross section function
    
Public Classes:
    CrossSection - abstract base class for part cross sections
        WeibullCrossSection - cross section handler for Weibull fits
    EEEPart - abstract base class for SEE-susceptible parts
        ProtonPart - computes SEE rate for proton-sensitive part
        IonPart - virtual base class for ion parts
            IonPartRPP - computes SEE rate for ion-sensitive part using RPP approximation
    RPPChordLengthDist - abstract base class for RPP chord length distributions
        BendelCLD - Bend RPP CLD
        LB88 - Luke and Buehler CLD

Glossary:
    slim - limiting cross section per defice or per bit (at E=inf), units: cm^2, scalar float
        Note slim has units of cm^2 and describes the saturating cross section at 
            normal incidence for a single sensitive volume. 
            **For devices with multiple sensitive volumes, slim is the cross section for a 
            single volume (e.g., one bit in a memory device)**
    L0 - LET threshold for Weibull, units MeV/(g/cm2), scalar float
    E0 - Energy threshold for Weibull, units MeV, scalar float
    W - scale LET or energy for Weibull, units MeV/(g/cm2) or MeV, scalar float
    S - shape parameter for Weibull, dimensionless, scalar float
    h - thickness (um) of part's sensitive volume
    a - middle dimension (um) of part's sensitive volume
    b - long dimension (um) of part's sensitive volume
    aspect - aspect ratio of middle to longest dimension of RPP (a/b, h < a < b, so a/b<1)    
    volumes - number of sensitive volumes (e.g., number of bits)
    default_fast - default run mode is fast (computing/storing matrices)
    epsrel: epsrel argument passed to quad, dblquad numerical integration (relative precision)
    ecludeH - used to exclude protons from LET spectrum and SEE rate for ion parts
        fixed at False for proton parts, defaults to True for ion parts, but can override
    
    
Notes:
    Integration is done by adaptive quadrature using scipy.integrate.quad and 
    dblquad to relative default precision of 1e-4. Unless fast is specified, 
    in which case trapezoidal integrals are used.
"""

import numpy as np
import inspect
from scipy.integrate import quad,dblquad
from degraded_spectra import smart_interp1d
from interpolating_matrix import trapz_dx

VERBOSE = False # print messages each time a new cache is created

def weibull_int_cross_section(E0,W,S,slim,E=None,return_function=False):
    """
    cs = weibull_int_cross_section(E0,W,S,slim,E=None,return_function=False)
    cumulative integral weibull cross section, in units of slim
    for E0,W,S,slim see glossary in module doc string
    E - energy at which to evaluate cross section, MeV scalar or array of float
    if return_fuction is True then a function handle is returned such that
        cs(E) = weibull_int_cross_section(E0,W,S,slim,E)
    cs has units of cm^2 (slim's units)
    """
    if return_function:
        return lambda L: weibull_int_cross_section(E0,W,S,slim,L)
    if E is None:
        raise Exception('Cannot evaluate weibull_int_cross_section without numeric E input')
    isScalar = np.isscalar(E)
    E = np.atleast_1d(E)
    s = np.zeros(E.shape)
    i = E>=E0
    s[i] = slim*(1-np.exp(-((E[i]-E0)/W)**S))
    if isScalar:
        s = s[0]
    return s

def weibull_diff_cross_section(L0,W,S,slim,L=None,return_function=False):
    """
    cs = weibull_diff_cross_section(L0,W,S,slim,L=None,return_function=False)
    differential weibull cross section, in units of slim
    for L0,W,S,slim see glossary in module doc string
    L - LET at which to evaluate cross section, MeV/(g/cm2) scalar or array of float
    if return_fuction is True then a function handle is returned such that
        cs(L) = weibull_diff_cross_section(L0,W,S,slim,L)
    cs has units of cm^2/(MeV/(g/cm2))
    For protons, L0->E0, L->E, MeV/(g/cm2) -> MeV
    """
    if return_function:
        return lambda L: weibull_diff_cross_section(L0,W,S,slim,L)
    if L is None:
        raise Exception('Cannot evaluate weibull_diff_cross_section without numeric L input')
    isScalar = np.isscalar(L)
    L = np.atleast_1d(L)
    s = np.zeros(L.shape)
    i = L>=L0
    # cdf = slim*(1-np.exp(-((L-L0)/W)**S))
    s[i] = slim*np.exp(-((L[i]-L0)/W)**S)*S/W*((L[i]-L0)/W)**(S-1)
    if isScalar:
        s = s[0]
    return s

class CrossSection(object):
    """
    cs = CrossSection(minX,slim)
    Public Methods:
        diff_cross_section - differential cross section
        int_cross_Section - cumulative integral cross section
    Public Data:
        minX - minium value of x for which there's any cross section
        slim - limiting cross section
    """
    def __init__(self,minX,slim):
        self.minX = minX
        self.slim = slim
    def diff_cross_section(self,x):
        """ cs = diff_cross_section(x)
        returns differential cross section at x
        abstract method - must be overloaded
        """
        raise NotImplementedError # virtual base class
    def int_cross_section(self,x):
        """ cs = int_cross_section(x)
        returns integral cross section less than or equal to x
        abstract method - must be overloaded
        """
        raise NotImplementedError  # virtual base class

class WeibullCrossSection(CrossSection):
    """
    cs = WeibullCrossSection(X0,W,S,slim)
    X0 is E0 or L0 (threshold)
    E0,L0,W,S,slim - see glossary in module doc string
    Public Data:
        same as CrossSection plus  X0,W,S,
    Public Methods:
        same as CrossSection
    """
    def __init__(self,X0,W,S,slim):
        super().__init__(X0,slim)
        self.X0 = X0
        self.W = W
        self.S = S
    def __str__(self):
        return 'Webuill with threshold %g, W=%g, S=%g, slim=%g' % (self.X0,self.W,self.S,self.slim)
    def diff_cross_section(self,x):
        """ cs = diff_cross_section(x)
        returns differential cross section at x
        """
        return weibull_diff_cross_section(self.X0,self.W,self.S,slim=self.slim,L=x)
    def int_cross_section(self,x):
        """ cs = int_cross_section(x)
        returns integral cross section less than or equal to x
        """
        return weibull_int_cross_section(self.X0,self.W,self.S,slim=self.slim,E=x)

class EEEPart(object):
    """ EEEpart - abstract base class for ProtonPart and IonPart, etc
    EEEPart(cs,is_proton_part=False,volumes=1,default_fast=True,epsrel=1e-4,excludeH=None) 
    cs CrossSection object (e.g., WeibullCrossSection)
    is_proton_part - is this a proton-susceptible part? (affecs behavior of other functions)
    excludeH - True to exclude protons from LET spectrum 
        applies to direct ionization calculation
        defaults to (not is_proton_part)
    for other arguments see glossary in module docstring
    Public Data:
        is_proton_part, default_fast,epsrel as defined in glossary in module doc string
        cs - CrossSection object (e.g., WeibullCrossSection)
    Public Methods:
        see_rate - compute see rate
        get_cache - retrieve/compute precomputed integration weights
    """
    def __init__(self,cs,is_proton_part=False,volumes=1,default_fast=True,epsrel=1e-4,excludeH=None):
        self.cs = cs
        self.is_proton_part = is_proton_part
        self.default_fast = default_fast
        self.volumes = volumes
        self.epsrel = 1e-4
        if excludeH is None:
            excludeH = not is_proton_part
        self.excludeH = excludeH
        self._cache = None
    def see_rate(self):
        """
        r = part.see_rate(...) - upset rate SEE/s (including multiple volumes)
        abstract method - must be overloaded
        """
        raise NotImplementedError
    def get_cache(self):
        """
        cache = part.get_cache(...) - get cache for rapidly computing SEE rate
        abstract method - must be overloaded
        """
        raise NotImplementedError
        
        
class ProtonPart(EEEPart):
    """
    ProtonPart - object to handle proton part SEE calculations
    part = ProtonPart(cs=None,E0=None,W=None,S=None,slim=None,volumes=1,default_fast=True,epsrel=1e-4)
    Supply either cs or E0,W,S,slim (assumes Weibull)
    cs - CrossSection object (e.g., WeibullCrossSection)
    for other inputs, see glossary module doc string
    Public Data:
        same as EEEPart
    Public Methods:
        same as EEEPart
    """
    def __init__(self,cs=None,E0=None,W=None,S=None,slim=None,volumes=1,default_fast=True,epsrel=1e-4):
        if cs is None:
            cs = WeibullCrossSection(E0,W,S,slim)
        super().__init__(cs,is_proton_part=True,volumes=volumes,default_fast=default_fast,epsrel=epsrel)

    def __str__(self):
        return 'Proton Part with E0=%g, W=%g MeV, S=%g, slim=%g cm^2 x %d volumes' % (
                self.cs.X0,self.cs.W,self.cs.S,self.cs.slim,self.volumes)
    def see_rate(self,flux,energy=None,energyRange=None,dEoverE=0.1,fast=None):
        """
        r = see_rate(flux,energy=None,energyRange=None,dEoverE=0.1,fast=None)
        r - upset rate SEE/s (including multiple volumes)
        flux - two options:
            function handle flux(Energy) returns differential flux vs energy. units #/cm2/s/MeV
            array differential flux vs Energy in #/cm2/s/MeV, requires energy input
        energy - array of energies at which flux is specified. MeV
        energyRange - option tuple (minE,maxE) gives domain of flux callable
        dEoverE: relative spacing of working energy grid when generated from energyRange
        fast: build/store/use precomputed matrices
            if fast is None, use self.default_fast
            matrices will be recomputed any time a new energy grid provided or energyRange/dEoverE specified
            But if same energy grid is provided/specified, integration weights will be re-used
        """
        
        if fast is None:
            fast = self.default_fast
            
        if fast:
            cache = self.get_cache(energy,energyRange,dEoverE)
            if callable(flux):
                flux = flux(cache['energy'])
            return np.dot(cache['g'],flux)
    
        if energy is not None:
            energyRange = (energy[0],energy[-1])
    
        if energyRange is None:
            energyRange = (0,np.inf)
    
        if callable(flux):
            flux_func = flux
        else:
            # flux is array, LET is provided, create interpolating function, zero outside limits
            if energy is None:
                raise Exception('If flux is not a callable, energy array must be provided')
            flux = np.array(flux)
            energy = np.array(energy)
            if len(flux) != len(energy):
                raise Exception('flux and energy must be same length')
            flux_func = lambda E: smart_interp1d(np.log(energy),flux,np.log(E),'linear',copy=False,assume_sorted=True,bounds_error=False,fill_value=0)
    
        # will integrate log(E)
        minE = np.maximum(energyRange[0],self.cs.minX) # do not integrate below threshold
        maxE = energyRange[1]
    
        func = lambda logE : self.cs.int_cross_section(np.exp(logE))*flux_func(np.exp(logE))*np.exp(logE) # func(y,x), extra E factor due to dlogE
        r = quad(func, np.log(minE), np.log(maxE),epsrel=self.epsrel)[0]
        return r*self.volumes
    
    def get_cache(self,energy=None,energyRange=None,dEoverE=0.1):
        """
        cache = part.get_cache(energy=None,energyRange=None,dEoverE=0.1)
        see part.see_rate for description of arguments
        cache is a dict with fields:
            energy: the energy grid provided by energy or specified by energyRange/dEoverE
            g: the weights for numerical integration: rate = g*flux
        """
        
        if energy is None:
            if energyRange is None:
                raise Exception('energy or energyRange must be specified')
            NE = np.maximum(100,int(np.ceil(25*0.1/dEoverE*np.log10(energyRange[1]/energyRange[0]))))
            energy = np.exp(np.linspace(np.log(energyRange[0]),np.log(energyRange[1]),NE))
            
        energy = energy.ravel() # (NE,)
        if self._cache is not None:
            if not np.array_equal(self._cache['energy'],energy):
                self._cache = None
        if self._cache is None:
            if VERBOSE:
                print('New g cache for %s' % self)
            self._cache = {}
            self._cache['energy'] = energy
            dE = trapz_dx(energy)
            g = self.cs.int_cross_section(energy)*dE
            self._cache['g'] = g*self.volumes
    
        return self._cache

class IonPart(EEEPart):
    """
    virtual base class for ion parts
        part = IonPart(cs,Ap,xnorm,xmax,volumes=1,default_fast=True,epsrel=1e-4,excludeH=True)
        cs CrossSection object (e.g., WeibullCrossSection)
        Ap - projected area accounting for all directions (cm^2)
        xnorm - effective chord length through sensitive volume at normal incidence (cm)
        xmax - maximum chord length through sensitive volume (cm)
        ecludeH - exclude protons from LET spectrum and SEE rate
        other args same meaning as EEEPart
        Public Data:
            same as EEEPart plus...
            Ap, xnorm, xmax
        Public Methods:
            same as EEEPart
    """
    def __init__(self,cs,Ap,xnorm,xmax,volumes=1,default_fast=True,epsrel=1e-4,excludeH=True):
        super().__init__(cs,is_proton_part=False,volumes=volumes,
             default_fast=default_fast,epsrel=epsrel,excludeH=excludeH)
        self.cs = cs
        self.Ap = Ap
        self.xnorm = xnorm
        self.xmax = xmax
        self._cache = None
    def chord_length_ccdf(self,x):
        """
        ccdf = chord_length_ccdf(x)
        complementary cumulative chord length distribution
        x - point at which to evaluate ccdf in cm
        abstract method - must be overloaded
        """
        raise NotImplementedError
    def see_rate(self,flux,LET=None,LETrange=None,dLoverL=0.1,fast=None):
        """
        r = see_rate_ions(flux,LET=None,LETrange=None,dLoverL=0.1,fast=None)
        r - upset rate SEE/s (including multiple volumes)
        flux - two options:
            function handle flux(LET) returns differential flux vs LET. units #/cm2/s/(MeV/(g/cm2))
            array differential flux vs LET in #/cm2/s/(MeV/(g/cm2)), requires LET input
        for L0,W,S,slim,h,a,b,aspect see glossary in module docstring
        if slim is provided then also provide aspect
        if slim is not provided, provide a and b
        The sides of RPP in um must satisfy h < a < b
        LET - array of LETs at which flux is specified. MeV/(g/cm2)
        LETrange - option tuple (minLET,maxLET) gives domain of flux callable
        dLoverL: relative spacing of working LET grid when generated from LETrange
        Note flux, L0,W,LET, and LETrange nominally have LET in units of MeV/(g/cm2)
          but the code also works if all LET are in MeV/(mg/cm2)
        fast: build/store/use precomputed matrices
            if fast is None, use self.default_fast
            matrices will be recomputed any time a new energy grid provided or energyRange/dEoverE specified
            But if same energy grid is provided/specified, integration weights will be re-used
        """
        
        if fast is None:
            fast = self.default_fast
            
        if fast:
            cache = self.get_cache(LET,LETrange,dLoverL)
            if callable(flux):
                flux = flux(cache['LET'])
            return np.dot(cache['g'],flux)
    
        if LET is not None:
            LETrange = (LET[0],LET[-1])
    
        if LETrange is None:
            LETrange = (0,np.inf)
    
        if callable(flux):
            flux_func = flux
        else:
            # flux is array, LET is provided, create interpolating function, zero outside limits
            if LET is None:
                raise Exception('If flux is not a callable, LET array must be provided')
            flux = np.array(flux)
            LET = np.array(LET)
            if len(flux) != len(LET):
                raise Exception('flux and LET must be same length')
            flux_func = lambda L: smart_interp1d(np.log(LET),flux,np.log(L),'linear',copy=False,assume_sorted=True,bounds_error=False,fill_value=0)
    
            
        # will integrate log(LET)
        minLET = LETrange[0]
        maxLET = LETrange[1]
        
        # integrate dlogLw dlogL
        logL0 = np.log(self.cs.minX)
        logLmin = np.maximum(np.log(minLET),np.log(self.xnorm/self.xmax)+logL0)
        logLmax = np.log(maxLET)
        logLmin = np.minimum(logLmin,logLmax) # insure integral is 0 when xnorm/xmax*L0 > maxLET
        logLwMin = lambda logL : np.maximum(logL,logL0)
        logLwMax = lambda logL : np.maximum(np.log(self.xmax/self.xnorm)+logL,logL0)
        
        
        # dblquad: int func(y,x) dy dx
        # y = dlogLw, x = logL
        # L*Lw factors due to dlogL dlogLw
        func = lambda logLw,logL : self.Ap/self.cs.slim*self.cs.diff_cross_section(np.exp(logLw))*flux_func(np.exp(logL))*self.chord_length_ccdf(np.exp(logLw-logL)*self.xnorm)*np.exp(logL+logLw) 
        r = dblquad(func, logLmin,logLmax,logLwMin,logLwMax,epsrel=self.epsrel)[0]
        
        return r*self.volumes

    def get_cache(self,LET=None,LETrange=None,dLoverL=0.1):
        """
        cache = part.get_cache(LET=None,LETrange=None,dLoverL=0.1)
        see part.see_rate for description of arguments
        cache is a dict with fields:
            LET: the LET grid provided by energy or specified by LETrange/dLoverL
            g: the weights for numerical integration: rate = g*flux
        """
        
        if LET is None:
            if LETrange is None:
                raise Exception('LET or LETrange must be specified')
            NL = np.maximum(100,int(np.ceil(25*0.1/dLoverL*np.log10(LETrange[1]/LETrange[0]))))
            LET = np.exp(np.linspace(np.log(LETrange[0]),np.log(LETrange[1]),NL))
        else:
            LET = np.array(LET).ravel() # force 1-D
            
        if self._cache is not None:
            if not np.array_equal(self._cache['LET'],LET):
                self._cache = None
        if self._cache is None:
            if VERBOSE:
                print('New g cache for %s' % self)
            self._cache = {}
            self._cache['LET'] = LET

            dLET = trapz_dx(LET)
            g = np.zeros(dLET.shape)
            
            # integrate dlogLw, hence extra exp(logLw) factor
            func = lambda logLw: self.cs.diff_cross_section(np.exp(logLw))*self.chord_length_ccdf(np.exp(logLw)/L*self.xnorm)*np.exp(logLw)
            for (i,L) in enumerate(LET):
                LET1 = np.maximum(L,self.cs.minX)
                LET2 = np.maximum(L*self.xmax/self.xnorm,self.cs.minX)
                g[i] = self.Ap/self.cs.slim*quad(func,np.log(LET1),np.log(LET2),epsrel=self.epsrel)[0]*dLET[i]
            self._cache['g'] = g*self.volumes
    
        return self._cache

class RPPChordLengthDist(object):
    """
    Virtual Base Class for RPP chord length calculations
    cld = RPPChordLengthDist(a,b,c)
    a, b, c dimensions in cm of RPP
    a <= b <= c
    Public Data:
        a,b,c
        slim - limiting cross section in cm^2
        xmax- longest chord length
    Public Methods:
        chord_length_ccdf
    """
    
    def __init__(self,a,b,c):
        # sort the inputs
        (a,b,c) = sorted((a,b,c)) # ensure a<=b<=c
        self.a = a
        self.b = b
        self.c = c
        self.slim = b*c
        self.xmax = np.sqrt(a**2 + b**2 + c**2)
    def chord_length_ccdf(self,x):
        """
        ccdf = chord_length_ccdf(x)
        complementary cumulative chord length distribution (RPP)
        x - point at which to evaluate ccdf in cm
        abstract method - must be overloaded
        """
        raise NotImplementedError

class BendelCLD(RPPChordLengthDist):
    """
    cld = BendelCLD(a,b,c)
    a, b, c dimensions in cm of RPP
    a <= b <= c
    computes RPP chord length distribuiton using Bendel 1984 method:
    Bendel, W.L., Length distribution of chords through a rectangular volume,
    NRL Memorandum 5369, Naval Research Laboratory, July 1984.    
    """
    
    def __init__(self,a,b,c):
        super().__init__(a,b,c)
    def __str__(self):
        return 'Bendel RPP CLD with a=%g cm' % self.a
    def chord_length_ccdf(self,x):
        """
        ccdf = chord_length_ccdf(x)
        complementary cumulative chord length distribution (RPP)
        x - point at which to evaluate ccdf in cm
        """
        isScalar = np.isscalar(x)
        x = np.atleast_1d(x)
        ccdf = np.zeros(x.shape)
        i = x<=self.a
        ccdf[i] = 1-x[i]/4/self.a
        i = np.logical_not(i)
        ccdf[i] = 3/4*(self.a/x[i])**2.2
        ccdf[x>=self.xmax] = 0
        if isScalar:
            ccdf = ccdf[0]
        return ccdf
        
    
class LB88(RPPChordLengthDist):
    """
    cld = LB88(a,b,c)
    a, b, c dimensions in cm of RPP
    a <= b <= c
    computes RPP chord length distribuiton using LB88 exact method
    LB88: Luke and Buehler, An exact, closed-form expression of the integral chord-length 
    distribution for the calculation of single-event upsets induced by cosmic rays,
    J. App. Phys. 64, 5132, 1988. doi:10.1063/1.342420.
    """
    
    def __init__(self,a,b,c):
        super().__init__(a,b,c)
        self.V = a*b*c # volume, equation A1
        self.S = 2*(b*c+c*a+a*b) # A2
        self.y = 6*self.V*np.pi # A3
        self.z = a+b+c # A4
        self.r = (b*c)**2 + (c*a)**2 + (a*b)**2 # A5
        self.w = self.xmax # A6
        self.va = np.sqrt(b**2+c**2) # A7
        self.vb = np.sqrt(a**2+c**2) # A8
        self.vc = np.sqrt(a**2+b**2) # A9
        self.x0 = np.maximum(c,self.vc) # A20
        
        # handle possible equation swap when vc>c
        if self.vc <= self.c: # default
            self.C4 = self.C4vc # vc < x <= v
        else: # fc>c: swap vc<->c
            self.C4 = self.C4c  # c < x <= vc
            
    def __str__(self):
        return 'Exact RPP CLD with %g x %g x %g cm' % (self.a,self.b,self.c)
    def s(self,x,gamma):
        """
        Equation A10 of LB88
        s = s(x,gamma) = sqrt(x**2-gamma**2)
        """
        return np.sqrt(x**2-gamma**2)

    def C1(self,x):
        """
        equation 1 from LB88
        C1(x) = 1-F(x;0)
        """
        return 1-self.F0(x)
    def C2(self,x):
        """
        equation 1 from LB88
        C2(x) = C1(a)-F(x;a)
        """
        return self.C1(self.a)-self.Fa(x)
    def C3(self,x):
        """
        equation 1 from LB88
        C3(x) = C2(b)-F(x;b)
        """
        return self.C2(self.b)-self.Fb(x)
    def C4vc(self,x):
        """
        equation 1 from LB88 when vc<=c
        C4(x) = C3(vc)-F(x;vc)
        """
        return self.C3(self.vc)-self.Fvc(x)
    def C4c(self,x):
        """
        equation 1 from LB88 when vc>c
        C4(x) = C3(c)-F(x;c)
        """
        return self.C3(self.c)-self.Fc(x)
    def C5(self,x):
        """
        equation 1 from LB88
        C5(x) = C4(x0)-F(x;x0)
        """
        return self.C4(self.x0)-self.Fx0(x)
    def C6(self,x):
        """
        equation 1 from LB88
        C6(x) = C5(vb)-F(x;vb)
        """
        return self.C5(self.vb)-self.Fvb(x)
    def C7(self,x):
        """
        equation 1 from LB88
        C7(x) = C6(va)-F(x;va) = F(w;va)-F(x;va)
        1st from is from LB88, but is badly affected by round off
        use 2nd form with better floating point numerics
        """
        return self.Fva(self.w)-self.Fva(x)
    def G3(self,x,gamma):
        """
        equation A13 from LB88
        G3(x,gamma)
        """
        s = self.s(x,gamma)
        return 8*s-2*gamma**2*s/x**2-6*gamma*np.arctan(s/gamma)
    def T3(self,x,a=None,b=None):
        """
        equation A14 from LB88
        T3 = T3(x,a,b)
        """
        # here vc is replaced by vb = sqrt(h**2+a**2)
        if (a is None) and (b is None):
            vc = self.vc # defaults, already know vc
        else:
            vc = None # compute it after getting a,b
        
        if a is None:
            a = self.a
        if b is None:
            b = self.b
        if vc is None:
            vc = np.sqrt(a**2+b**2)            
            
        sqrt_x2mvc2 = np.sqrt(x**2-vc**2) # save for use 3 times
        return 12*self.V*(
                (-a/2/x**2+1/2/a)*np.arctan(sqrt_x2mvc2/b)
                +(-b/2/x**2+1/2/b)*np.arctan(sqrt_x2mvc2/a)
                -(1/2/vc)*(b/a+a/b)*np.arctan(sqrt_x2mvc2/vc)
                          )
    def F0(self,x):
        """
        equation A15 from LB88
        F(x;0)
        """
        assert (0<=x) and (x<=self.a),'F0: 0<=x<=a (0<=%g<=%g)' % (x,self.a)
        return 2/3/np.pi/self.S*(
                -9*x**2/2 + 8*self.z*x
                )
    def Fa(self,x):
        """
        equation A16 from LB88
        F(x;a)
        """
        assert (self.a<=x) and (x<=self.b),'Fa: a<x<=b'
        return 2/3/np.pi/self.S*(
                (self.a*self.y - self.a**4)/2*(1/self.a**2 - 1/x**2) 
                + 8*(self.b + self.c)*(x - self.a)
                - (self.b + self.c)*self.G3(x,self.a))
    def Fb(self,x):
        """
        equation A17 from LB88
        F(x;b)
        """
        assert (self.b<=x) and (x<=np.minimum(self.c,self.vc)),'Fb: b<x<=min(c,vc) (%g<=%g<=%g)' % (self.b,x,np.minimum(self.c,self.vc))
        return 2/3/np.pi/self.S*(
                ((self.a + self.b)*self.y-(self.a**4 + self.b**4))/2*(1/self.b**2 - 1/x**2)
                +9*(x**2 - self.b**2)/2 + 8*self.c*(x - self.b)
                -(self.b + self.c)*(self.G3(x,self.a)-self.G3(self.b,self.a))
                -(self.c + self.a)*self.G3(x,self.b)                
                )
    def Fvc(self,x):
        """
        equation A18 from LB88
        F(x;vc) 
        """
        assert (self.vc<=x) and (x<=self.c),'Fvc: vc<x<=c (%g<%g<=%g)' % (self.vc,x,self.c)
        return 2/3/np.pi/self.S*(
                ((self.a+self.b)*self.y - 6*self.a**2*self.b**2)/2*(1/self.vc**2 - 1/x**2)
                +8*self.c*(x - self.vc) + self.c*self.G3(x,self.vc)
                -self.c*(self.G3(x,self.a)-self.G3(self.vc,self.a))
                -self.c*(self.G3(x,self.b)-self.G3(self.vc,self.b))
                -self.T3(x)
                )
    def Fc(self,x):
        """
        equation A19 from LB88
        F(x;c)
        """
        assert (self.c<=x) and (x<=self.vc),'Fc: c<x<=vc (%g<%g<=%g)' % (self.c,x,self.vc)
        return 2/3/np.pi/self.S*(
                (self.z*self.y - (self.a**4+self.b**4+self.c**4))/2*(1/self.c**2 - 1/x**2)
                +9*(x**2 - self.c**2)
                -(self.b+self.c)*(self.G3(x,self.a)-self.G3(self.c,self.a))
                -(self.c+self.a)*(self.G3(x,self.b)-self.G3(self.c,self.b))
                -(self.a+self.b)*self.G3(x,self.c)
                )
    def Fx0(self,x):
        """
        equation A20 from LB88
        F(x;x0) 
        """
        assert (self.x0<=x) and (x<=self.vb),'Fx0: x0<x<=vb (%g<%g<=%g)' % (self.x0,x,self.vb)
        return 2/3/np.pi/self.S*(
                (self.z*self.y - self.c**4 - 6*self.a**2*self.b**2)/2*(1/self.x0**2 - 1/x**2)
                +9*(x**2 - self.x0**2)/2
                +self.c*(self.G3(x,self.vc)-self.G3(self.x0,self.vc))
                -self.c*(self.G3(x,self.a)-self.G3(self.x0,self.a))
                -self.c*(self.G3(x,self.b)-self.G3(self.x0,self.b))
                -(self.a+self.b)*(self.G3(x,self.c)-self.G3(self.x0,self.c))
                -(self.T3(x)-self.T3(self.x0))
                )
    def Fvb(self,x):
        """
        equation A21 from LB88
        F(x;vb)
        """
        assert (self.vb<=x) and (x<=self.va),'Fvb: vb<x<=va (%g<%g<=%g)' % (self.vb,x,self.va)
        return 2/3/np.pi/self.S*(
                (self.z*self.y + self.a**4 - 6*self.a**2*self.va**2)/2*(1/self.vb**2 - 1/x**2)
                +self.b*self.G3(x,self.vb)
                +self.c*(self.G3(x,self.vc)-self.G3(self.vb,self.vc))
                -self.b*(self.G3(x,self.c)-self.G3(self.vb,self.c))
                -self.c*(self.G3(x,self.b)-self.G3(self.vb,self.b))
                -(self.T3(x)-self.T3(self.vb))
                -self.T3(x,self.a,self.c)
                )
    def Fva(self,x):
        """
        equation A22 from LB88
        F(x;va)
        NOTE! the 9(x^2-va^2)/2 term is positive in LB88 but 
        here it is negative. This change was required to get the 
        ccdf not to go negative (i.e., it makes it go to exactly 0 at w)
        """
        assert (self.va<=x) and (x<=self.w),'Fva: va<x<=w (%g<%g<=%g)' % (self.va,x,self.w)
        return 2/3/np.pi/self.S*(
                (self.z*self.y + self.w**4 - 8*self.r)/2*(1/self.va**2 - 1/x**2)
                -9*(x**2 - self.va**2)/2 # sign reversed from LB88!!!!
                +self.a*self.G3(x,self.va)
                +self.b*(self.G3(x,self.vb)-self.G3(self.va,self.vb))
                +self.c*(self.G3(x,self.vc)-self.G3(self.va,self.vc))
                -self.T3(x,self.b,self.c)
                -(self.T3(x,self.c,self.a)-self.T3(self.va,self.c,self.a))
                -(self.T3(x)-self.T3(self.va))
                )
    def chord_length_ccdf(self,x):
        """
        ccdf = chord_length_ccdf(x)
        compelmentary cumulative chord length distribution (RPP)
        x - point at which to evaluate ccdf in cm
        """
        if not np.isscalar(x):
            return np.vectorize(self.chord_length_ccdf)(x)
        
        # equation 1 of LB88
        C = 0
        if x <= self.a:
            C = self.C1(x)
        elif x <= self.b:
            C = self.C2(x)
        elif x <= np.minimum(self.c,self.vc):
            C = self.C3(x)
        elif x <= self.x0: # max(c,vc)
            C = self.C4(x)
        elif x <= self.vb:
            C = self.C5(x)
        elif x <= self.va:
            C =self.C6(x) # problem is here
        elif x <= self.w:
            C = self.C7(x)
        if C < -10*np.finfo(float).eps: # check for *significantly* negative C
            raise Exception('chord_length_ccdf %g < 0 at x=%g' % (C,x))
        else:
            C = np.maximum(C,0) # correct for cumulative round off errors
            
        return C
    
class IonPartRPP(IonPart):
    """
    part = IonPartRPP(h,cs=None,L0=None,W=None,S=None,slim=None,aspect=1,a=None,b=None,cld=None,volumes=1,default_fast=True,epsrel=1e-4,excludeH=True)

    Ion part modeled as right parallel piped (RPP) using chord length distribution (default LB88 method)
    
    if cs is provided, it replaces L0, W, S, slim (supply aspect ratio)
    aspect - aspect ratio of longest to middle dimension of RPP (b/a)
    for other arguments see glossary in module doc string
    if slim is provided then aspect is used to determine a, b
    if slim is not provided, provide a and b
    cld: RPP chord length distribution (if supplied no need for a,b,aspect,slim)
    The sides of RPP in um must satisfy h < a < b
    Public Data:
        same as IonPart plus...
        h,a,b as in glossary in module doc string
        hcm, acm, bcm - h,a,b in cm
        
    Public Methods:
        same as IonPart
    """
    def __init__(self,h,cs=None,L0=None,W=None,S=None,slim=None,aspect=1,a=None,b=None,cld=None,volumes=1,default_fast=True,epsrel=1e-4,excludeH=True):
        
        if cs is not None:
            slim = cs.slim            

        if cld is not None:
            if slim is None:
                slim = cld.slim            

        # figure out if slim or a,b given and create h,a,b variants in cm
        if slim is None:
            if (a is None) or (b is None):
                raise Exception('If slim not provided, a and b must be')
            if a > b:
                (a,b) = (b,a) # swap to ensure h < a < b
            acm = a/1e4 # um -> cm
            bcm = b/1e4 # um -> cm
            slim = acm*bcm
            aspect = 1
        else:
            if aspect>1:
                aspect = 1/aspect # user entered aspect ratio upside down
            acm = np.sqrt(slim*aspect) # a < b for aspect < 1
            bcm = np.sqrt(slim/aspect)
            a = acm*1e4
            b = bcm*1e4
            # acm*bcm = slim
            # bcm/acm = aspect
        # above logic ensures acm <= bcm

        hcm = h/1e4 # um->cm    
        if hcm>acm:
            raise Exception('Thickness is larger than 2nd dimension (a). (h < a < b) or (h^2 < slim) required')
            
        if cld is None:
            cld = LB88(hcm,acm,bcm)
    
        if cs is None:
            cs = WeibullCrossSection(L0,W,S,slim)

        assert np.abs(cld.slim-cs.slim)<(cld.slim+cs.slim)/2/1000,'inconsistent slim cld->%g,cs->%g' % (cld.slim,cs.slim)

        Ap = (hcm*acm + hcm*bcm + acm*bcm)/2 # projected surface area
        xmax = np.sqrt(hcm**2 + acm**2 + bcm**2) # maximum path length
        # xnorm == hcm
        super().__init__(cs,Ap,hcm,xmax,volumes=volumes,default_fast=default_fast,epsrel=epsrel,excludeH=excludeH)
        
        self.cld = cld
        self.h = h
        self.hcm = hcm
        self.a = a
        self.acm = acm
        self.b = b
        self.bcm = bcm
        
    def __str__(self):
        return 'Ion Part with L0=%g, W=%g MeV/(g/cm^2), S=%g, RPP: %g x %g x %g um, x %d volumes' % (
                self.cs.X0,self.cs.W,self.cs.S,self.h,self.a,self.b,self.volumes)
        
    def chord_length_ccdf(self,x):
        return self.cld.chord_length_ccdf(x)


def dict_to_part(part,protonClass=ProtonPart,ionClass=IonPartRPP):
    """
    Part = dict_to_part(part,protonClass=ProtonPart,ionClass=IonPartRPP)
    part is a dict produced by e.g., creme/read_report
    Part is an EEEPart
    part has fields:
        is_proton_part - True for proton part
        E0 or L0, W, S, slim, volumes - Weibull parameters as defined in glossary in module doc string
        for ion parts: h, aspect        
    protonClass - class used for proton parts
    ionClass - class used for ion parts        
    """
    
    is_proton_part = ('is_proton_part' in part) and part['is_proton_part']
    if is_proton_part:
        c = protonClass
    else:
        c = ionClass
        
    # build arguments by keeping only keys in parts that are arguments to class
    args = {}
    for (k,info) in inspect.signature(c).parameters.items():
        if info.default == info.empty:
            continue # not a keyword arg
        if k in part:
            args[k] = part[k]

    if is_proton_part:
        return protonClass(**args)
    else:
        return ionClass(part['h'],**args)
    

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import os
    from creme import read_report
    
    _TOP_PATH = os.path.dirname(os.path.abspath(__file__))    
    
    data_path = os.path.join(_TOP_PATH,'data')

    # prepare figure path
    figpath = os.path.join(_TOP_PATH,'figures')
    if not os.path.exists(figpath):
        os.mkdir(figpath)
    
    pups = read_report(os.path.join(data_path,'q800a90i100.PUP'))
    hups = read_report(os.path.join(data_path,'q800a90i100Z2.HUP'))
    
    # make plot of chord length dist, raw weibull, rpp cross section
    # examples from LB88: fig 1 and fig 3
    for (ifig,(a,b,c)) in enumerate([(1,4,6),(1,4,40)]):
        # a,b,c in um even though normally use cm for CLDs
        bendel = BendelCLD(a,b,c)
        lb88 = LB88(a,b,c)
        plt.figure()
        x = np.linspace(0,np.sqrt(a**2+b**2+c**2)*1.2,10000)
        plt.semilogy(x,bendel.chord_length_ccdf(x),'k-',x,lb88.chord_length_ccdf(x),'b-')
        plt.xlabel('Chord Length, x (um)')
        plt.ylabel('Complementary CDF, C(x)')
        plt.ylim([1e-7,1])
        plt.legend(['Bendel','Luke and Buehler'])
        plt.savefig(os.path.join(figpath,'LB88_fig%d.png' % (ifig+1)))

    # test cases, protons
    for partd in pups.values():
        part = dict_to_part(partd)
        Egrid = np.exp(np.linspace(np.log(1e-1),np.log(1e5),1002))
        flux = np.exp(-Egrid/10)
        rslow = part.see_rate(flux,energy=Egrid,fast=False)
        rfast = part.see_rate(flux,energy=Egrid,fast=True)
        cs = WeibullCrossSection(partd['E0'],partd['W'],partd['S'],partd['slim'])
        part2 = ProtonPart(cs=cs)
        rslow2 = part2.see_rate(flux,energy=Egrid,fast=False)
        rfast2 = part2.see_rate(flux,energy=Egrid,fast=True)
        print('part: %s' % part)
        print('proton SEE rate slow %g(%g), fast = %g(%g)' % (rslow,rslow2,rfast,rfast2))
    
    # test cases, ions
    for partd in hups.values():
        part = dict_to_part(partd)
        LETgrid = np.exp(np.linspace(np.log(1),np.log(110000),1002))
        fluxLET = (LETgrid/20)**-3
        rslow = part.see_rate(fluxLET,LET=LETgrid,fast=False)
        rfast = part.see_rate(fluxLET,LET=LETgrid,fast=True)
        cs = WeibullCrossSection(partd['L0'],partd['W'],partd['S'],partd['slim'])
        part2 = IonPartRPP(partd['h'],cs=cs)
        rslow2 = part2.see_rate(fluxLET,LET=LETgrid,fast=False)
        rfast2 = part2.see_rate(fluxLET,LET=LETgrid,fast=True)
        print('part: %s' % part)
        print('ion SEE rate slow %g(%g), fast = %g(%g)' % (rslow,rslow2,rfast,rfast2))
    
    plt.show()
    
