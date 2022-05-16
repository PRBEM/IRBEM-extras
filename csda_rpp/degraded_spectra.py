""" degraded_spectra.py
Routines for computing degraded proton or electron spectra for slab geometries
Continuous Slowing Down Approximation (CSDA) only
Paul O'Brien

v1.0 - single-species tables (multiple materials)
v2.0 - add support for multiple species in same table

Note: while CSDA is a very good approximation for protons, it is not 
nearly as good for electrons. Electrons scatter significantly.

SCREAM proton and electron range tables provided by Dr. Scott Messenger. SCREAM
proton range data come from SRIM and electron range data from NIST/ESTAR.

NIST/PSTAR proton range tables and SRIM ion range tables are also provided.

Public Globals:
    VERBOSE - report new cache computations

Public functions:
    locate_file - find file in working directory or module folder
    read_elements - read table of elements
    standardize_species - standardize a species string (e.g., 'e-' to 'e')
    tocm - convert length unit to centimeters
    smart_interp1d - 1-d inerpolation that tries in log y and fails over to linear in y

Public Classes:
    RangeTable - abstract base class for range tables
        ScreamTable - range table class for SCREAM and NIST data (single species, multiple materials)
        SRIMTable - range table class for SRIM data (multiple species, single material)

Notes:
    
By default, LET spectra omit Z=1 because it takes strong evidence to conclude a part
is sentivie to direct ionization by protons. And otherwise, the 
constant-LET-along-path assumption in the RPP integrals can break for protons.
    
Interpolation is done with power-law first, falling over to semi-log when zeros are present

Differentiation is done in a power law sense: e.g., dR/dE = dlogR/dlogE*R/E

For energies beyond range table (1 GeV for SCREAM, 10 GeV for NIST, 10 MeV/nuc for SRIM), 
particle is assumed to traverse any depth of shielding with negligible 
energy loss.

Integration is done by adaptive quadrature using scipy.integrate.quad and 
dblquad to relative default precision of 1e-4. Unless fast is specified, 
in which case trapezoidal integrals are used.


"""

import os
import re
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
from interpolating_matrix import interpolating_matrix,trapz_dx

import warnings
warnings.filterwarnings('ignore','invalid value encountered in')
warnings.filterwarnings('ignore','divide by zero encountered in')

global _HIGH_ENERGY_WARNED # have we already warned about the high energy limit?
_HIGH_ENERGY_WARNED = False


VERBOSE = False # report new cache computations

_TOP_PATH = os.path.dirname(os.path.abspath(__file__))
def locate_file(filename):
    """
    filename = locate_file(filename)
    if filename is not present in working directory
        look for it in the folder containing this module file
        and its data subdirectory
    return the filename (unchanged if not found)
    """    
    if not os.path.exists(filename):
        testname = os.path.join(_TOP_PATH,filename)
        if os.path.exists(testname):
            return testname
        testname = os.path.join(_TOP_PATH,'data',filename)
        if os.path.exists(testname):
            return testname
    return filename

chemical_elements = {}
def read_elements(filename=None):
    """
    chemical_elements = read_elements(filename=None)
    filename defaults to elements.csv
    sets global chemical_elements
    chemical_elements[symbol].Z, .A
    symbol - e (electron), H(proton/hydrogen), He, Li,...
    Z - atomic number
    A - atomic mass (abundance-weighted average of isotopes)
    MAI - integer mass number of most abundant ion
    if filename=None, and chemical_elements already loaded, has no effect
    returns chemical_elements (which is also a global)
    """
    global chemical_elements
    if filename is None:
        if bool(chemical_elements):
            return chemical_elements # already loaded, don't replace
        filename = 'elements.csv'
    filename = locate_file(filename)
    if not os.path.exists(filename):
        raise Exception('%s does not exist' % filename)
    # Z,A,MAI
    data = np.genfromtxt(filename,delimiter=',',names=True,dtype=None,encoding='UTF-8')
    for (symbol,Z,A,MAI) in data:
        chemical_elements[symbol] = {'A':A,'Z':Z,'MAI':MAI}
    return chemical_elements
        
def standardize_species(species):
    """
    convert various alternatives to standard species name: 'e' for electrons, 'H','He',... for ions
    """
    if species is None:
        raise Exception('Cannot standardize species "None"')
    species0 = species
    species = species.replace('-','')
    species= species.replace('+','')
    if species.lower() in ('e','e-','electron','beta'):
        species = 'e'
    elif species.lower() in ('p+','h+','proton','h','hydrogen'):
        species = 'H'    
    elif species == 'p': # lowercase p is proton. uppercase P is phosphorous
        species = 'H'
    else:
        read_elements()
        if len(species) == 2:
            species = species[0].upper() + species[1].lower() # chemical symbol
        if species not in chemical_elements:
            raise Exception('Unknown species: %s (%s) ' % (str(species0),species))
    return species
    
def tocm(thickness,units,density=None):
    """ cm = tocm(thickness,units,density=None)
    converts thickness to cm
    thickness: numerical thickness
    units - string, one of: mils, in, cm, mm, um, A, g/cm^2, mg/cm^2
    density - numerical density in g/cm^3, only used by g/cm^2 and mg/cm^2 units
    """
    units = units.replace('cm2','cm^2') # standardize
    if units == 'mils':
        return thickness/1000*2.54
    elif units == 'in':
        return thickness*2.54
    elif units == 'cm':
        return thickness
    elif units == 'mm':
        return thickness/10
    elif units == 'um':
        return thickness/1e4
    elif units == 'A': # angstroms
        return thickness/1e5
    elif units == 'g/cm^2':
        if density is not None:
            return thickness/density
        else:
            raise Exception('Density required for units %s' % units)
    elif units == 'mg/cm^2':
        if density is not None:
            return thickness/density*1e-3
        else:
            raise Exception('Density required for units %s' % units)
    else:
        raise Exception('unknown units %s' % units)

def smart_interp1d(x,y,xI,*args,**kwargs):
    """
    yI = smart_interp1d(x,y,xI,*args,**kwargs)
    try to log-interpolate in y ont xI and fail over to linear wherever y=0
    arg and kwargs are passed to interp1d
    fill_value, if specified, is log'd in the log interp, appropriately
    """
    fill_value = np.nan
    if 'fill_value' in kwargs:
        fill_value = kwargs['fill_value']
        del kwargs['fill_value']
    isScalar = np.isscalar(xI)    
    xI = np.atleast_1d(xI)
    yI = np.exp(interp1d(x,np.log(y),*args,**kwargs,fill_value=np.nan)(xI))
    j = np.logical_not(np.isfinite(yI))
    if np.any(j):
        yI[j] = interp1d(x,y,*args,**kwargs,fill_value=fill_value)(xI[j])
    if isScalar:
        yI = yI[0]
    return yI

class RangeTable(object):
    """
    RangeTable abstract base class
        provides methods for accessing and utilizing range-energy data
        table = RangeTable(materials,Species,species,Z2sym)
    Public Methods:
        material_check - check material is in table
        species_check - check species is in table
        EtoR - energy (MeV) to range (g/cm^2)
        RtoE - range (g/cm^2) to energy (MeV)
        dRdE - derivative of range w.r.t energy, g/cm^2/MeV
        LET - stopping power, linear energy transfer, derivative of energy w.r.t range, MeV/(g/cm^2)
        density - return density of material
        Rlimits - return range limits as [min,max] in g/cm^2
        Rlimit_check - check an energy is in table's range
        lookupRange - look up range list
        Elimits - valid energy limits as [min,max], MeV
        Elimit_check - check an energy is in table's range
        lookupEnergy- look up energy list
        get_segments - get monotonic segments of LET list
        LET_spectrum - compute LET spectrum from one or more ion energy spectra
        get_Mcache - compute/retrieve matrix that converts energy spectrum to LET spectrum
    Public Data:
        materials - list of materials
        Species - list of species
        species - selected spicies (None for multi-species tables)
        Z2sym - dict converting from integer Z to symbol

    """
    def __init__(self,materials,Species,species,Z2sym):
        self.materials = materials
        self.Species = Species
        self.species = species
        self.Z2sym = Z2sym
        # these dicts are private because they may not be used by all descendents
        self._Mcache = {} # cached conversions from energy spectrum to LET spectrum. keys are (material,species)
        self._dRdE = {} # cached gradients, keys are material,species
        pass # base class. this must be overloaded in derived class
    def material_check(self,material):
        """
        material_check(material) - check material is known, raise Exception otherwise
        if material is None and table only has one material, returns that
        returns material
        """
        if (material is None) and (len(self.materials)==1):
            return self.materials[0]
        if material not in self.materials:
            raise Exception('Unrecognized material "%s"' % material)
        return material
    def species_check(self,species):
        """
        species_check(species) - check species is known, raise Exception otherwise
        if species is None and table only has one species, returns that
        returns species
        """
        if species is None:
            species = self.species
        if (species is None) and (len(self.Species)==1):
            return self.Species[0]
        species = standardize_species(species)
        if species not in self.Species:
            raise Exception('Unrecognized species "%s"' % species)
        return species
    def density(self,material):
        """
        density(material) - return material density, g/cm^3
        """
        pass # base class. this must be overloaded in derived class
    def Elimits(self,material,species=None):
        """
        Elimits(material,species=None)
        return [min,max] in MeV allowed for material,species
        if species is omitted, all species are used to find bounding energies
        """
        if species is None:
            species = self.Species
        species = np.atleast_1d(species)
        if len(species)==1:
            E = self.lookupEnergy(material,species[0])
            Elimits = np.array([E[0],E[-1]])
        else:
            Elimits = [np.inf,-np.inf]
            for s in species:
                E = self.lookupEnergy(material,species=s)
                Elimits = [np.minimum(Elimits[0],E[0]),np.maximum(Elimits[-1],E[-1])]
        return Elimits
    def Elimit_check(self,material,E,species=None):
        """
        Elimit_check(material,E,species=None) check E MeV against Elimits for material. Raise exception if out of bounds
        """
        Elimits = self.Elimits(material,species) # this will validate material and species
        if np.any(E<Elimits[0]):
            raise Exception('Requested energy %g MeV out of bounds %g-%g' % (np.min(E),Elimits[0],Elimits[1]))
        if np.any(E>Elimits[1]):
            raise Exception('Requested energy %g MeV out of bounds %g-%g' % (np.max(E),Elimits[0],Elimits[1]))
        return True
    def Rlimits(self,material,species=None):
        """
        Rlimits(material,species)
        return [min,max] in g/cm^2 allowed for material, species
        if species is omitted, all species are used to find bounding ranges
        """
        if species is None:
            species = self.Species
        species = np.atleast_1d(species)
        if len(species)==1:
            R = self.lookupRange(material,species[0])
            Rlimits = np.array([R[0],R[-1]])
        else:
            Rlimits = [np.inf,-np.inf]
            for s in species:
                R = self.lookupRange(material,species=s)
                Rlimits = [np.minimum(Rlimits[0],R[0]),np.maximum(Rlimits[-1],R[-1])]
            
        return Rlimits
    def Rlimit_check(self,material,R,species=None):
        """
        Rlimit_check(material,range,species=None) check range g/cm^2 against Rlimits for material. Raise exception if out of bounds
        """
        Rlimits = self.Rlimits(material,species) # this will validate material and species
        if np.any(R<Rlimits[0]):
            raise Exception('Requested range %g g/cm^2 out of bounds %g-%g' % (np.min(R),Rlimits[0],Rlimits[1]))
        if np.any(R>Rlimits[1]):
            raise Exception('Requested range %g g/cm^2 out of bounds %g-%g' % (np.max(R),Rlimits[0],Rlimits[1]))
        return True
    def lookupRange(self,material,species=None):
        """
        R = lookupRange(material,species=None) - return range array in g/cm^2, given material and species
        """
        pass # base class. this must be overloaded in derived class
    def lookupEnergy(self,material,species=None):
        """
        E = lookupEnergy(material,species=None) - return energy array in MeV, given material and species
        """
        pass # base class. this must be overloaded in derived class
        
    def EtoR(self,material,E,species=None):
        """
        range = EtoR(material,E,species=None) - return range in g/cm^2 given material and energy in MeV
        """
        self.Elimit_check(material,E,species)
        y = self.lookupRange(material,species)
        x = self.lookupEnergy(material,species)
        
        f = interp1d(np.log(x),np.log(y),'linear',copy=False,assume_sorted=True)
        return np.exp(f(np.log(E)))
    def RtoE(self,material,R,species=None):
        """
        E = RtoE(material,range,species=None) - return energy in MeV, given material and range in g/cm^2 
        """
        self.Rlimit_check(material,R,species)
        x = self.lookupRange(material,species)
        y = self.lookupEnergy(material,species)
        f = interp1d(np.log(x),np.log(y),'linear',copy=False,assume_sorted=True)
        return np.exp(f(np.log(R)))
    def dRdE(self,material,E,species=None):
        """
        grad = dRdE(material,E,species=None) - returns dRdE in g/cm^2/MeV given material and energy
        (derivative of range w.r.t energy)
        """
        if species is None:
            species = self.species
        species = standardize_species(species)
        _dRdE = None
        x = self.lookupEnergy(material,species)
        if (material,species) in self._dRdE: # retrieve precomputed numerical gradient
            _dRdE = self._dRdE[material,species]
        else: # compute numerical gradient
            y = self.lookupRange(material,species)
            _dRdE = np.gradient(np.log(y),np.log(x))*y/x # derivative with power-law interpolation
            assert np.all(_dRdE>0),'log-log derivative did not work for material %s' % material
            self._dRdE[material,species] = _dRdE
         # linearly interpolate numerical gradient
        f = interp1d(np.log(x),np.log(_dRdE),'linear',copy=False,assume_sorted=True,fill_value='extrapolate')
        return np.exp(f(np.log(E)))
    def LET(self,material,E,species=None):
        """
        grad = LET(material,E,species=None) - returns dE/dR in MeV/(g/cm^2 given material and energy
        (derivative of energy w.r.t range, also LET = linear energy transfer, also stopping power)
        """
        return 1/self.dRdE(material,E,species)
    def get_segments(self,LET):
        """
        segments = get_segments(LET)
        return list of monotonic segments of LET array
        segments = [[ileft,iright],...]
        ileft is index of first element in segment
        iright is index of last element in segment
        segment is ileft:(iright+1) in python slice notation    
        """
        LET = np.array(LET) 
        NE = len(LET)
        segments = [] # index lists for monotonic segments
        ileft = 0
        iright = 1
        s = np.sign(LET[iright]-LET[ileft])
        while iright < NE:
            news = np.sign(LET[iright]-LET[iright-1])
            if (news != s): # iright is in new segment
                segments.append([ileft,iright-1])
                ileft = iright-1
            s = news
            iright = iright+1
        segments.append([ileft,iright-1]) # final segment
        return segments    
    
    def LET_spectrum(self,flux,E,LETgrid,Z=None,species=None,material=None,LET=None,fast=False,excludeH=True):
        """
        fluxLET = LET_spectrum(flux,E,LETgrid,Z=None,species=Non,material=None,LET=None,fast=False,excludeH=True)
        convert flux(E) to flux(LET)
        flux - two options:
            (NE,) flux (differential in E) vs energy array
            (NE,NZ) flux for each species
        E - two options:
            (NE,) energy grid array for all species, MeV
            (NE,NZ) energy for each species
        LETgrid - (NL,) output LET grid array
        LET - two options:
            (NE,) array of LET at E for single species
            (NE,NZ) array of LET at E for each species
        if LET not provided, speficy either species or Z:
            species - str, single species call. flux is (NE,) or (NE,1)
            Z - (NZ,) array of Zs. flux is (NE,NZ)
            material - str, target material (None OK for single-material tables)
        fast - use/store precomputed LET to flux matrix
        excludeH - exclude protons? (protons very rarely should be included in
                   LET spectrum. Usually protons do their SEE via knock-on nuclei)
        fluxLET - (NL,) flux (differential in LET) vs LET array
        if LET is specified, units of E, LET are not important, so long as 
            flux is differential in same energy units as E
            LET and LETgrid use the same units
        flux[i] - flux at energy E[i]
        LET[i] - LET at energy E[i]
        fluxLET[j] - flux at LETgrid[j]
        
        """
        
        E = np.atleast_1d(E)
        NE = E.shape[0]
        flux = np.atleast_1d(flux)
        assert np.array_equal(flux.shape,E.shape), 'flux and E inputs must have same shape'
        
        if flux.ndim==1:
            flux = flux.reshape((NE,1))
            E = E.reshape((NE,1))
        NZ = flux.shape[1]
        # now flux is NE,NZ

        if fast or (LET is None): # will need material and species
            material = self.material_check(material) # replaces None if possible
            assert (species is None) or (Z is None), 'Cannot specify both species and Z. At most one'
            if species is not None:
                species = [self.species_check(sp) for sp in np.atleast_1d(species)] # standardize
            if Z is not None:
                species = [self.Z2sym[z] for z in np.atleast_1d(Z)]
            assert len(species) == NZ, 'flux shape does not match number of species/Z'
        # now use species when iterating
        Z = None # force errors to be thrown if Z is accessed
        
        if excludeH and ('H' in species):
            #exclude protons
            iZ = [s != 'H' for s in species]
            flux = flux[:,iZ]
            if E.shape[1] == len(species):
                E = E[:,iZ]
            NZ = flux.shape[1]
            species = [s for s in species if s != 'H']

        LETgrid = np.atleast_1d(LETgrid).ravel() # (NL,)
        fluxLET = np.zeros(LETgrid.shape) # (NL,)
        
        if E.shape[1] < NZ:
            dE = 0 # do not advance energy index
        else:
            dE = 1  # advance energy index
        
        if fast: # use matrix/cache
            for [iz,sp] in enumerate(species):
                if LET is None:
                    M = self.get_Mcache(E[:,iz*dE],LETgrid,species=sp,material=material)
                else:
                    M = self.get_Mcache(E[:,iz*dE],LETgrid,species=sp,material=material,LET=LET[:,iz])
                fluxLET += np.dot(M,flux[:,iz])
            return fluxLET
                
        
        if LET is None: # build LET for each species
            LET = np.full((NE,NZ),np.nan)
            for [iz,sp] in enumerate(species):
                LET[:,iz] = self.LET(material,E[:,iz*dE],species=sp)
        else:
            LET = np.atleast_1d(LET)
            if LET.ndim == 1:
                LET = LET.reshape((len(LET),1))
            assert np.array_equal(LET.shape,flux.shape), 'flux and LET inputs must have same shape'
        
        # accumulate fluxLET
        for [iz,sp] in enumerate(species):
            segments = self.get_segments(LET[:,iz])
            
            for seg in segments:
                Es = E[seg[0]:(seg[1]+1),iz*dE]
                Ls = LET[seg[0]:(seg[1]+1),iz]
                js = flux[seg[0]:(seg[1]+1),iz]
                s = np.sign(Ls[-1]-Ls[0])
                if s==0:
                    continue # flat LET vs E segments do not contribute to fluxLET
                
                # multiply abscissa by s so that s*x or s*log(Ls) is monotonically increasing
                if np.all(js>0): # log-log interp
                    jfunc = lambda Li: np.exp(interp1d(s*np.log(Ls),np.log(js),kind='linear',copy=False,assume_sorted=True,bounds_error=False,fill_value=-np.inf)(s*np.log(Li)))
                else: # semilogx interp
                    jfunc = lambda Li: interp1d(s*np.log(Ls),js,kind='linear',copy=False,assume_sorted=True,bounds_error=False,fill_value=0)(s*np.log(Li))
        
                dEdLET = np.abs(np.gradient(np.log(Es),np.log(Ls)))*Es/Ls # |dE/dLET| in a power-law sense
                # log-log interpolating function for dE/dLET
                gfunc = lambda Li: np.exp(interp1d(s*np.log(Ls),np.log(dEdLET),kind='linear',copy=False,assume_sorted=True,bounds_error=False,fill_value=0)(s*np.log(Li)))
        
                # j(LET) = j(E(LET))*dE/dLET
                fluxLET += jfunc(LETgrid)*gfunc(LETgrid) # both funcs zero out of bounds 
        
        return fluxLET
        
    
    def get_Mcache(self,E,LETgrid,LET=None,species=None,material=None):
        """
        M = LET_spectrum_matrix(E,LETgrid,LET=None,species=None,material=None)
        E - (NE,) energy grid array
        LET - (NE,) LET vs energy array
        LETgrid - (NL,) output LET grid array
        M - (NL,NE) matrix such that flux(LET) = M*flux(energy)
        units of E, LET are not important, so long as 
            flux is differential in same energy units as E
            LET and LETgrid use the same units
        LET[i] - LET at energy E[i]
        flux[i] - flux (differential in E) at energy E[i]
        (M*flux)[j] - flux (differential in LET) at LETgrid[j], * is dot product
        species - str, species of interest
        material - str, material of interest
        """
        
        species = self.species_check(species)
        material = self.material_check(material)
        
        E = np.atleast_1d(E).ravel() # (NE,)
        LETgrid = np.atleast_1d(LETgrid).ravel() # (NL,)
        
        if (((material,species) in self._Mcache) and
            (np.array_equal(E,self._Mcache[material,species]['E'])) and
            (np.array_equal(LETgrid,self._Mcache[material,species]['LETgrid']))):
            return self._Mcache[material,species]['M'] # cache is still good
        
        if LET is None:
            LET = self.LET(material,E,species=species)
        else:
            LET = np.atleast_1d(LET).ravel() # (NE,)
            assert len(E) == len(LET), 'E and LET inputs must have same length'

        segments = self.get_segments(LET) # get monotonic segments of LET vs E
        
        M = np.zeros((len(LETgrid),len(LET)))
        for seg in segments:
            Es = E[seg[0]:(seg[1]+1)]
            Ls = LET[seg[0]:(seg[1]+1)]
            s = np.sign(Ls[-1]-Ls[0])
            if s==0:
                continue # flat LET vs E segments do not contribute
    
            dEdLET = np.abs(np.gradient(np.log(Es),np.log(Ls)))*Es/Ls # |dE/dLET| in a power-law sense
            # log-log interpolating function for dE/dLET
            # multiply abscissa by s so that s*x or s*log(Ls) is monotonically increasing
            gfunc = lambda Li: np.exp(interp1d(s*np.log(Ls),np.log(dEdLET),kind='linear',copy=False,assume_sorted=True,bounds_error=False,fill_value=0)(s*np.log(Li)))
    
            # j(LET) = j(E(LET))*dE/dLET
            M[:,seg[0]:(seg[1]+1)] += np.dot(np.diag(gfunc(LETgrid)),interpolating_matrix(np.log(Ls),np.log(LETgrid),out_of_bounds='ZERO'))
        self._Mcache[material,species] = {'E':E,'LETgrid':LETgrid,'M':M}
        return M

    
    
class ScreamTable(RangeTable):
    """ table = ScreamTable(filename=None)
    multi-material, single-species table of Scream or NIST ranges
    filename = is the relative or absolute path to the range table file
    default filename is ScreamRangeData_protons.csv (proton range table)
    for electrons, supply ScreamRangeData_electrons.csv
    NIST proton range can be used via NISTRangeData_protons.csv
    Public Methods:
        same as RangeTable
    Public Data:
        same as RangeTable plus:
        filename - as above
        ranges - dict containing range-energy data
    """
    def __init__(self,filename=None):
        
        if filename is None:
            filename = 'ScreamRangeData_protons.csv'
        filename = locate_file(filename)
        if not os.path.exists(filename):
            raise Exception('%s does not exist' % filename)
        if 'electron' in filename:
            Species = ['e']
            species = 'e'
            Z2sym = {-1:'e'}
        else: # all others assume protons
            if 'proton' not in filename:
                warnings.warn('guessing protons for %s' % filename)
            Species = ['H']
            species = 'H'
            Z2sym = {1:'H'}
    
        ranges = {}
        """ranges is a dict
            ranges['Energy (MeV)'] is a list of energies in MeV
            all other keys are materials with values that are themselves dicts
            value dicts have two fields: density and range
            density is a scalar in g/cm^3
            range is alist of g/cm^2 corresponding to entries in Energy list"""
        
            
        """density,2.32,2.32,...
           Energy,,,...
           (MeV),Si,Si2,...
           1.000E-04,6.264E-07,"""
    
        with open(filename,'r') as f:
            headers = []
            for i in range(3):
                headers.append(f.readline().strip().split(','))
            densities = np.array([float(x) for x in headers[0][1:]])
            materials = [(a+' '+b).strip() for (a,b) in zip(headers[1],headers[2])]
            data = []
            for line in f.readlines():
                line = line.strip();
                if line == ','*len(line): # all commas
                    continue
                line = [float(x) for x in line.split(',')]
                data.append(line)
            data = np.array(data)
        for i in range(len(materials)):
            if i == 0:
                ranges[materials[0]] = data[:,0] # energy
            else:
                rec = {'density':densities[i-1]}
                rec['range'] = data[:,i]
                ranges[materials[i]] = rec
        ranges['materials'] = materials[1:] # remove first one, which is energy
        
        super().__init__(ranges['materials'],Species,species,Z2sym)
        self.ranges = ranges # shorthand
        self.filename = filename
        
    def __str__(self):
        return "Scream %s Range Table" % self.species
    def density(self,material):
        """
        density(material) - return material density, g/cm^3
        """
        self.material_check(material)
        return self.ranges[material]['density']

    def lookupRange(self,material,species=None):
        """
        R = lookupRange(material,species=None) - return range array in g/cm^2, given material and species
        """
        self.material_check(material)
        self.species_check(species)
        R = self.ranges[material]['range']
        return R
    def lookupEnergy(self,material,species=None):
        """
        E = lookupEnergy(material,species=None) - return energy array in MeV, given material and species
        """
        self.material_check(material)
        self.species_check(species)
        E = self.ranges['Energy (MeV)']
        return E

class SRIMTable(RangeTable):
    """ table = SRIMTable(filename)
    multi-species, single-material table of SRIM ranges
    filename = is the relative or absolute path to the range table file (SRIM_Al.csv)
    file provides target material and its density
    filename is expected to be comma separated values with three headers:
        One of those headers is "# Density = X g/cm^3" where X is a numeric density
        One of those headers is "SRIM range/LET in X ..." where X is the material
        After that, there are 4 lines per element:
            line 1 gives symbol,name,Z,A
            line 2 gives the Energies in MeV (has a 1-column label)
            line 3 gives the Range in um (has a 1-column label)
            line 4 gives the LET in (MeV/(mg/cm^2)) (has a 1-column label)
    Public Methods:
        same as RangeTable
    Public Data:
        same as RangeTable plus:
            filename - as above
            data which is a dict with fields:
            H,He,... - entry for each element is a dict with fields:
                symbol - element symbol
                name - element common name
                Z - atomic number
                A - atomic mass (appears to be abundance-weighted mass)
                MAI - most abundant isotope (integer mass number)
                MeV - energy list in MeV
                range - range in g/cm^2
                LET - LET in MeV/(g/cm^2)
    """
    def __init__(self,filename=None):
        
        if filename is None:
            filename = 'SRIM_Al.csv'
        filename = locate_file(filename)
        if not os.path.exists(filename):
            raise Exception('%s does not exist' % filename)
        density = None
        data = {}
        with open(filename,'r') as f:
            for i in range(3): # read 3 header lines
                line = f.readline()
                m = re.match('.*SRIM\s+range/LET\s+in\s+(\S+)',line)
                if m:
                    material = m.group(1)
                m = re.match('.*Density\s*=\s*(\S+)',line)
                if m:
                    density = float(m.group(1))
            if material is None: 
                raise Exception('Material not specified in %s' % filename)
            if density is None: 
                raise Exception('Density not specified in %s' % filename)
            iline = 0;
            for line in f:
                line = line.strip()
                iline = (iline % 4) + 1
                if iline == 1:
                    (symbol,name,Z,A,MAI) = re.split('\s*,\s*',line)
                    data[symbol] = {'symbol':symbol,'name':name,'Z':int(Z),'A':float(A),'MAI':int(MAI)}
                elif iline == 2:
                    values = re.split('\s*,\s*',line)
                    data[symbol]['MeV'] = np.array([float(v) for v in values[1:]])
                elif iline == 3:
                    values = re.split('\s*,\s*',line)
                    data[symbol]['range'] = np.array([float(v) for v in values[1:]]) # um
                    data[symbol]['range'] *= 1e-4*density # um -> g/cm^2
                    assert len(data[symbol]['range']) == len(data[symbol]['MeV']),'Element %s in %s has inconsistent MeV and range lists' % (symbol,filename)
                elif iline == 4:
                    values = re.split('\s*,\s*',line)
                    data[symbol]['LET'] = np.array([float(v) for v in values[1:]]) # MeV/(mg/cm^2)
                    data[symbol]['LET'] *= 1e3 # MeV/(mg/cm^2) -> MeV/(g/cm^2)
                    assert len(data[symbol]['LET']) == len(data[symbol]['MeV']),'Element %s in %s has inconsistent MeV and LET lists' % (symbol,filename)

        Z2sym = {d['Z']:d['symbol'] for d in data.values()}
        Species = list(data.keys())
        species = None # All calls should specify species
        super().__init__([material],Species,species,Z2sym) 
        self.filename = filename
        self.data = data
        self._density = density

    def __str__(self):
        return "Scream Range Table for material %s" % self.materials[0]
    def dRdE(self,material,E,species=None):
        """
        grad = dRdE(material,E,species=None) - returns dRdE in g/cm^2/MeV given material and energy
        (derivative of range w.r.t energy)
        """
        return 1/self.LET(material,E,species) # SRIM table has LET, use that
    def LET(self,material,E,species=None):
        """
        grad = LET(material,E,species=None) - returns dE/dR in MeV/(g/cm^2 given material and energy
        (derivative of energy w.r.t range, also LET = linear energy transfer, also stopping power)
        """
        species = standardize_species(species)
        x = self.lookupEnergy(material,species)
        _LET = self.data[species]['LET']
         # linearly interpolate numerical gradient
        f = interp1d(np.log(x),np.log(_LET),'linear',copy=False,assume_sorted=True,fill_value='extrapolate')
        return np.exp(f(np.log(E)))
    def density(self,material):
        """
        density(material) - return material density, g/cm^3
        """
        self.material_check(material)
        return self._density # only one material

    def lookupRange(self,material,species=None):
        """
        R = lookupRange(material,species=None) - return range array in g/cm^2, given material and species
        """
        self.material_check(material)
        self.species_check(species)
        species = standardize_species(species)
        R = self.data[species]['range']
        return R
    def lookupEnergy(self,material,species=None):
        """
        E = lookupEnergy(material,species=None) - return energy array in MeV, given material and species
        """
        self.material_check(material)
        self.species_check(species)
        E = self.data[species]['MeV']
        return E
        

class ShieldLayer(object):
    """
    ShieldLayer(table,material,thickness,units='mils',geometry='UNKNOWN',one_sided_correction=1)
    Shield layer abstract base class. Manages shield layer properties and methods
    table - object with base class RangeTable
    material - string naming material supported by RangeTable
    thickness - numeric thickness
    units - string naming supported units (mils, in, cm, mm, um, A, g/cm^2, mg/cm^2)
    overload the K function to apply to specific geometry
    one_sided_correction - set to 1/2 for one-sided layers, like slabs or hemisphere
    geometry: 'Sphere','Slab', etc
    Public Methods:
        tocm - convert shielding depth to centimeters
        Eout - output energy given incident energy
        Ein - incident energy given output energy
        dRdE - derivative of range with energy
        K - path length distribution (abstract:overload)
        makeQ - make range spectrum conversion matrix
        get_cache - get the cached D matrix that degrades energy spectrum
        slow - compute degraded spectrum via adaptive quadrature
        degraded_spectrum - compute degraded spectrum
        
    Public Data:
        table - RangeTable
        material - str, shielding material
        thickness - float,shielding thickness
        units - str, units of shielding thickness
        density - float, density of shielding material
        cm - float, thickness of shielding in cm
        gcm2 - float, thickness of shielding in g/cm^2
        Tmin - float, minimum shielding thickness along any ray
        one_sided_correction - bool, whether shield is considered one-sided (2pi incidence)
        geometry - str, short geometry identifier (e.g., 'Slab')
    """
    def __init__(self,table,material,thickness,units='mils',geometry='UNKNOWN',one_sided_correction=1):
        self.table = table
        self.material = material
        self.thickness = thickness
        self.units = units
        self.density = table.density(material)
        self.cm = self.tocm()
        self.gcm2 = self.cm*self.density
        self.Tmin = self.gcm2 # minimum thickness of shielding, K(T<Tmin) = 0
        self.one_sided_correction = one_sided_correction
        self.geometry = geometry
        self._cache = None # degraded_spectrum cache
    def __repr__(self):
        return self.__class__.__name__ + ' ' + self.__str__() # ShieldLayer 5 mils SiO2
    def __str__(self): # e.g., 5 mils SiO2
        return '%s %s %s %s' % (self.geometry,self.thickness,self.units,self.material)
    def tocm(self,thickness=None,units=None,density=None):
        """tocm(thickness=None,units=None,density=None)
        thickness numeric shield thickness
        units is string naming any of the units supported by module function tocm
        returns thickness in cm
        """
        if thickness is None:
            thickness = self.thickness
        if units is None:
            units = self.units
        if density is None:
            density = self.density
        return tocm(thickness,units,density)
    def Eout(self,Ein,angle=0,species=None,T=None):
        """
        Eout = Eout(Ein,angle=0,species=None,T=None)
        Compute exit energy for particle incident with energy Ein in MeV
        Returns exit energy in MeV
        angle in radians from normal
        T = shielding thickness in g/cm^2 (used instead of self.gcm2)
        species - species of interest (str) defaults to range table's species
        Returns 0 if absorbed
        angle and species are scalars
        """
        if T is None:
            T = self.gcm2
        if not (np.isscalar(Ein) and np.isscalar(angle) and np.isscalar(T)):
            return np.vectorize(self.Eout)(Ein,angle,species,T)
        
        dR = T/np.cos(angle) # secant adjusted path length through layer
        Rin = self.table.EtoR(self.material,Ein,species=species)
        Rout = Rin-dR
        if Rout<=self.table.Rlimits(self.material,species=species)[0]: # below Rlimit, treat as absorbed
            return 0
        else:
            return self.table.RtoE(self.material,Rout,species=species)
    def Ein(self,Eout,angle=0,species=None,T=None):
        """
        Ein = Ein(Eout,angle=0,species=None,T=None)
        Compute incident energy for particle exiting with energy Eout in MeV
        Returns incident energy in MeV
        angle in radians from normal
        T = shielding thickness in g/cm^2 (used instead of self.gcm2)
        species - species of interest (str) defaults to range table's species
        returns inf if above range table Elimits
        angle and species are scalars
        """
        if T is None:
            T = self.gcm2
        if not (np.isscalar(Eout) and np.isscalar(angle) and np.isscalar(T)):
            return np.vectorize(self.Ein)(Eout,angle,species,T)
        dR = T/np.cos(angle) # secant adjusted path length through layer
        Rout = self.table.EtoR(self.material,Eout,species=species)
        Rin = Rout+dR
        if Rin > self.table.Rlimits(self.material,species=species)[1]:
            return np.inf
        else:
            return self.table.RtoE(self.material,Rin,species=species)
        
    def dRdE(self,E,species=None):
        """
        grad = dRdE(E,species=None)
        derivative of range with respect to energy at E
        E in MeV
        grad in g/cm^2/MeV
        species: species of interest (str) defaults to range table's species
        """
        return self.table.dRdE(self.material,E,species=species)
    def K(self,T):
        """
        K(T) - probability of thickness T for rays passing through shield
        abstract method - must be overloaded
        """
        raise NotImplementedError
    def makeQ(self,dEoverE=0.1,epsrel=1e-4,species=None,Tfast=True):
        """
        makeQ(dEoverE=0.1,epsrel=1e-4,species=None,Tfast=True)
        add T, Q to cache as _T_ and _Q_
        T - list of working ranges
        Q - matrix that converts from incident range-differential flux 
            to exiting rangle-differential flux
        see get_cache for dEoverE,epsrel,Tfast definitions
        """
        if VERBOSE:
            print('New %s Q cache for %s' % (species,self))
        
        if self._cache is None:
            self._cache = {}
        Rlimits = self.table.Rlimits(self.material,species=None) # Rlimits over all species
        # number of points in range integral 100 or enough to have dEoverE spacing
        NT = np.maximum(100,int(np.ceil(25*0.1/dEoverE*np.log10(Rlimits[1]/Rlimits[0]))))
        T = np.exp(np.linspace(np.log(Rlimits[0]),np.log(Rlimits[1]),NT)) # working range grid
        self._cache['_T_'] = T
        Q = np.zeros((NT,NT))
        if Tfast:
            dT = trapz_dx(T) # integration weights
            for i in range(NT):
                # A = interpolating_matrix(x,y). u(y) = A(x,y)*u(x)
                # interpolates from x to y
                # here we interpolate from T[i]+T to T (incident range onto range grid)
                Q[i,:] = np.dot(interpolating_matrix(np.log(T),np.log(T+T[i])),self.K(T)*dT)
        else: # not-so-fast, fast==True or fast == 1
            # do integral version
            for (m,Tm) in enumerate(T): # exiting range index m
                for (mpp,Tmpp) in enumerate(T): # entering range index m''
                    for side in ['left','right']: # integrate left, right side of trapezoidal rule
                        if side == 'left':
                            if mpp==0:
                                continue # no neighbor on left side
                            m0 = mpp-1 # index where weight is 0
                        else: # right
                            if mpp == NT-1:
                                continue # no neighbor on right side
                            m0 = mpp+1
                        # weight is zero at Tm
                        # weight is one at Tmpp
                        R0 = T[m0] # Range where weight is 0
                        def func(R):
                            """ R is particle range """
                            Tp = R-Tm # thickness of material along current ray, T'
                            weight = (R-R0)/(Tmpp-R0) # trapezoidal weight: 0 at R0, 1 at Tmpp
                            return weight*self.K(Tp)
                    
                        R1 = np.minimum(Tmpp,R0)
                        R2 = np.maximum(Tmpp,R0)
                        if R2>R1:
                            (y,abserr) = quad(func,R1,R2,epsrel=epsrel)
                            Q[m,mpp] += y# integral A*K*dT'
        self._cache['_Q_'] = Q
        
    def get_cache(self,dEoverE=0.1,epsrel=1e-4,species=None,Tfast=True):
        """cache = get_cache(dEoverE=0.1,epsrel=1e-4,species=None,Tfast=True)
        return the cache that stores precomputed transform
        cache[species] - D degrading matrix, (NE,NE) on the species' range-energy grid
        cache['_T_'] - working range grid (NT,)
        cache['_Q_'] - species-independent transform on working grid (NT,NT)
        exitflux = dot(transform,influx) where fluxes are on the range table's energy grid
        options:
            dEoverE: relative spacing of working grid (dT/T = dE/E)
            epsrel: epsrel argument passed to quad numerical integration (relative precision)
            species: species of interest (str) defaults to range table's species
            Tfast : perform KdT integral using quadrature (False) or trapezoidal (True) method            
        """
        # build a transform that, for each exiting energy and incident energy, integrates over just the relevant incident angles
        if species is None:
            species = self.table.species
        species = standardize_species(species)
        self.table.species_check(species)
        
        if self._cache is None:
            self._cache = {}
        if species in self._cache:
            cache = self._cache[species]
        else:
            cache = None            
        if cache is None: # prepare cache
            if VERBOSE:
                print('New %s cache for %s' % (species,self))
            if '_Q_' not in self._cache:
                # make species independent Q
                self.makeQ(dEoverE=dEoverE,epsrel=epsrel,Tfast=Tfast)
                
            E = self.table.lookupEnergy(self.material,species)
            R = self.table.lookupRange(self.material,species)
            T = self._cache['_T_']
            Q = self._cache['_Q_']
            # interpolate from working grid T to rangel tables grid R
            A = interpolating_matrix(np.log(T),np.log(R))
            At = interpolating_matrix(np.log(R),np.log(T)) # At ~ transpose(A)
            dRdE = self.dRdE(E,species=species)
            cache = np.dot(np.diag(dRdE),np.dot(A,np.dot(Q,np.dot(At,np.diag(1/dRdE))))) # dRdE*A*Q*A'/dRdE
            self._cache[species] = cache
        return cache

    def slow(self,energy,flux,epsrel=1e-4,species=None):
        """
        exitflux = slow(energy,flux,epsrel=1e-4,species=None)
        compute degraded spectrum via adaptive quadrature
        energy - array of energies
        flux - incidient flux at energy grid points
        epsrel - relative precision of numerical integration, passed to quad
        species - species of interest (str) defaults to range table's species
        exitflux - degraded flux at energy grid points
        """
        exitflux = np.full(energy.shape,np.nan)        
        Elimits = self.table.Elimits(self.material,species=species)
        influx = lambda E: smart_interp1d(np.log(energy),flux,np.log(E),'linear',copy=False,assume_sorted=True,bounds_error=False,fill_value=0)
        
        def func(lnT,Eout,dRdEout):
            T = np.exp(lnT)
            Ein = self.Ein(Eout,species=species,T=T)
            if (Ein < energy[0]) or (Ein > energy[-1]):
                return 0
             # extra factor because integrating dlnT
            return influx(Ein) * dRdEout / self.dRdE(Ein,species=species) * self.K(T)*T
        
        Rlimits = self.table.Rlimits(self.material,species=species)
        Tmin = np.maximum(self.Tmin,Rlimits[0]) # lower bound of integration
        for (iE,Eout) in enumerate(energy):
            if Eout <= Elimits[-1]:
                dRdEout = self.dRdE(Eout,species=species)
                # integral dlnT
                (y,abserr)= quad(func,np.log(Tmin),np.log(Rlimits[-1]),args=(Eout,dRdEout),epsrel=epsrel)
                exitflux[iE] = y
        self._correct_high_energy(flux,energy,exitflux,Elimits[-1])
        return exitflux
        
    def _correct_high_energy(self,influx,energy,exitflux,Elimit=None):
        """
            _correct_high_energy(exitflux,influx,energy,Elimit=None)
            correct flux at high energies assuming no effect from shield
            influx - entering flux (array)
            energy - energy array
            exitflux - exiting flux (array)
            Elimit - apply correction where energy>Elimit and next grid below
        """
        global _HIGH_ENERGY_WARNED
        iHigh = energy>Elimit
        if np.any(iHigh):
            if not _HIGH_ENERGY_WARNED:
                warnings.warn('Warning: %g MeV requested. Assuming negligible energy loss for all energies > %g MeV' 
                              % (np.min(energy[iHigh]),Elimit))
                _HIGH_ENERGY_WARNED = True
            exitflux[iHigh] = influx[iHigh]*self.one_sided_correction # assume highest energies penetrate unaffected

    def degraded_spectrum(self,energy,flux,fast=True,epsrel=1e-4,dEoverE=0.1,species=None,Tfast=False):
        """
        exitflux = degraded_spectrum(energy,flux,fast=True,epsrel=1e-4,dEoverE=0.1,species=None,Tfast=False)
        compute degraded spectrum
        energy - array of incident energies
        flux - array of incident fluxes (integrated over 4pi sr, differential #/cm^2/MeV)
        exitflux - array of transmitted fluxes (integrated over emerging 2pi sr, differential #/cm^2/MeV)
        options:
            fast - True/1/2/False
                if False performs quad adaptive integral on every call (slowest)
                if True precomputes transform matrix and uses that on all subsequent calls:
                if 1 also forces Tfast = False 
                if 2 also forces Tfast = True (fastest)
            Tfast - False/True
                if False uses quadrature integral in T' for pre-computed transform
                if True uses trapezoidal integral in T' and forces fast=2 (fastest)
            epsrel - relative precision of numerical integration, passed to quad
            dEoverE - dEoverE passed to self.get_cache, realtive spacing of working grid
            species - species of interest (str) defaults to range table's species
        """
        
        # handle fast, Tfast forcing
        if fast == 1:
            Tfast = False
        elif fast == 2:
            Tfast = True

        if Tfast and not fast:
            fast = True
        if fast: #use precomputed cache
            influx = lambda E: smart_interp1d(np.log(energy),flux,np.log(E),'linear',copy=False,assume_sorted=True,bounds_error=False,fill_value=0)
            Egrid = self.table.lookupEnergy(self.material,species)
            cache = self.get_cache(dEoverE=dEoverE,epsrel=epsrel,species=species,Tfast=Tfast)
            tmpflux = np.dot(cache,influx(Egrid))
            exitflux = smart_interp1d(np.log(Egrid),tmpflux,np.log(energy),'linear',copy=False,assume_sorted=True,bounds_error=False,fill_value=0)
            self._correct_high_energy(flux,energy,exitflux,Egrid[-2]) # include 2nd to last energy grid point
        else: # do numerical integral (quad)
            exitflux = self.slow(energy,flux,epsrel=epsrel,species=species)
        return exitflux

class SlabLayer(ShieldLayer):
    """
    SlabLayer(table,material,thickness,units='mils',one_sided=True)
    manages infinite slab layer properties and methods
    table - object with base class RangeTable
    material - string naming material supported by RangeTable
    thickness - numeric thickness
    units - string naming supported units (mils, in, cm, mm, um, A, g/cm^2, mg/cm^2)
    one_sided - represent one-sided slab. Add two slab outputs to get flux between them.
        if not one_sided, output flux is doubled as if target is sandwiched between two slabs
    Public Data:
        same as ShieldLayer plus
        one_sided - as above
    Public Methods:
        same as ShieldLayer
    """
    def __init__(self,table,material,thickness,units='mils',one_sided=True):
        super().__init__(table,material,thickness,units=units,geometry='Slab',one_sided_correction = 1/2 if one_sided else 1)
        self.one_sided = one_sided
    def __str__(self): # e.g., 5 mils SiO2 one-sided
        return super().__str__() + ' ' + ('one-sided' if self.one_sided else 'two-sided')
    def K(self,T):
        """
        K(T) - probability of thickness T for rays passing through shield
        """
        isscalar = np.isscalar(T)
        T = np.atleast_1d(T)
        K = self.one_sided_correction*self.gcm2/T**2
        K[T<self.gcm2] = 0
        if isscalar:
            return K[0]
        else:
            return K

class SphereLayer(ShieldLayer):
    """
    SphereLayer(table,material,thickness,units='mils')
    manages spherical shell layer properties and methods
    table - object with base class RangeTable
    material - string naming material supported by RangeTable
    thickness - numeric thickness
    units - string naming supported units (mils, in, cm, mm, um, g/cm^2)
    Public Data:
        same as ShieldLayer
    Public Methods:
        same as ShieldLayer
    """
    def __init__(self,table,material,thickness,units='mils'):
        super().__init__(table,material,thickness,units=units,geometry='Sphere')
    def get_cache(self,dEoverE=0.1,epsrel=1e-4,species=None,Tfast=True):
        """cache = get_cache(dEoverE=0.1,epsrel=1e-4,species=None,Tfast=True)
        return the cache that stores precomputed transform
        cache[species] - (NE,NE) on the species' range-energy grid
        exitflux = dot(transform,influx) where fluxes are on the range table's energy grid
        options:
            dEoverE: relative spacing of working grid (dT/T = dE/E)
            epsrel: epsrel argument passed to quad numerical integration (relative precision)
            species: species of interest (str) defaults to range table's species
            Tfast : not used
        """
        # build a transform that, for each exiting energy and incident energy, integrates over just the relevant incident angles
        if species is None:
            species = self.table.species
        species = standardize_species(species)
        self.table.species_check(species)
        
        
        if self._cache is None:
            self._cache = {}
        if species in self._cache:
            cache = self._cache[species]
        else:
            if VERBOSE:
                print('New %s cache for %s' % (species,self))
            E = self.table.lookupEnergy(self.material,species)
            Ein = self.Ein(E,species=species)
            i = (Ein>=E[0]) & (Ein<=E[-1])
            # interpolate from working grid T to rangel tables grid R
            dRdEout = self.dRdE(E[i],species=species)
            dRdEin = self.dRdE(Ein[i],species=species)
            A = np.zeros((len(E),len(E)))
            A[i,:] = interpolating_matrix(E,Ein[i]) # interpolate input flux onto Ein
            A[i,:] = np.dot(np.diag(dRdEout/dRdEin),A[i,:])
            cache = A
            self._cache[species] = cache
        return cache

    def slow(self,energy,flux,epsrel=1e-4,species=None):
        """
        exitflux = slow(energy,flux,epsrel=1e-4,species=None)
        compute degraded spectrum via adaptive quadrature
        inputs/outputs same as ShieldLayer.slow,except:
            epsrel is not used
        """
        
        Elimits = self.table.Elimits(self.material,species)
        influx = lambda E: smart_interp1d(np.log(energy),flux,np.log(E),'linear',copy=False,assume_sorted=True,bounds_error=False,fill_value=0)

        j = (energy>=Elimits[0]) & (energy<=Elimits[1]) # in-range input energies
        dRdEout = self.dRdE(energy[j],species=species)
        Ein = self.Ein(energy[j],species=species)
        i = (Ein>=Elimits[0]) & (Ein<=Elimits[1]) # in range incident energies
        j[j] = i # combine energy and Ein constraints
        dRdEin = self.dRdE(Ein[i],species=species)
        exitflux = np.zeros(len(energy))
        exitflux[j] = influx(Ein[i])*dRdEout[i]/dRdEin
        self._correct_high_energy(flux,energy,exitflux,Elimits[1])
        return exitflux

if __name__ == '__main__':
    # perform three demos
    import matplotlib.pyplot as plt
    plt.close('all')
    
    # prepare figure path
    _TOP_PATH = os.path.dirname(os.path.abspath(__file__))
    figpath = os.path.join(_TOP_PATH,'figures')
    if not os.path.exists(figpath):
        os.mkdir(figpath)
    
    
    # load range tables
    stable = ScreamTable()
    ntable = ScreamTable('NISTRangeData_protons.csv')
    etable = ScreamTable('ScreamRangeData_electrons.csv')

    al_table = SRIMTable() # Al by default
    si_table = SRIMTable('SRIM_Si.csv')
    
    # proton fluence
    pMeV = 10**np.linspace(0,3.5,51) # 1 MeV to ~3 GeV
    pFluence = pMeV**-2 # falling spectrum
    
    # electron fluence
    eMeV = 10**np.linspace(-1,0.5,51) # 0.1 MeV to ~3 MeV
    eFluence = eMeV**-3 # falling spectrum


    # HE fluence
    HeMeV = 10**np.linspace(0,3.5,51)*4 # 1 MeV/nuc to ~3 GeV/nuc
    HeFluence = 0.1*HeMeV**-2 # 10 % of protons at same MeV/nuc
    ions = {'H':{'MeV':pMeV,'fluence':pFluence},'He':{'MeV':HeMeV,'fluence':HeFluence}}
    
    # demo LET spetrum
    
    LET = np.exp(np.linspace(np.log(1),np.log(110000),1002))
    Z = [si_table.data[sp]['Z'] for sp in ions.keys()]
    ionFluence = np.stack([ions[sp]['fluence'] for sp in ions.keys()],axis=1)
    ionEnergies = np.stack([ions[sp]['MeV'] for sp in ions.keys()],axis=1)
    LETflux_slow = si_table.LET_spectrum(ionFluence,ionEnergies,LET,Z=Z,fast=False)
    LETflux_fast = si_table.LET_spectrum(ionFluence,ionEnergies,LET,Z=Z,fast=True)
    plt.figure()
    plt.loglog(LET,LETflux_slow,'k-',LET,LETflux_fast,'r.')
    plt.xlabel('LET, MeV/(g/cm$^2$)')
    plt.ylabel('#/cm$^2$/(MeV/(g/cm$^2$))')
    plt.legend(['slow','fast'])
    plt.title('LET Spectrum for %s' % (si_table))
    plt.savefig(os.path.join(figpath,'SRIM_let_slow_fast.png'))

    # demo sphere ------------------------------
    
    # compare ions SRIM slow to fast method
    srim_sphere = SphereLayer(al_table,'Al',100) # 100 mils Al
    plt.figure()
    colors = ['k','darkgreen']
    for (isp,sp) in enumerate(ions):
        plt.loglog(ions[sp]['MeV'],ions[sp]['fluence'],'-',color=colors[isp],label='%s Incident' % sp)
        ions[sp]['slow'] = srim_sphere.degraded_spectrum(ions[sp]['MeV'],ions[sp]['fluence'],fast=False,species=sp)
        plt.loglog(ions[sp]['MeV'],ions[sp]['slow'],'b--',label='%s Slow' % sp)
        ions[sp]['fast'] = srim_sphere.degraded_spectrum(ions[sp]['MeV'],ions[sp]['fluence'],fast=True,species=sp)
        plt.loglog(ions[sp]['MeV'],ions[sp]['fast'],'r.',label='%s Fast' % sp)
    plt.xlabel('MeV')
    plt.ylabel('ions/cm$^2$/MeV')
    plt.title('SRIM Comparison of Fast and Slow methods %s' % repr(srim_sphere),fontsize='smaller')
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(figpath,'SRIM_degraded_slow_fast.png'))
    
    # compare protons SCREAM slow to fast method
    sphere = SphereLayer(stable,'Al',100) # 100 mils Al
    sph_slow = sphere.degraded_spectrum(pMeV,pFluence,fast=False)
    sph_fast = sphere.degraded_spectrum(pMeV,pFluence,fast=True)
    plt.figure()
    plt.loglog(pMeV,pFluence,'k-',pMeV,sph_slow,'b-',pMeV,sph_fast,'r--')
    plt.loglog(ions['H']['MeV'],ions['H']['slow'],'.',color='darkgreen')
    plt.xlabel('MeV')
    plt.ylabel('protons/cm$^2$/MeV')
    plt.title('Comparison of Scream Fast and Slow methods %s' % repr(sphere),fontsize='smaller')
    plt.legend(['Incident','Slow','Fast','SRIM'])
    plt.grid(True)
    plt.savefig(os.path.join(figpath,'SCREAM_SRIM_degraded_slow_fast_sphere.png'))
    
    # demo slabs ------------------------------
    
    # compare protons SCREAM slow to fast method
    front = SlabLayer(stable,'SiO2',15) # 15 mils of SiO2
    slow = front.degraded_spectrum(pMeV,pFluence,fast=False)
    fast = front.degraded_spectrum(pMeV,pFluence,fast=1)
    fast2 = front.degraded_spectrum(pMeV,pFluence,fast=2)
    plt.figure()
    plt.loglog(pMeV,pFluence/2,'k-',pMeV,slow,'b-',pMeV,fast,'r--',pMeV,fast2,'g.')
    plt.xlabel('MeV')
    plt.ylabel('protons/cm$^2$/MeV')
    plt.title('Comparison of Fast and Slow methods %s' % repr(front),fontsize='smaller')
    plt.legend(['Incident/2','Slow','Faster','Fastest'])
    plt.grid(True)
    plt.savefig(os.path.join(figpath,'SCREAM_degraded_slow_fast_slab.png'))
    
    # compare protons SCREAM to NIST (slow method)
    nfront = SlabLayer(ntable,'SILICON DIOXIDE',15) # 15 mils of SiO2
    nslow = nfront.degraded_spectrum(pMeV,pFluence,fast=False)
    plt.figure()
    plt.loglog(pMeV,pFluence/2,'k-',pMeV,slow,'b-',pMeV,nslow,'r--')
    plt.xlabel('MeV')
    plt.ylabel('protons/cm$^2$/MeV')
    plt.title('Comparison of Scream and NIST tables %s' % repr(front),fontsize='smaller')
    plt.legend(['Incident/2','Scream','NIST'])
    plt.grid(True)
    plt.text(1e0,4e-8,'NIST drops to 0 @ ~3e3 MeV b/c end\n'
             'of degraded incident spectrum w/in table.\n' +
             'Scream continues b/c assumes non-degraded\n'+
             'spectrum beyond end of table.',verticalalignment='bottom',
             backgroundcolor='w')
    plt.savefig(os.path.join(figpath,'NIST_degraded_slow_fast_SCREAM_slab.png'))
    
    # demo electrons (slow method)
    efront = SlabLayer(etable,'SiO2',5) # 5 mils of SiO2
    eslow = efront.degraded_spectrum(eMeV,eFluence,fast=False)
    plt.figure()
    plt.loglog(eMeV,eFluence/2,'k-',eMeV,eslow,'b-')
    plt.xlabel('MeV')
    plt.ylabel('electrons/cm$^2$/MeV')
    plt.title('Electrons %s' % repr(efront),fontsize='smaller')
    plt.legend(['Incident/2','Slow Method'])
    plt.grid(True)
    plt.savefig(os.path.join(figpath,'SCREAM_electrons_degraded_slab.png'))
    
    plt.show()