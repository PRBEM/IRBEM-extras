"""
shiedled_part.py - provides wrapper class for a shielded part that allows
sharing of degraded spectra and energy->LET spectral conversion results 
& precomputed matrices

Public Globals:
    VERBOSE - report new cache computations
    DEFAULT_LET_GRID - default LET grid used by ShieldedPart

Public Functions:
    None

Public Classes:
    ShieldedPart - object that takes incident flux and produces SEE rate
"""

import numpy as np
from interpolating_matrix import interpolating_matrix

VERBOSE = True # print new cache calculations
_DEBUG = False

DEFAULT_LET_GRID = np.exp(np.linspace(np.log(1),np.log(110000),1002))


class ShieldedPart(object):
    """
    ShieldedPart(part,targetTable,shield=None,fast=True)
    part = EEEPart
    shield = ShieldLayer (None = no shielding)
    targetTable = RangeTable for target material
    fast - whether to use fast options for shielding and rate calculations
    
    The class can use extensive caching to generate and store precomputed 
     matrices, even allowing to skip steps in the calculation from:
      incident flux -A> degraded flux -B> LET flux -C> SEE rate (ions)
      incident flux -A> degraded flux -C> SEE rate (protons)
      degraded_spectrum provides the -A> operation
      LET_spectrum provides the -B> operation (and -A> if needed)
      see_rate provides the -C> operation (and -A> and/or -B> if needed)
    The functions and cache account for whether the part has ion or proton 
     SEEs using part.is_proton_part
    For maximum speed, call the earliest function in the chain with skip=True.
     This will allow internal tracking to avoid unnecessary array comparisons 
     and cache rebuilding. If intermediate results are needed, set skip=False.
    
    methods:
        degraded_spectrum - compute degraded spectrum for part's shield
        LET_spectrum - compute LET spectrum for target material
        see_rate = compute SEE rate for part
        clear_cache - clear/reset the cache (use when energy,Z,LET grids change)
    
    data:
        cache: dict (so that intermediate results can be accessed and re-used for other parts)
            (all fields will be present, uninitialized/cleared will be None)
            iflux -incident flux from last call to see_rate (NE,NZ)
            dflux - degraded flux from last call to see_rate (NE,NZ)
            LETflux - LET spectrum from last call to see_rate (NLET,)
            (entries below will be reatined between calls unless clear_cache is called. 
            This can result in the cache not being used if a request is made with 
            an incompatible Z, energy or LET argument)
            fluxZ - Z values used in most recent calculation
            icache - mapping from fluxZ to Z
            Z - array of Zs used in energy matrix and other matrices
            energy - energy matrix (NE,NZ)
            LET - LET grid (NLET,NLET)
            D - incident spectrum to degraded spectrum matrices(NE,NE,NZ)
            M - degraded spectrum to LET spectrum matrices (NLET,NE,NZ) (ions)
            MD - incident energy spectrum to LET spectrum matrices (NLET,NE,NZ) (ions)
            g - LET spectrum to rate vector. (NLET,) (ions)or protons (NE,) (protons)
            gM - degraded ion spectrum to SEE rate (NE,NZ) (ions)
            gMD - incident ion spectrum to SEE rate (NE,NZ) (ions)
            gD - incident proton spectrum to SEE rate (NE,NZ) NZ=1 (protons)
            

    Note: proton fluxes and energy will be broadcast to (NE,1) with Z=[1]
    """
    def __init__(self,part,targetTable,shield=None,fast=True):
        self.part = part
        self.targetTable = targetTable
        self.shield = shield
        self.fast = True
        self.cache = {}
        self.clear_cache()
    def __str__(self):
        return "%s behind %s" % (self.part,self.shield)
    def clear_cache(self,fluxes=True,grids=True,matrices=True):
        """
        clear cached last* data
        clear_cache(fluxes=True,grids=True,matrices=True)
        fluxes - iflux, dflux, LETflux 
        grids - energy, Z, LET, fluxZ, icache
        matrices - D, M, MD, g, gM, gMD, gD
        """
        if fluxes:
            for key in ['iflux','dflux','LETflux']:
                self.cache[key] = None
        if grids:
            for key in ['energy','Z','LET','fluxZ','icache']:
                self.cache[key] = None
                
        if matrices:
            for key in ['D','M','MD','g','gM','gMD','gD']:
                self.cache[key] = None
                
    def _makeDcache(self):
        """
            D = _makeDcache()
            make or return cache['D'] using cache['energy'] and cache['Z']
        """
        if self.cache['D'] is None:
            if VERBOSE:
                print('New D cache for %s' % self)
            energy = self.cache['energy']
            assert energy is not None,'energy is not in cache yet'
            Z = self.cache['Z']
            assert Z is not None,'Z is not in cache yet'
            if self.shield is None:
                D = np.stack([np.eye(len(energy))]*len(Z),axis=2) # identity matrix
            else:
                D = np.full((energy.shape[0],energy.shape[0],len(Z)),np.nan)
                for (iz,z) in enumerate(Z):
                    s = self.targetTable.Z2sym[z]
                    cache = self.shield.get_cache(species=s,Tfast=True) # on shield.table's energy grid
                    # interpolate onto energy grid used in this call, use log energy
                    Etable = self.shield.table.lookupEnergy(self.shield.material,species=s) # table's energy grid
                    cache = np.dot(interpolating_matrix(np.log(Etable),np.log(energy[:,iz])),
                                   np.dot(cache,interpolating_matrix(np.log(energy[:,iz]),np.log(Etable))))
                    D[:,:,iz] = cache
            self.cache['D'] = D
        return self.cache['D']
    def _makeMcache(self,DMkey='MD'):
        """
        matrix = _makeMcache(DMkey='MD')
        make or return cache['M'], cache['MD'] using cache['energy'], cache['Z'], and cache['LET']
        DMkey - 'M' or 'MD' which matrix to make/request
        M = _makeDcache('M') (also saves MD in cache)
        MD = _makeDcache('MD') (also saves M in cache)
        """
        if self.fast and (self.cache[DMkey] is not None):
            return self.cache[DMkey]
        D = self._makeDcache() # always get D cache, even if not needed
        if VERBOSE:
            print('New M & DM cache for %s' % self)

        energy = self.cache['energy'] # already confirmed in _makeDcache
        Z = self.cache['Z'] # already confirmed in _makeDcache
        LET = self.cache['LET']
        assert LET is not None,'LET is not in cache yet'
            
        # build matrices that convert from energy spectrum to LET spectrum
        # one per Z
        M = np.zeros((energy.shape[0],len(LET),len(Z)))
        MD = np.zeros(M.shape)
        for (iz,z) in enumerate(Z):
            M[:,:,iz] = self.targetTable.get_Mcache(energy[:,iz],LET,species=self.targetTable.Z2sym[z])
            MD[:,:,iz] = np.dot(M[:,:,iz],D[:,:,iz])
        self.cache['M'] = M
        self.cache['MD'] = MD
        
        return self.cache[DMkey]
    def _makegcache(self):
        """
        _makegcache() 
        make all the g, gM, gMD, gD caches possible for this part
        """
        if self.part.is_proton_part:
            if self.cache['g'] is None:
                if VERBOSE:
                    print('New g cache for %s' % self)
                cache = self.part.get_cache(energy=self.cache['energy'])
                self.cache['g'] = cache['g']
            self._makeDcache() # will try to make D
            if self.cache['D'] is not None:
                self.cache['gD'] = np.tensordot(self.cache['g'],self.cache['D'],([0],[0])) # (NE,NZ)
        else: # ion pat
            if self.cache['g'] is None:
                if VERBOSE:
                    print('New g cache for %s' % self)
                cache = self.part.get_cache(LET=self.cache['LET'])
                self.cache['g'] = cache['g']
            if self.cache['energy'] is not None:
                self._makeMcache() # will try to make D, M, and MD
                if (self.cache['gM'] is None)and (self.cache['M'] is not None) and (self.cache['g'] is not None):
                    if VERBOSE:
                        print('New gM cache for %s' % self)
                    self.cache['gM'] = np.tensordot(self.cache['g'],self.cache['M'],axes=([0],[0])) # dot along LET axis, (NE,NZ)
                if (self.cache['gMD'] is None)and (self.cache['D'] is not None) and (self.cache['gM'] is not None):
                    if VERBOSE:
                        print('New gMD cache for %s' % self)
                    gMD = np.zeros(self.cache['energy'].shape)
                    for iz in range(self.cache['energy'].shape[1]):
                        gMD[:,iz] = np.dot(self.cache['gM'][:,iz],self.cache['D'][:,:,iz])
                    self.cache['gMD'] = gMD
        
    
    def _prep_JZE(self,flux,Z,energy,useCache,DMkey,makeCache=None,icache=None):
        """
        (flux,Z,energy,useCache,icache) = self._prep_JZE(flux,Z,energy,useCache,DMkey,makeCache=None,icache=None)
        prepare flux, Z, energy arguments relative to cache and dimension conventions
        flux = (NE,) or (NE,NZ) ion spectrum, units are #/cm^2/s/MeV
        Z - scalar or (NZ,) array of Z values used to interpret ion flux spectrum        
            Z can be subset or re-ordered cache['Z'] and will still try to use cache
            for a proton part, omitted Z defaults to [1]
        energy - (NE,) or (NE,NZ) array of energies in MeV
        Energy and Z will be pulled from cache if None
        useCache - try to use cache
        icache - three options:
            None - energy and Z *have not* been prepared use with the cache
            True - energy and Z *have* been prepared use with the cache (retrieve actual icache from the cache)
            list - index array to relate columns of Z/energy to cache e.g., Z[i] = cache['Z'][icache[i]]
        DMkey - 'D', 'MD', or 'M' - cache key for matrix that will be used in fast mode
        makeCache - callable that builds the  cache matrix/matrices
        on output:
        flux - (NE,NZ)
        energy - (NE,) or (NE,NZ)
        useCache - can the cache be used (False if energy supplied but doesn't agree with cache)
        NZ - (NZ,), removes Z=1 if part.excludeH
        # for proton parts, Z will be cropped to [1], flux and energy will be cropped to shape (NE,1)
        """
        if makeCache is None:
            makeCache = lambda *args: None # do nothing function
            
        if useCache and (icache is None):
            icache = self.cache['icache'] # try to retrieve icache from cache, None otherwise
            
        iprot = None # index of protons in Z, energy, and flux arrays as 1-element list
        if self.part.is_proton_part and (Z is not None):
            iprot = np.where(np.atleast_1d(Z).ravel()==1)[0]
            assert len(iprot)==1,'Z=1 not provided (uniquely)'
            
        flux = np.atleast_1d(flux)        
        if self.part.excludeH and (Z is not None):
            # remove protons
            Z = np.atleast_1d(Z).ravel()
            if 1 in Z: # remove protons
                iz = (Z!=1)
                if (flux.ndim == 2) and (flux.shape[1] == len(Z)):
                    flux = flux[:,iz]
                if energy is not None:
                    energy = np.atleast_1d(energy)
                    if (energy.ndim == 2) and (energy.shape[1] == len(Z)):
                        energy = energy[:,iz]
                Z = Z[iz]
            
            
        if useCache and (icache is None): # do all the icache/fluxZ stuff
            if Z is None:
                Z = self.cache['Z']
                if (Z is None) and (self.part.is_proton_part):
                    Z = [1]
                assert Z is not None,'Z must be supplied when cache is empty'
            else:
                Z = np.atleast_1d(Z).ravel()
                if self.part.is_proton_part: # shrink wrap to Z=1 only
                    Z = np.array([1]) # NZ=1
                    
                # assemble energy if cache Z not equal to argument Z
                if (useCache and (self.cache['energy'] is not None) and 
                    (self.cache[DMkey] is not None) and 
                    (self.cache['Z'] is not None) and 
                    not np.array_equal(Z,self.cache['Z'])):
                    cenergy = np.zeros((self.cache['energy'].shape[0],len(Z))) # energies from cache
                    icache = [None]*len(Z)
                    for (iz,z) in enumerate(Z):
                        jz = np.where(Z==z)
                        if len(jz[0]) == 0:
                            raise Exception('No cached energy for Z=%d' % z)
                        jz = jz[0][0]
                        icache[iz]= jz
                        cenergy[:,iz] = self.cache['energy'][:,jz]
                    if (energy is None):
                        energy = cenergy
                    elif not np.array_equal(energy,cenergy):
                        useCache = False # cannot use cache, energies are different
            
            if energy is None:
                energy = self.cache['energy']
                icache = list(range(len(Z)))
                assert energy is not None,'energy must be supplied when cache is empty'
                if useCache and (self.cache[DMkey] is None):
                    # make cache now to avoid extra array checking later
                    if self.cache['Z'] is None:
                        self.cache['Z'] = Z 
                        makeCache()
                    elif np.array_equal(self.cache['Z'],Z):
                        makeCache()
            energy = np.atleast_1d(energy)
            
            if energy.ndim == 1:
                energy = energy.reshape((len(energy),1)) # (NE,1)
            elif iprot is not None:
                energy = energy[:,iprot] # (NE,1)
        else:
            if Z is None:
                Z = self.cache['fluxZ']
            energy = self.cache['energy'][:,icache] # correct shape &c
            makeCache()

        if useCache and (self.cache[DMkey] is None):
            if ((self.cache['energy'] is None) and (self.cache['Z'] is None)):
                # needed cache grids and matrices are empty
                if (energy.shape[1]==1) and (len(Z)>1):
                    energy = np.hstack([energy]*len(Z)) # expand energy to length of Z
                icache = list(range(len(Z)))
                self.cache['energy'] = energy
                self.cache['Z'] = Z
                makeCache()
            elif (np.array_equal(energy,self.cache['energy']) and np.array_equal(Z,self.cache['Z'])):
                # cache matches energy,Z
                icache = list(range(len(Z)))
                makeCache()

        if useCache and (icache is None):
            icache = self.cache['icache']
            

        flux = np.atleast_1d(flux)
        if flux.ndim == 1:
            flux = flux.reshape((len(flux),1)) # NE,1
        elif iprot is not None:
            flux = flux[:,iprot] # (NE,1)
            
        assert np.array_equal(flux.shape,energy.shape), 'flux and energy must have same length'

        if (iprot is not None) and (len(Z)>1):
            Z = Z[iprot]
            
        return (flux,Z,energy,useCache,icache)
        
    def _prep_LET(self,LET,useCache):
        """
        (LET,useCache) = self._prep_LET(LET,useCache)
        prepare LET and/or cache
        if None gets LET from cache or generates from default array
        useCache - can the cache be used (false if LET supplied by doesn't match cache)
        """
        if useCache and (LET is not None) and (self.cache['LET'] is not None):
            useCache = np.array_equal(LET,self.cache['LET'])
        
        if LET is None:
            LET = self.cache['LET']
        if LET is None: # supply from default array
            LET = DEFAULT_LET_GRID

        if useCache and (self.cache['LET'] is None):
            self.cache['LET'] = LET # store
            
        return (LET,useCache)

    def degraded_spectrum(self,flux,energy=None,Z=None,icache=None):
        """
        dflux = degraded_spectrum(flux,energy=None,Z=None,icache=None)
        compute degraded flux
        flux = (NE,) or (NE,NZ) ion spectrum, units are #/cm^2/s/MeV
        energy - (NE,) or (NE,NZ) array of energies in MeV
        Z - scalar or (NZ,) array of Z values used to interpret ion flux spectrum        
            for a proton part, omitted Z defaults to [1]
        Energy and Z will be pulled from cache if None
        dflux - same shape as flux, degraded flux in #/cm^2/s/MeV
            (except if proton part, then dflux will never have more than 1 column)
        icache - three options:
            None - energy and Z *have not* been prepared use with the cache
            True - energy and Z *have* been prepared use with the cache (retrieve actual icache from the cache)
            list - index array to relate columns of Z/energy to cache e.g., Z[i] = cache['Z'][icache[i]]

        if self.fast, will update cache fluxes and build caches if possible
        (e.g., if energy and Z are None or match cached values or there's no cache)
        """        

        indims = np.atleast_1d(flux).ndim

        self.clear_cache(grids=False,matrices=False) # clear fluxes only
        (flux,Z,energy,useCache,icache) = self._prep_JZE(flux,Z,energy,self.fast,'D',self._makeDcache,icache)

        if energy.shape[1] < len(Z):
            dE = 0 # do not advance energy index
        else:
            dE = 1  # advance energy index
        
        dflux = np.zeros(flux.shape)
        if useCache:
            self._makeDcache() # ensure D is in cache
            for (i,z) in enumerate(Z):
                dflux[:,i] = np.dot(self.cache['D'][:,:,icache[i]],flux[:,i])
        else:
            for (i,z) in enumerate(Z):
                sp = self.shield.table.Z2sym[z]
                dflux[:,i] = self.shield.degraded_spectrum(energy[:,i*dE],flux[:,i],fast=True,species=sp)
            
        self.clear_cache(grids=False,matrices=False)
        self.cache['iflux'] = flux
        self.cache['dflux'] = dflux
        self.cache['fluxZ'] = Z
        self.cache['icache'] = icache
        # fix dimensions of dflux to match dimensions of flux
        if indims == 0:
            dflux = np.asscalar(dflux)
        elif indims == 1:
            dflux = dflux.reshape((len(dflux),))
        # flux came in with mulitple Zs, leave alone
        return dflux
        
    def LET_spectrum(self,flux,energy=None,LET=None,Z=None,degraded=False,skip=True,icache=None):
        """
        LETflux = LET_spectrum(flux,energy=None,LET=None,Z=None,degraded=False,skip=True,icache=None)
        convert energy spectra into LET spectrum
        flux = ion spectrum, units are #/cm^2/s/MeV
        energy - (NE,) or (NE,NZ) array of energies in MeV  (None to use cache)
        LET - (NLET,) array of LETs in MeV/(g/cm^2) (None to use cache)
        Z - (NZ,) array of Z values used to interpret ion flux spectrum (None to use cache)
            for a proton part, omitted Z defaults to [1]
        degraded - has the proton/ion spectrum already been passed through shielding?
        skip - only compute LETflux rate skipping explicit middle steps if possible (only works if self.fast)
        icache - three options:
            None - energy and Z *have not* been prepared use with the cache
            True - energy and Z *have* been prepared use with the cache (retrieve actual icache from the cache)
            list - index array to relate columns of Z/energy to cache e.g., Z[i] = cache['Z'][icache[i]]
        LETflux - (NLET,) LET spectrum, #/cm^2/s/(MeV/(g/cm^2))        

        if self.fast, will update cache fluxes and build caches if possible
        (e.g., if energy and Z are None or match cached values or there's no cache)

        """

        
        if degraded:
            DMkey = 'M' # for degraded flux, only need cache M
        else:
            DMkey = 'MD' # for incident flux, need cache D and M or DM
        makeCache = lambda *args: self._makeMcache(DMkey)
        
        iflux = None
        dflux = None
        LETflux = None

        self.clear_cache(grids=False,matrices=False) # clear fluxes only
        (LET,useCache) = self._prep_LET(LET,self.fast)
        (flux,Z,energy,useCache,icache) = self._prep_JZE(flux,Z,energy,useCache,DMkey,makeCache,icache)
        
        if degraded:
            dflux = flux
        else:
            iflux = flux
        
        if useCache:
            makeCache() # ensure DMkey is in cache
            LETflux = np.zeros((len(LET),))
            if skip or degraded:
                for (iz,z) in enumerate(Z):
                    LETflux += np.dot(self.cache[DMkey][:,:,iz],flux[:,iz]) # flux is iflux or dflux, as needed by DMkey
            else:
                if not degraded:
                    dflux = self.degraded_spectrum(iflux,energy=energy,Z=Z,icache=icache)
                for (iz,z) in enumerate(Z):
                    LETflux += np.dot(self.cache['M'][:,:,iz],dflux[:,iz])
        else:
            if not degraded:
                dflux = self.degraded_spectrum(iflux,energy=energy,Z=Z,icache=icache)
            LETflux = self.targetTable.LET_spectrum(dflux,E,LET,Z=Z,fast=self.fast)
                
        self.clear_cache(grids=False,matrices=False)
        self.cache['iflux'] = iflux
        self.cache['dflux'] = dflux
        self.cache['fluxZ'] = Z
        self.cache['icache'] = icache
        self.cache['LET'] = LET
        self.cache['LETflux'] = LETflux
        return LETflux
        
    def see_rate(self,flux,energy=None,LET=None,Z=None,degraded=False,skip=True,icache=None):
        """
        rate = see_rate(flux,energy=None,LET=None,Z=None,degraded=False,skip=True)
        flux = proton, ion or LET spectrum
            for proton and ion spectra, units are #/cm^2/s/MeV
            for LET spectrum, units are #/cm^2/s/(MeV/(g/cm^2))
        energy - (NE,) or (NE,NZ) array of energies in MeV
        LET - (NLET,) array of LETs in MeV/(g/cm^2)
        Z - (NZ,) array of Z values used to interpret ion flux spectrum (if omitted, assumes Z=[1])
        degraded - has the proton/ion spectrum already been passed through shielding?
        skip- only compute SEE rate skipping explicit middle steps if possible (only works if self.fast)
        icache - three options:
            None - energy and Z *have not* been prepared use with the cache
            True - energy and Z *have* been prepared use with the cache (retrieve actual icache from the cache if needed)
            list - index array to relate columns of Z/energy to cache e.g., Z[i] = cache['Z'][icache[i]]

        rate: SEE rate: #/s 
        
        several ways to call:
        rate = see_rate(flux,energy=energy) - flux is (NE,) or (NE,1) incident proton flux
        rate = see_rate(flux,energy=energy,Z=Z) - flux is (NE,NZ) incident ion flux
        rate = see_rate(flux,energy=energy,degraded=True) - flux is (NE,) proton flux behind shielding
        rate = see_rate(flux,energy=energy,Z=Z,degraded=True) - flux is (NE,NZ) ion flux behind shielding
        rate = see_rate(flux,LET=LET) - flux is LET spectrum (NLET,) *behind shielding* (degarded ignored)
        """
        
        if (energy is None) == (LET is None):
            raise Exception('Must supply exactly one of energy or LET')
            
            
        is_energy_spectrum = (LET is None) # is flux an energy spectrum?

        (LET,useCache) = self._prep_LET(LET,self.fast)
        
        iflux = None
        dflux = None
        LETflux = None
        r = None
        
        if self.part.is_proton_part: # compute r from dflux
            assert is_energy_spectrum,'Proton parts require energy spectrum, not LET spectrum'
            (flux,Z,energy,useCache,icache) = self._prep_JZE(flux,Z,energy,useCache,'g',self._makegcache,icache=icache)
            if degraded:
                dflux = flux
                if useCache:
                    r = np.dot(self.cache['g'].ravel(),dflux.ravel())
            else:
                iflux = flux
                if useCache and skip:
                    r = np.dot(self.cache['gD'].ravel(),iflux.ravel())
                else:
                    dflux = self.degraded_spectrum(flux,energy=energy,Z=Z,icache=icache)
            if r is None:
                if useCache:
                    r = np.dot(self.cache['g'].ravel(),dflux.ravel())
                else:
                    r = self.part.see_rate(dflux[:,0],energy=energy[:,0],fast=self.fast)
        else: # ion part, compute r from LETflux
            if is_energy_spectrum: # make LETflux
                (flux,Z,energy,useCache,icache) = self._prep_JZE(flux,Z,energy,useCache,'g',self._makegcache,icache=icache)
                if useCache and skip:
                    if degraded:
                        dflux = flux
                        r = (self.cache['gM']*dflux).sum() # gM .* dflux, sum all elements
                    else:
                        iflux = flux
                        r = (self.cache['gMD']*iflux).sum() # gMD .* iflux, sum all elements
                else:
                    LETflux = self.LET_spectrum(flux,energy=energy,LET=LET,Z=Z,degraded=degraded,skip=skip,icache=icache)
                    # retrieve intermediate results from cache, if they were generated
                    iflux = self.cache['iflux']
                    dflux = self.cache['dflux']
            else: # format LETflux
                LETflux = np.atleast_1d(flux).ravel() # (NLET,)
                assert np.array_equal(LET.shape,LETflux.shape),'LET and flux have inconsistent shapes'
        if r is None: # compute r from LETflux
            if useCache:
                self._makegcache()
                r = np.dot(self.cache['g'].ravel(),LETflux.ravel())
            else:
                r = self.part.see_rate(LETflux,LET=LET,fast=self.fast)
            
        self.clear_cache(grids=False,matrices=False)
        self.cache['iflux'] = iflux
        self.cache['dflux'] = dflux
        self.cache['fluxZ'] = Z
        if useCache:
            if icache is not None:
                self.cache['icache'] = icache
            self.cache['LET'] = LET
        self.cache['LETflux'] = LETflux
        return r

if __name__ == '__main__':
    
    import datetime as dt
    import matplotlib.pyplot as plt
    from rpp import ProtonPart, IonPartRPP
    from degraded_spectra import SRIMTable,SphereLayer,SlabLayer
    
    import degraded_spectra
    degraded_spectra.VERBOSE = VERBOSE
    
    import rpp
    rpp.VERBOSE = VERBOSE
    
    AlTable = SRIMTable('SRIM_Al.csv')
    SiTable = SRIMTable('SRIM_Si.csv')
    # assume AlTable and SiTable have same species
    syms = AlTable.Species
    Z = np.array([AlTable.data[s]['Z'] for s in syms])
    if _DEBUG:
        Z = Z[:8] # only do subset of Zs (faster)
    
    slab = SlabLayer(AlTable,'Al',100,units='mils')
    sphere = SphereLayer(AlTable,'Al',100,units='mils')
    
    
    # make up fake flux
    # protons
    pEgrid = np.exp(np.linspace(np.log(1e-1),np.log(1e5),1002))
    pflux = np.exp(-pEgrid/10)
    t = np.linspace(0,1,20)    
    Egrid = np.zeros((len(pEgrid),len(Z)))
    MeVperNuc = np.stack([pEgrid]*len(Z),axis=1)  # all have same MeV/nuc
    flux = np.zeros((len(pEgrid),len(Z),len(t)))
    for (iz,z) in enumerate(Z):
        s = syms[iz]
        A = AlTable.data[s]['A']
        Egrid[:,iz] = MeVperNuc[:,iz]*A
        for (it,ti) in enumerate(t):
            flux[:,iz,it] = (1+np.sin(ti*2*np.pi))*pflux/z**2
    
    ppart = ProtonPart(E0=30,W=10,S=1,slim=1e-11,volumes=1)
    ipart = IonPartRPP(1,L0=500,W=30e3,S=1,a=44.7124,b=44.7124,volumes=100)
    
    markers = ['s','d','o','+','.','x',None]
    
    for shield in [sphere,slab]:
        #p0 rate = see_rate(flux,energy=energy,skip=False) - flux is (NE,) or (NE,1) incident proton flux
        #p1 rate = see_rate(flux,energy=energy,Z=Z,skip=False) - flux is (NE,NZ) incident ion flux
        #p2 rate = see_rate(flux,energy=energy,degraded=True,skip=False) - flux is (NE,) proton flux behind shielding
        #p3 rate = see_rate(flux,energy=energy,Z=Z,degraded=True,skip=False) - flux is (NE,NZ) ion flux behind shielding
        #p4 rate = see_rate(flux,energy=energy,skip=True) - flux is (NE,NZ) ion flux behind shielding
        #p5 rate = see_rate(flux,energy=energy,degraded=True,skip=True) - flux is (NE,NZ) ion flux behind shielding
        #p6 - separate calles to degraded_spectrum, then see_rate
        
        #i0 rate = see_rate(flux,energy=energy,Z=Z,skip=False) - flux is (NE,NZ) incident ion flux
        #i1 rate = see_rate(flux,energy=energy,Z=Z,degraded=True,skip=False) - flux is (NE,NZ) ion flux behind shielding
        #i2 rate = see_rate(flux,LET=LET,skip=False) - flux is LET spectrum (NLET,) *behind shielding* (degarded ignored)
        #i3 rate = see_rate(flux,energy=energy,Z=Z,skip=true) - flux is (NE,NZ) incident ion flux
        #i4 rate = see_rate(flux,energy=Egrid,Z=Z,degraded=True,skip=True)
        #i5 rate = see_rate(LETflux,LET=LET,skip=True)
        #i6 - separate calles to degraded_spectrum, LET_spectrum then see_rate
    
        p0 = ShieldedPart(ppart,SiTable,shield)
        p1 = ShieldedPart(ppart,SiTable,shield)
        p2 = ShieldedPart(ppart,SiTable,shield)
        p3 = ShieldedPart(ppart,SiTable,shield)
        p4 = ShieldedPart(ppart,SiTable,shield)
        p5 = ShieldedPart(ppart,SiTable,shield)
        p6 = ShieldedPart(ppart,SiTable,shield)
    
        i0 = ShieldedPart(ipart,SiTable,shield)
        i1 = ShieldedPart(ipart,SiTable,shield)
        i2 = ShieldedPart(ipart,SiTable,shield)
        i3 = ShieldedPart(ipart,SiTable,shield)
        i4 = ShieldedPart(ipart,SiTable,shield)
        i5 = ShieldedPart(ipart,SiTable,shield)
        i6 = ShieldedPart(ipart,SiTable,shield)
        
        prates = np.full((len(t),7),np.nan)
        irates = np.full((len(t),7),np.nan)
        
        iprot = np.where(Z)[0][0]
        for (it,ti) in enumerate(t):

            print('time %d/%d i0 %s' % (it+1,len(t),dt.datetime.now()))
            irates[it,0] = i0.see_rate(flux[:,:,it],energy=Egrid,Z=Z,skip=False)
            dflux = i0.cache['dflux']
            print('time %d/%d i1 %s' % (it+1,len(t),dt.datetime.now()))
            irates[it,1] = i1.see_rate(dflux,energy=Egrid,Z=Z,degraded=True,skip=False)
            LETflux = i1.cache['LETflux']
            LET = i1.cache['LET']
            print('time %d/%d i2 %s' % (it+1,len(t),dt.datetime.now()))
            irates[it,2] = i2.see_rate(LETflux,LET=LET,skip=False)
            print('time %d/%d i3 %s' % (it+1,len(t),dt.datetime.now()))
            
            irates[it,3] = i3.see_rate(flux[:,:,it],energy=Egrid,Z=Z,skip=True)
            print('time %d/%d i4 %s' % (it+1,len(t),dt.datetime.now()))
            irates[it,4] = i4.see_rate(dflux,energy=Egrid,Z=Z,degraded=True,skip=True)
            print('time %d/%d i5 %s' % (it+1,len(t),dt.datetime.now()))
            irates[it,5] = i5.see_rate(LETflux,LET=LET,skip=True)
            print('time %d/%d i6 %s' % (it+1,len(t),dt.datetime.now()))
            dflux = i6.degraded_spectrum(flux[:,:,it],energy=Egrid,Z=Z)
            LETflux = i6.LET_spectrum(dflux,energy=Egrid,LET=LET,Z=Z,degraded=True)
            irates[it,6] = i6.see_rate(LETflux,LET=LET)
            
            print('time %d/%d p0 %s' % (it+1,len(t),dt.datetime.now()))
            prates[it,0] = p0.see_rate(flux[:,iprot,it],energy=Egrid[:,iprot],skip=False)
            print('time %d/%d p1 %s' % (it+1,len(t),dt.datetime.now()))
            prates[it,1] = p1.see_rate(flux[:,:,it],energy=Egrid,Z=Z,skip=False)
            dflux = p1.cache['dflux']
            E = p1.cache['energy']
            z = p1.cache['Z']
            jprot = np.where(z==1)[0][0]
            print('time %d/%d p2 %s' % (it+1,len(t),dt.datetime.now()))
            prates[it,2] = p2.see_rate(dflux[:,jprot],energy=E[:,jprot],degraded=True,skip=False)
            print('time %d/%d p3 %s' % (it+1,len(t),dt.datetime.now()))
            prates[it,3] = p3.see_rate(dflux,energy=E,Z=z,degraded=True,skip=False)
            print('time %d/%d p4 %s' % (it+1,len(t),dt.datetime.now()))
            
            prates[it,4] = p4.see_rate(flux[:,iprot,it],energy=Egrid[:,iprot],skip=True)
            print('time %d/%d p5 %s' % (it+1,len(t),dt.datetime.now()))
            prates[it,5] = p5.see_rate(dflux[:,jprot],energy=E[:,jprot],degraded=True,skip=True)
            print('time %d/%d p6 %s' % (it+1,len(t),dt.datetime.now()))
            dflux = p6.degraded_spectrum(flux[:,:,it],energy=Egrid,Z=Z)
            prates[it,6] = p6.see_rate(dflux,energy=Egrid,Z=Z,degraded=True)
    
    
        plt.figure()
        H = plt.semilogy(t,prates)
        for (ih,h) in enumerate(H):
            h.set_marker(markers[ih])
        plt.legend(['p0: H only','p1: all Z','p2: H only, degraded','p3: all Z degraded','p4: rate only'])
        plt.xlabel('time')
        plt.ylabel('SEE/s')
        plt.title('Proton SEE rates %s ' % p0,fontsize='x-small')
    
        plt.figure()
        H = plt.semilogy(t,irates)
        for (ih,h) in enumerate(H):
            h.set_marker(markers[ih])
        plt.legend(['i0: incident','i1: degraded','i2: LET','i3: rate_only'])
        plt.xlabel('time')
        plt.ylabel('SEE/s')
        plt.title('Ion SEE rates %s ' % i0,fontsize='xx-small')
        
    plt.show()
