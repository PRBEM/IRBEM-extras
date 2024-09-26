"""
approximate translation of IRBEM open diffusion code (odc) utility library from Matlab
"""

# todo:
# deal with species
# 180/pi and pi/180 -> np.radians, np.degrees

import numpy as np
# global shared constants
#SI/mks constants
mks = {}
mks['e']= 1.602176487e-19 # Coulombs per fundamental charge
mks['c'] = 299792458 # m/s
mks['eV'] = mks['e'] # Joules
mks['keV'] = mks['eV']*1e3 # Joules
mks['MeV'] = mks['eV']*1e6 # Joules
mks['GeV'] = mks['eV']*1e9 # Joules
mks['epsilon0'] = 8.854187817e-12 # % F/m = C^2 s^2  / kg / m^3 - permitivity of free space
mks['mu0'] = 4*np.pi*1e-7 # % H/m = N/A^2 = kg m/C^2 - permeability of free space
mks['R_E'] = 6371.2e3 # m, surface-area-averaged reference value of Earth Radius
mks['electron'] = {}
mks['electron']['q'] = -mks['e'] # charge
mks['electron']['m0'] = 9.10938215e-31 # rest mass, kg
mks['proton'] = {}
mks['proton']['q'] = mks['e'] # charge
mks['proton']['m0'] = 1.672621637e-27 # rest mass, kg

# Schulz & Lanzerotti constants
SL = {};
  
SL['T0']  = 1+1/(2*np.sqrt(3))*np.log(2+np.sqrt(3)) # S&L 1.28a
SL['T1']  = np.pi/6*np.sqrt(2) # S&L 1.28b
SL['B0']  = 31e3*1e-9 # T dipole field at lambda=0,L=1 (S&L value)
SL['a']   = 6371.2e3 # Earth Radius, meters (S&L value)
SL['Q0']  = -27.1266694 # S&L 1.78c
SL['Q1']  = -90*SL['T1'] # S&L 1.78b
SL['Qp1'] = (15/2)*(9*SL['T0'] - 41*SL['T1']) # S&L 1.78d
SL['Y0']  = 2*SL['T0'] # S&L 1.31 (limiting case in subsequent text)

constants = {'mks':mks,'SL':SL}

def alphaL2K(alpha,L,Bunit='nT'):
    """
    returns K in RE*sqrt(nT) for given equatorial pitch angle
    alpha, in deg, and L shell, in RE
    Bunit - optional, if 'G' assumes K in RE*sqrt(G). Otherwise K in RE*sqrt(nT)
    global odc_constants
    """
    
    y = np.sin(np.radians(alpha))
    if Bunit == 'nT':
        Bm = SL['B0']*1e9/L**3./y**2 # mirror field, nT
    elif Bunit=='G':
        Bm = (SL['B0']*1e9/L**3./y**2 )/1e5 # nT to G
    K = SL_Y(y)*L*np.sqrt(Bm)
    return K

global _KL2alpha_table # private persistent variable
_KL2alpha_table = None # private persistent variable
def KL2alpha(K,L,Bunit='nT'):
    """
    returns equatorial pitch angle alpha in degrees
    for given K, L. K in RE*sqrt(nT)
    Bunit - optional, if 'G' returns K in RE*sqrt(G). Otherwise K in
    RE*sqrt(nT)
    Not sure if table should be a dictionary like I've set it up for now, or not. 
    """
    global _KL2alpha_table
    if _KL2alpha_table is None:
        table = {}
        dy = 1e-4
        #Start at dy:go in steps of dy: up to 1
        table['y'] = np.arange(dy,1+dy/2,dy)
        table['Yy'] = SL_Y(table['y'])/table['y']
    if Bunit=='G':
        K = K*np.sqrt(1e5) # RE*sqrt(G) to RE*sqrt(nT)
    # K*sqrt(L)/sqrt(B0) = Y(y)/y
    KLB = K*np.sqrt(L)/np.sqrt(SL['B0']*1e9)
    y = np.interp(KLB,_KL2alpha_table['y'],_KL2alpha_table['Yy'],KLB)
    alpha = np.arcsin(y)
    return alpha

def fce2B(fce):
    """
    returns B in nT for electron gyro given in Hz
    shorthand
    """
    electron = mks['electron']
    B = fce*2*np.pi*electron['m0']/np.abs(electron['q'])*1e9 # nT
    return B

def SL_T(y, returndT = False):
    """
    Schulz & Lanzerotti's T(y)
    and dT = dT/dy if requested
    """
    T = SL['T0']-0.5*(SL['T0']-SL['T1'])*(y+np.sqrt(y)) # 1/4 bounce integral of 1
    if returndT:
        dT = -0.5*(SL['T0']-SL['T1'])*(1+1./np.sqrt(y)/2) # dT/dy
        return T, dT
    else:
        return T

def SL_Q(y):
    """
    computes Q(y) from S&L, 1.79
    """
    Q = SL['Q0']+(2*SL['Q1'] - 2*SL['Q0'] - (1/4)*SL['Qp1'])*y**4 + (SL['Q0']-SL['Q1']+(1/4)*SL['Qp1'])*y**8 # Q(y) from S&L, 1.79
    return Q

def SL_D(y):
    """
    computes D(y) from S&L, 1.36
    """
    D = SL_T(y)/2-SL_Y(y)/12
    return D

def SL_Y(y):
    """
    computes Y(y) from S&L 1.31
    """
    Y = 2*(1-y)*SL['T0'] + (SL['T0']-SL['T1'])*(y*np.log(y) + 2*y-2*y**0.5) # S&L 1.31
    if np.isscalar(y):
        if y == 0:
            Y = SL['Y0']
    else:
        Y[y==0] = SL['Y0']
    return Y

def dipole_mirror_latitude(alpha0,units = 'deg'):
    """
    mirror_lat = dipole_mirror_latitude(alpha0)
    compute dipolar mirror latitude of particle with equatorial pitch angle
    alpha0, angles in degrees
    mirror_lat = dipole_mirror_latitude(alpha0,'rad')
    angles in radians
    """
    
    if np.isscalar(alpha0):
        alpha0 = np.array([alpha0])
        return dipole_mirror_latitude(alpha0,units=units)[0]
    
    if units.lower().startswith('d'):
        torad = np.pi/180
    elif units.lower().startswith('r'):
        torad = 1
    else:
        raise Exception('do not know what units were ment, the accpetible unit calls are deg or rad')

    sina0 = np.sin(alpha0*torad)
    if sina0==0:
        mirror_lat = np.inf
        return mirror_lat

    # note, below fixes error in Shprits thesis, eqn  F13
    # which has sina0**2. Instead use sina0**4 from Shprits 2006 eqn 10.
    mirror_lat = np.zeros(len(sina0))*np.nan
    for i in range(len(sina0)):
        Px = [1, 0, 0, 0, 0, 3*sina0[i]**4, -4*sina0[i]**4]
        xroots = np.roots(Px)
        #xroot = xroots((imag(xroots)==0) & (real(xroots)>0))
        xroot = xroots[(np.abs(xroots.imag)<1e-30) & (xroots.real>0)]        
        mirror_lat[i] = np.degrees(np.arccos(np.sqrt(xroot))) # mirror latitude
    return mirror_lat

global _maglat_table
_maglat_table = None
def BB0toMagLat(BB0, unit = 'deg'):
    """
    maglat = BB0toMagLat(BB0) # returns maglat in deg
    maglat = BB0toMagLat(BB0,'rad') # returns maglat in radians
    
    BB0 is B / Bmin
    where B is the local magnetic field strength
    and Bmin is the minimum (i.e., equatorial) magnetic field strength on the
    same field line.
    
    maglat is the unsigned dipole latitude
    """
    
    # we invert the B/Bmin vs maglat equation via a look-up table
    # we use a persistent variable so we don't have to regenerate the table
    # for each call to this function
    
    global _maglat_table
    if _maglat_table is None:
        _maglat_table = {}
        _maglat_table['deg'] = np.arange(0,90.0005,0.001)# matlab: (0:0.01:90)
        _maglat_table['rads'] = np.radians(_maglat_table['deg'])
        # this next bit is the expression for B/Bmin vs magnetic latitude for a dipole
        _maglat_table['bb0'] = (1+3*np.sin(_maglat_table['rads'])**2)**(1/2)/np.cos(_maglat_table['rads'])**6

    if unit.lower().startswith('d'):
        maglat = np.interp(BB0,_maglat_table['deg'], _maglat_table['bb0'])
    elif unit.lower().startswith('r'):
        maglat = np.interp(BB0,_maglat_table['rads'], _maglat_table['bb0'])
    else:
        raise Exception('Unknown unit ' + unit + ' please use either deg or rad');

    return maglat

def DriftPeriod(Species,Energy,PitchAngle,L, unit = 'deg'):
    """ 
    calculates particle drift period (seconds), dipole
    Species: 'e' for electrons, 'p' for protons
    Energy: in MeV
    PitchAngle: in degrees (equatorial pitch angle)
    L: dimensionless dipole L value
    """

    Species = SelectSpecies(Species)
    m0 = mks[Species]['m0']
    q = mks[Species]['q']
    c = mks['c'] # m/s    
    a = SL['a'] # RE, meters
    B0 = SL['B0'] # T

    gamma,v,m = MeVtogamma(Energy,Species)
    if unit.lower().startswith('d'):
        y = np.sin(np.radians(PitchAngle))
    elif unit.lower().startswith('r'):
        y = np.sin(PitchAngle)
    Ty = SL_T(y)
    Dy = SL_D(y)

    # S&L eq 1.35, w/o minus sign, assume always want positive drift velocity
    f = (3*L/2/np.pi/gamma)*(gamma**2-1)*(c/a)**2*(m0*c/abs(q)/B0)*(Dy/Ty)/c # extra 1/c for SI
    Td = 1.0/f # seconds
    return Td

def BouncePeriod(Species,Energy,PitchAngle,L,unit='deg'):
    """
    function Tb = BouncePeriod(Species,Energy,PitchAngle,L,unit='deg')
    calculates particle bounce period Tb (seconds), dipole
    Species: 'e' for electrons, 'p' for protons
    Energy: in MeV
    PitchAngle: in degrees or radians (equatorial pitch angle)
    unit: specify 'rad' for PitchAngle in radians
    L: dimensionless dipole L value
    """

    a = SL['a'] # Earth Radius, meters
    #y = sind(PitchAngle);
    if unit.lower().startswith('d'):
        y = np.sin(np.radians(PitchAngle))
    elif unit.lower().startswith('r'):
        y = np.sin(PitchAngle)
    (gamma,v,m) = MeVtogamma(Energy,Species)
    Ty = SL_T(y)
    Tb = 4*L*a/v*Ty
    return Tb

def GyroPeriod(Species,Energy,MagLat,L):
    """
    function Tg = GyroPeriod(Species,Energy,MagLat,L,unit='deg')
    calculates particle gyro period Tg (seconds), dipole
    Species: 'e' for electrons, 'p' for protons
    Energy: in MeV
    MagLat: in degrees
    L: dimensionless dipole L value
    
    Calculation done in .* so that vector and matrix input can work,
    if all inputs are the same size (or scalar)
    """

    Species = SelectSpecies(Species)
    q = mks[Species]['q'] # C
    (gamma,v,m) = MeVtogamma(Energy,Species)

    B = dipoleB(L,MagLat)/1e9 # T

    f = abs(q)*B/(2*np.pi*m) # no "c" in denominator in SI units
    Tg = 1./f
    return Tg

def dipoleB(L,MagLat,phi_deg,nargout=1):
    """
    B = dipoleB(L,MagLat,phi_deg,nargout=1)
    (B,Bx,By,Bz,X,Y,Z) = dipoleB(L,MagLat,phi_deg,nargout=7)
    computes magnitude of dipole field, nT
    returns components and positions (in RE) if requested
    MagLat: in degrees
    L: dimensionless dipole L value
    phi_deg: azimuth angle, degrees
    """

    Beq = SL['B0']/L**3*1e9 # T to nT
    smlat = np.sin(np.radians(MagLat))
    cmlat = np.cos(np.radians(MagLat))
    cmlat6 = cmlat**6
    B = Beq*np.sqrt(1+3*smlat**2)/cmlat6
    if nargout==1:
        return B

    phi = np.radians(phi_deg)
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    
    # angular part of Bx, By, Bz
    Btmp = Beq/cmlat6
    Bx = -3*cphi*cmlat*smlat*Btmp
    By = -3*sphi*cmlat*smlat*Btmp
    Bz = -1.*(3*smlat**2 - 1)*Btmp
    
    R = L*cmlat**2
    X = R*cmlat*cphi
    Y = R*cmlat*sphi
    Z = R*smlat
    return B,Bx,By,Bz,X,Y,Z

def SelectSpecies(Species):
    """
    convert various alternatives to standard species name
    """
    if Species.lower() in ('e','e-','electron','beta'):
        Species = 'electron'
    elif Species.lower() in ('p','p+','h+','proton','h','hydrogen'):
        Species = 'proton'
    else:
        print('Unknown Species')
        
    return Species

def MeVtogamma(MeV,species):
    """
    compute gamma, v, m
    given energy in MeV and species 'e','p', etc
    v in m/s
    m in kg
    """
    species = SelectSpecies(species)
    m0 = mks[species]['m0']

    W = MeV*mks['MeV'] # Energy, Joules

    gamma = 1+W/(m0*mks['c']**2) # relativistic factor   
    v=mks['c']*np.sqrt(1-gamma**-2)
    m = m0*gamma

    return gamma,v,m

def dipoleIJ(L,y,MeV=None,species=None):
    """
    I = util.dipoleIJ(L,y);
    (I,J) = util.dipoleIJ(L,y,MeV,species)
    computes I = L*Y(y) (in RE)
    and, if requested J = 2*p*I (in RE*MeV/c)
    for species ('e' or 'p')
    energy in MeV
    L: dimensionless dipole L value
    and y = sin(alpha_equatorial)
    """
    
    I = L*SL_Y(y)
    
    if MeV is None:
        return I
    
    (gamma,v,m) = MeVtogamma(MeV,species)
    p = m*v/mks['MeV']*mks['c'] # MeV/c
    J = 2.*p*I # RE*MeV/c
        
    return (I,J)

def  MBtoMeV(M,B,alpha,species):
    """
    given M in MeV/G, and B in nT, alpha in degrees, and species
    returns energy, in MeV
    """
    Bm = B/np.sin(np.radians(alpha))**2 # Bmirror
    species = SelectSpecies(species)
    m0 = mks[species]['m0']
    c = mks['c']

    p2 = 2*m0*Bm*M # (kg * nT * MeV/G) = (1e5)*(kg MeV)
    p2 = 1e-5*p2*mks['MeV'] # kg J = kg**2 m**2 / s**2
    gamma = np.sqrt(1+p2/m0**2/c**2)
    E = (gamma-1)*m0*c**2 # J
    MeV = E/mks['MeV']

    return MeV

def rigidity(MeV,species):
    """
    returns rigidity (p/q) in GV/c (/c is usually dropped)
    """

    gamma,v,m = MeVtogamma(MeV,species);
    # v in m/s
    # m in kg

    R = m*v/abs(mks[species]['q'])*mks['c'] # kg m**2 /s**2 / C = V 
    R = R*1e-9 # R in GV
    return R

def flux2psd(flux,energy,species,energy_unit):
    """
    # psd = flux2psd(flux,energy,species,energy_unit)
    # psd = flux2psd(flux,energy,'e','MeV');
    # species is:
    # electrons: 'electron','e','e-'
    # protons: 'p','H+','p+'
    flux is expected to be a matrix Nt x NE
    flux is in units of #/cm**2/s/sr/(energy_unit)
    energy is expected to be a vector of length NE
    psd is in units of (MeV s)**(-3)
    """

    c_cm = mks['c']*100 # speed of light in cm/s

    species = SelectSpecies(species)
    m0 = mks[species]['m0']
    m0c2 = m0*mks['c']**2/mks['MeV']

    inMeV = EnergyUnitInMeV(energy_unit)

    W = energy*inMeV # energy, MeV
    gamma = W/m0c2+1 # relativistic factor
    p2 = (gamma**2-1)*m0c2**2/c_cm**2 # p**2 = (gamma**2-1)*m0**2*c**2; units of (MeV/cm*s)**2
    psd = flux/inMeV/p2 # psd = flux/p**2
    # #/cm**2/s/sr/MeV / (MeV/cm*s)**2
    # cm**2 / (cm**2 s sr MeV MeV**2 s**2)
    # # / (MeV**3 s**3)

    return psd

def EnergyUnitInMeV(energy_unit):
    """
    inMeV = EnergyUnitInMeV(energy_unit)
    returns the MeV equivalent of 1 energy_unit
    e.g., EnergyUnitInMeV('GeV') = 1000
    """
    if energy_unit == 'eV':
        inMeV = 1e-6
    elif energy_unit == 'keV':
        inMeV = 1e-3
    elif energy_unit == 'MeV':
        inMeV = 1
    elif energy_unit == 'GeV':
        inMeV = 1e3
    else:
        print ('Unknown energy unit, please use either eV, keV, MeV, GeV')
    return inMeV


def GyroRadius(Species,Energy,B):
    """
    # r = GyroRadius(Species,Energy,B)
    # Species - 'e', 'p', etc
    # Energy - particle kinetic energy MeV
    # B - local magnetic field strength, nT
    # r = gyroradius in m
    """

    species = SelectSpecies(Species)
    q = mks[species]['q']
    Bsi = B/1e9 # B in TeSL['a']
    gamma,vmag,m = MeVtogamma(Energy,species) # gamma, speed (m/s), m (kg)
    r = m*vmag/np.abs(q)/Bsi # kg * m/s / C / T = m
    return r

def EBalpha2M(E,B,alpha,species):
    """
    # for E in MeV, alpha in degrees
    # returns M in MeV/G or MeV/nT, using same B unit as input
    """

    mks = constants['mks']
    species = SelectSpecies(species)
    m0 = mks[species]['m0']
    
    c = mks['c'] # m/s
    EJ = E*mks['MeV'] # J = kg (m/s)**2
    p2 = (EJ**2+2*EJ*m0*c**2)/c**2  # (kg m/s)**2
    M = p2*np.sin(alpha*np.pi/180.)**2./(2*m0*B)/mks['MeV'] # J/G -> MeV/G
    return M

def Ealpha2MK(E,alpha,L,species,Bunit = 'G'):
    """
    for E in MeV, alpha in degrees
    Bunit is 'nT' or 'G' (default is 'G')
    returns dipole M,K
    M in MeV/G (or MeV/nT if Bunit='nT')
    K in RE*sqrt(G)(or RE*sqrt(nT) if Bunit='nT')
    """
    
    mks = constants['mks']
    species = SelectSpecies(species)
    m0 = mks[species]['m0']
    K = alphaL2K(alpha,L,Bunit)
    BnT = dipoleB(L,0,0)
    if Bunit =='G':
        Beq = BnT/1e5 # 1 nT = 1E-5 G
    else:
        Beq = BnT

    c = mks['c'] # m/s
    EJ = E*mks['MeV'] # J = kg (m/s)**2
    p2 = (EJ**2+2*EJ*m0*c**2)/c**2  # (kg m/s)**2
    M = p2*np.sin(alpha*np.pi/180.)**2./(2*m0*Beq)/mks['MeV'] # J/G -> MeV/G
    return M,K

def MK2Ealpha(M,K,L,species,Bunit='G'):
    """
    for M in MeV/G and K in RE*sqrt(G)
    (or MeV/nT and RE*sqrt(nT) if Bunit = 'nT')
    returns E in MeV, alpha in degrees
    using dipole field
    Bunit = 'G' by default
    """
    alpha = KL2alpha(K,L,Bunit)
    B = dipoleB(L,0,0) # nT
    if Bunit=='G':
        MG = M
    else:
        MG = M*1e5 # MeV/nT -> MeV/G
        

    E = MBtoMeV(MG,B,alpha,species)
    return E,alpha

def dbydL(f,E,alpha,L,species,dL= 0.001):
    """
    returns df/dL at constant M,K
    f is a function handle with 3 arguments: E,alpha,L
    performs numerical derivative with dL=0.001 or specified by user
    E in MeV, alpha in degrees
    f is evaluated at L and at L+dL, 
    at L+dL, E and alpha are adjusted to preserve M,K
    assumes dipole field for E,alpha <-> M,K conversions
    """

    f1 = f(E,alpha,L)
    M1,K1 = Ealpha2MK(E,alpha,L,species,'G')
    L2 = L+dL
    E2,alpha2 = MK2Ealpha(M1,K1,L2,species,'G')
    f2 = f(E2,alpha2,L2)
    dfdL = (f2-f1)/(L2-L)
    return dfdL
