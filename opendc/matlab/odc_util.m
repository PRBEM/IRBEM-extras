function util = odc_util
% util = odc_util
% provides a set of constants and handles for utility functions
%
% mirror_lat = util.dipole_mirror_latitude(alpha0)
% compute dipolar mirror latitude of particle with equatorial pitch angle
% alpha0, angles in degrees
% mirror_lat = util.dipole_mirror_latitude(alpha0,'rad')
% angles in radians
%
% maglat = util.BB0toMagLat(BB0) % returns maglat in deg
% maglat = util.BB0toMagLat(BB0,'rad') % returns maglat in radians
% BB0 is B / Bmin
% where B is the local magnetic field strength
% and Bmin is the minimum (i.e., equatorial) magnetic field strength on the
% same field line.%
% maglat is the unsigned dipole latitude
% 
% R = util.rigidity(MeV,species) % returns rigidity (p/q) in GV
%
% [gamma,v,m] = util.MeVtogamma(MeV,species)
% convert Energy in MeV to gamma factor
% and velocity in  m / s, if requested
% and relativistic mass in kg, if requested
% species is 'e' for electrons, 'p' for protons
% v in m/s
% m in kg
%
% T = util.T(y)
% [T,dT] = util.T(y)
% Schulz & Lanzerotti's T(y) = T(sin(alpha0))
% where alpha0 = equatorial pitch angle
% and dT/dy if requested
%
% D = util.D(y)
% Schulz & Lanzerotti's D(y) = D(sin(alpha0))
% where alpha0 = equatorial pitch angle
%
% Q = util.Q(y)
% Schulz & Lanzerotti's Q(y) = Q(sin(alpha0))
% where alpha0 = equatorial pitch angle
%
% Y = util.Y(y)
% Schulz & Lanzerotti's Y(y) = Y(sin(alpha0))
%
% B = util.fce2B(fce)
% computes B in nT given the electron gyrofrequency in Hz
%
% MeV = util.MBtoMeV(M,B,alpha,species)
% given M in MeV/G, and B in nT, alpha in degrees, and species
% returns energy, in MeV
%
% Tg = util.GyroPeriod(Species,Energy,MagLat,L)
%   calculates particle gyro period Tg (seconds), dipole
% Tb = util.BouncePeriod(Species,Energy,PitchAngle,L)
%   calculates particle bounce period Tb (seconds), dipole
% Td = util.DriftPeriod(Species,Energy,PitchAngle,L);
%   calculates particle drift period (seconds), dipole
%   Species: 'e' for electrons, 'p' for protons
%   Energy: in MeV
%   MagLat: in degrees
%   PitchAngle: in degrees (equatorial pitch angle)
%   L: dimensionless dipole L value
%
% r = util.GyroRadius(Species,Energy,B)
%   Species - 'e', 'p', etc
%   Energy - particle kinetic energy MeV
%   B - local magnetic field strength, nT
%   r = gyroradius in m
%
% B = util.dipoleB(L,MagLat,phi_deg)
% [B,Bvec] = util.dipoleB(L,MagLat,phi_deg)
% [B,Bvec,XYZ] = util.dipoleB(L,MagLat,phi_deg)
% [B,Bx,By,Bz,X,Y,Z] = util.dipoleB(L,MagLat,phi_deg)
% computes magnitude of dipole field, nT
% returns components and positions (in RE) if requested
% MagLat: in degrees
% L: dimensionless dipole L value
% phi_deg: azimuth angle, degrees
% Bvec = [Bx(:) By(:) Bz(:)]
% XYZ = [X(:) Y(:) Z(:)]
%
% I = util.dipoleIJ(L,y);
% [I,J] = util.dipoleIJ(L,y,MeV,species)
% computes I = L*Y(y) (in RE)
% and, if requested J = 2*p*I (in RE*MeV/c)
% for species ('e' or 'p')
% energy in MeV
% L: dimensionless dipole L value
% and y = sin(alpha_equatorial)
%
% K = alphaL2K(alpha,L,Bunit)
% returns K in RE*sqrt(nT) for given equatorial pitch angle
% alpha, in deg, and L shell, in RE
% Bunit - optional, if 'G' assumes K in RE*sqrt(G). Otherwise K in
%   RE*sqrt(nT)
%
% alpha = KL2alpha(K,L,Bunit)
% returns equatorial pitch angle alpha in degrees
% for given K, L. K in RE*sqrt(nT)
% Bunit - optional, if 'G' returns K in RE*sqrt(G). Otherwise K in
%   RE*sqrt(nT)
%
% M = EBalpha2M(E,B,alpha,species)
% for E in MeV, alpha in degrees
% species is 'e' for electrons, 'p' for protons
% returns M in MeV/G or MeV/nT, using same B unit as input
%
% [M,K] = Ealpha2MK(E,alpha,L,species,Bunit)
% for E in MeV, alpha in degrees
% species is 'e' for electrons, 'p' for protons
% Bunit is 'nT' or 'G' (default is 'G')
% returns dipole M,K
% M in MeV/G (or MeV/nT if Bunit='nT')
% K in RE*sqrt(G)(or RE*sqrt(nT) if Bunit='nT')
%
% [E,alpha] = MK2Ealpha(M,K,L,species,Bunit)
% for M in MeV/G and K in RE*sqrt(G)
%  (or MeV/nT and RE*sqrt(nT) if Bunit = 'nT')
% returns E in MeV, alpha in degrees
% using dipole field
% species is 'e' for electrons, 'p' for protons
% Bunit = 'G' by default
%
% dfdL = dbydL(f,E,alpha,L,species)
% dfdL = dbydL(f,E,alpha,L,species,dL)
% returns df/dL at constant M,K
% f is a function handle with 3 arguments: E,alpha,L
% performs numerical derivative with dL=0.001 or specified by user
% E in MeV, alpha in degrees
% assumes dipole field
% species is 'e' for electrons, 'p' for protons
%
% psd = flux2psd(flux,energy,species,energy_unit)
% psd = flux2psd(flux,energy,'e','MeV');
% species is:
% electrons: 'electron','e','e-'
% protons: 'p','H+','p+'
% flux is expected to be a matrix Nt x NE
% flux is in units of #/cm^2/s/sr/(energy_unit)
% energy is expected to be a vector of length NE
% psd is in units of (MeV s)^(-3)
%
% inMeV = EnergyUnitInMeV(energy_unit)
% returns the MeV equivalent of 1 energy_unit
% e.g., EnergyUnitInMeV('GeV') = 1000
%
% SI/mks constants
% util.mks:
% mks.e = 1.602176487e-19; % Coulombs per fundamental charge
% mks.c = 299792458; % m/s
% mks.eV = mks.e; % Joules
% mks.keV = mks.eV*1e3; % Joules
% mks.MeV = mks.eV*1e6; % Joules
% mks.GeV = mks.eV*1e9; % Joules
% mks.epsilon0 = 8.854187817e-12; % F/m = C^2 s^2  / kg / m^3 - permitivity of free space
% mks.electron.q = -mks.e; % charge, Coulombs
% mks.electron.m0 = 9.10938215e-31; % rest mass, kg
% mks.proton.m0 = 1.672621637e-27; % rest mass, kg
% mks.proton.q = +mks.e; % charge, Coulombs
% mks.R_E =  6371.2e3 % Earth Radius, m (IAU 1966 mean?, used by IRBEM lib GDZ)
% (mks.RE deprecated because it was in km which is not an mks unit)

%
% Schulz & Lanzerotti constants
% util.SL:
% SL.T0 = T(y=0)
% SL.T1 = T(y=1)
% SL.B0 = 3.1e-5 T, dipole moment of the Earth (0.31 G) (S&L value)
% SL.a = 6371.2e3; % Earth Radius, (S&L value)
% SL.Q0 = Q(y=0)
% SL.Q1 = Q(y=1)
% SL.Qp1 = Q'(y=1)
% SL.Y0 = Y(y=0)

global odc_constants
if isempty(odc_constants),
    
    % SI/mks constants
    mks.e = 1.602176487e-19; % Coulombs per fundamental charge
    mks.c = 299792458; % m/s
    mks.eV = mks.e; % Joules
    mks.keV = mks.eV*1e3; % Joules
    mks.MeV = mks.eV*1e6; % Joules
    mks.GeV = mks.eV*1e9; % Joules
    mks.epsilon0 = 8.854187817e-12; % F/m = C^2 s^2  / kg / m^3 - permitivity of free space
    mks.mu0 = 4*pi*1e-7; % H/m = N/A^2 = kg m/C^2 - permeability of free space
    mks.R_E =  6371.2e3; % m, IAU 1966? IBEM GDZ value
    
    
    mks.electron.q = -mks.e; % charge
    mks.electron.m0 = 9.10938215e-31; % rest mass, kg
    
    mks.proton.m0 = 1.672621637e-27; % rest mass, kg
    mks.proton.q = +mks.e; % charge
    
    % Schulz & Lanzerotti constants
    SL.T0 = 1+1/(2*sqrt(3))*log(2+sqrt(3)); % S&L 1.28a
    SL.T1 = pi/6*sqrt(2); % S&L 1.28b
    SL.B0 = 31e3*1e-9; % T dipole field at lambda=0,L=1 (S&L value)
    SL.a = 6371.2e3; % Earth Radius, meters (S&L value)
    SL.Q0 = -27.1266694; % S&L 1.78c
    SL.Q1 = -90*SL.T1; % S&L 1.78b
    SL.Qp1 = (15/2)*(9*SL.T0 - 41*SL.T1); % S&L 1.78d
    SL.Y0 = 2*SL.T0; % S&L 1.31 (limiting case in subsequent text)
    
    odc_constants = struct('mks',mks,'SL',SL);
end

util = odc_constants;
util.dipole_mirror_latitude = @dipole_mirror_latitude;
util.BB0toMagLat = @BB0toMagLat;
util.MeVtogamma = @MeVtogamma;
util.T = @SL_T;
util.D = @SL_D;
util.Q = @SL_Q;
util.Y = @SL_Y;
util.fce2B = @fce2B;
util.GyroPeriod = @GyroPeriod;
util.BouncePeriod = @BouncePeriod;
util.DriftPeriod = @DriftPeriod;
util.GyroRadius = @GyroRadius;
util.dipoleB = @dipoleB;
util.dipoleIJ = @dipoleIJ;
util.MBtoMeV = @MBtoMeV;
util.rigidity = @rigidity;
util.alphaL2K = @alphaL2K;
util.KL2alpha = @KL2alpha;
util.flux2psd = @flux2psd;
util.EnergyUnitInMeV = @EnergyUnitInMeV;
util.EBalpha2M = @EBalpha2M;
util.Ealpha2MK = @Ealpha2MK;
util.MK2Ealpha = @MK2Ealpha;
util.dbydL = @dbydL;

function K = alphaL2K(alpha,L,Bunit)
% returns K in RE*sqrt(nT) for given equatorial pitch angle
% alpha, in deg, and L shell, in RE
% Bunit - optional, if 'G' assumes K in RE*sqrt(G). Otherwise K in
%   RE*sqrt(nT)
global odc_constants

if nargin < 3,
    Bunit = 'nT';
end

y = sind(alpha);
Bm = odc_constants.SL.B0*1e9./L.^3./y.^2; % mirror field, nT
if any(Bunit=='G'),
    Bm = Bm/1e5; % nT to G
end
K = SL_Y(y).*L.*sqrt(Bm);

function alpha = KL2alpha(K,L,Bunit)
% returns equatorial pitch angle alpha in degrees
% for given K, L. K in RE*sqrt(nT)
% Bunit - optional, if 'G' returns K in RE*sqrt(G). Otherwise K in
%   RE*sqrt(nT)

global odc_constants
persistent table
if isempty(table),
    dy = 1e-4;
    table.y = dy:dy:1;
    table.Yy = SL_Y(table.y)./table.y;
end

if nargin < 3,
    Bunit = 'nT';
end

if any(Bunit=='G'),
    K = K*sqrt(1e5); % RE*sqrt(G) to RE*sqrt(nT)
end

% K*sqrt(L)/sqrt(B0) = Y(y)/y
KLB = K.*sqrt(L)./sqrt(odc_constants.SL.B0*1e9);
y = interp1(table.Yy,table.y,KLB,'linear');
alpha = asind(y);

function B = fce2B(fce)
% returns B in nT for electron gyro given in Hz
global odc_constants
electron = odc_constants.mks.electron; % shorthand
B = fce*2*pi*electron.m0/abs(electron.q)*1e9; % nT

function [T,dT] = SL_T(y)
% Schulz & Lanzerotti's T(y)
% and dT = dT/dy if requested

global odc_constants
SL = odc_constants.SL;
T = SL.T0-0.5*(SL.T0-SL.T1)*(y+sqrt(y)); % 1/4 bounce integral of 1
if nargout >= 2,
    dT = -0.5*(SL.T0-SL.T1)*(1+1./sqrt(y)/2); % dT/dy
end

function Q = SL_Q(y)
% computes Q(y) from S&L, 1.79
global odc_constants
SL = odc_constants.SL;
Q = SL.Q0+(2*SL.Q1 - 2*SL.Q0 - (1/4)*SL.Qp1)*y.^4 + (SL.Q0-SL.Q1+(1/4)*SL.Qp1)*y.^8; % Q(y) from S&L, 1.79

function D = SL_D(y)
% computes D(y) from S&L, 1.36

D = SL_T(y)/2-SL_Y(y)/12;

function Y = SL_Y(y)
% computes Y(y) from S&L 1.31
global odc_constants
SL = odc_constants.SL;
Y = 2*(1-y)*SL.T0 + (SL.T0-SL.T1)*(y.*log(y) + 2*y-2*y.^0.5); % S&L 1.31
Y(y==0) = SL.Y0;

function mirror_lat = dipole_mirror_latitude(alpha0,units)
% mirror_lat = dipole_mirror_latitude(alpha0)
% compute dipolar mirror latitude of particle with equatorial pitch angle
% alpha0, angles in degrees
% mirror_lat = dipole_mirror_latitude(alpha0,'rad')
% angles in radians

if nargin < 2,
    units = [];
end

if isempty(units),
    units = 'deg';
end

if lower(units(1))=='d',
    torad = pi/180;
else
    torad = 1;
end

sina0 = sin(alpha0*torad);
if sina0==0,
    mirror_lat = inf;
    return;
end
% note, below fixes error in Shprits thesis, eqn  F13
% which has sina0^2. Instead use sina0^4 from Shprits 2006 eqn 10.
mirror_lat = nan(size(sina0));
for i = 1:numel(sina0),
    Px = [1, 0, 0, 0, 0, 3*sina0(i)^4, -4*sina0(i)^4];
    xroots = roots(Px);
    xroot = xroots((imag(xroots)==0) & (real(xroots)>0));
    mirror_lat(i) = acos(sqrt(xroot))/torad; % mirror latitude
end
function maglat = BB0toMagLat(BB0,units)
% maglat = BB0toMagLat(BB0) % returns maglat in deg
% maglat = BB0toMagLat(BB0,'rad') % returns maglat in radians
%
% BB0 is B / Bmin
% where B is the local magnetic field strength
% and Bmin is the minimum (i.e., equatorial) magnetic field strength on the
% same field line.
%
% maglat is the unsigned dipole latitude

% we invert the B/Bmin vs maglat equation via a look-up table
% we use a persistent variable so we don't have to regenerate the table
% for each call to this function
persistent maglat_table
if isempty(maglat_table),
    maglat_table.deg = (0:0.001:90)';
    maglat_table.rads = maglat_table.deg*pi/180;
    % this next bit is the expression for B/Bmin vs magnetic latitude for a dipole
    maglat_table.bb0 = (1+3*sin(maglat_table.rads).^2).^(1/2)./cos(maglat_table.rads).^6;
end

if nargin < 2,
    units = 'deg';
end
switch(lower(units(1))),
    case 'd',
        maglat = interp1(maglat_table.bb0,maglat_table.deg,BB0,'linear');
    case 'r',
        maglat = interp1(maglat_table.bb0,maglat_table.rads,BB0,'linear');
    otherwise
        error('Unknown units "%s"',units);
end

function Td = DriftPeriod(Species,Energy,PitchAngle,L)
%   calculates particle drift period (seconds), dipole
%   Species: 'e' for electrons, 'p' for protons
%   Energy: in MeV
%   PitchAngle: in degrees (equatorial pitch angle)
%   L: dimensionless dipole L value
%

global odc_constants

species = SelectSpecies(Species);
m0 = odc_constants.mks.(species).m0;

q = odc_constants.mks.(species).q; % C
c = odc_constants.mks.c; % m/s
a = odc_constants.SL.a; % Earth Radius, meters
B0 = odc_constants.SL.B0; % T

gamma = MeVtogamma(Energy,Species);

y = sind(PitchAngle);
Ty = SL_T(y);
Dy = SL_D(y);

% S&L eq 1.35, w/o minus sign, assume always want positive drift velocity
f = (3*L/2/pi./gamma).*(gamma.^2-1)*(c./a).^2.*(m0*c./abs(q)./B0).*(Dy./Ty)/c; % extra 1/c for SI
Td = 1./f; % seconds

function Tb = BouncePeriod(Species,Energy,PitchAngle,L)
% function Tb = BouncePeriod(Species,Energy,PitchAngle,L)
% calculates particle bounce period Tb (seconds), dipole
% Species: 'e' for electrons, 'p' for protons
% Energy: in MeV
% PitchAngle: in degrees (equatorial pitch angle)
% L: dimensionless dipole L value

global odc_constants

a = odc_constants.SL.a; % Earth Radius, meters
y = sind(PitchAngle);
[gamma,v] = MeVtogamma(Energy,Species);

Ty = SL_T(y);
Tb = 4*L.*a./v.*Ty;

function Tg = GyroPeriod(Species,Energy,MagLat,L)
% function Tg = GyroPeriod(Species,Energy,MagLat,L)
% calculates particle gyro period Tg (seconds), dipole
% Species: 'e' for electrons, 'p' for protons
% Energy: in MeV
% MagLat: in degrees
% L: dimensionless dipole L value
%
% Calculation done in .* so that vector and matrix input can work,
% if all inputs are the same size (or scalar)

global odc_constants

species = SelectSpecies(Species);
q = odc_constants.mks.(species).q; % C

[gamma,v,m] = MeVtogamma(Energy,Species);

B = dipoleB(L,MagLat)/1e9; % T

f = abs(q)*B./(2*pi*m); % no "c" in denominator in SI units
Tg = 1./f;

function [B,Bx,By,Bz,X,Y,Z] = dipoleB(L,MagLat,phi_deg)
% B = dipoleB(L,MagLat,phi_deg)
% [B,Bvec] = util.dipoleB(L,MagLat,phi_deg)
% [B,Bvec,XYZ] = util.dipoleB(L,MagLat,phi_deg)
% [B,Bx,By,Bz,X,Y,Z] = util.dipoleB(L,MagLat,phi_deg)
% computes magnitude of dipole field, nT
% returns components and positions (in RE) if requested
% MagLat: in degrees
% L: dimensionless dipole L value
% phi_deg: azimuth angle, degrees
% Bvec = [Bx(:) By(:) Bz(:)]
% XYZ = [X(:) Y(:) Z(:)]

global odc_constants
Beq = odc_constants.SL.B0./L.^3*1e9; % T to nT
smlat = sind(MagLat);
cmlat = cosd(MagLat);
cmlat6 = cmlat.^6;
B = Beq.*sqrt(1+3*smlat.^2)./cmlat6;
if nargout >= 2, % Bvec or Bx, By, Bz
    cphi = cosd(phi_deg);
    sphi = sind(phi_deg);
    
    % angular part of Bx, By, Bz
    Btmp = Beq./cmlat6;
    Bx = -3*cphi.*cmlat.*smlat.*Btmp;
    By = -3*sphi.*cmlat.*smlat.*Btmp;
    Bz = -(3*smlat.^2 - 1).*Btmp;
    
    if nargout == 2, % [B,Bvec] style output
        Bx = [Bx(:), By(:), Bz(:)]; % Bvec output
    end
    
    if ismember(nargout,[3,7]), % compute XYZ or X Y Z
        
        R = L.*cmlat.^2;
        X = R.*cmlat.*cphi;
        Y = R.*cmlat.*sphi;
        Z = R.*smlat;
        
        if nargout==3, % [B, Bvec, XYZ] XYZ style output
            By = [X(:), Y(:), Z(:)]; % XYZ output
        end
    end
end



function species = SelectSpecies(Species)
% convert various alternatives to standard species name

switch(lower(Species)),
    case {'e','e-','electron'},
        species = 'electron';
    case {'p','p+','h+','proton'},
        species = 'proton';
    otherwise
        error('%s:Unknown Species [%s]',mfilename,Species);
end

function [gamma,v,m] = MeVtogamma(MeV,species)
% compute gamma, v, m
% given energy in MeV and species 'e','p', etc
% v in m/s
% m in kg

global odc_constants
mks = odc_constants.mks;

species = SelectSpecies(species);

m0 = mks.(species).m0;

W = MeV*mks.MeV; % Energy, Joules

gamma = 1+W/(m0*mks.c^2); % relativistic factor
if nargout >= 2,    
    v=mks.c.*sqrt(1-gamma.^-2);
end
if nargout >= 3,
    m = m0*gamma;
end

function [I,J] = dipoleIJ(L,y,MeV,species)
% I = util.dipoleIJ(L,y);
% [I,J] = util.dipoleIJ(L,y,MeV,species)
% computes I = L*Y(y) (in RE)
% and, if requested J = 2*p*I (in RE*MeV/c)
% for species ('e' or 'p')
% energy in MeV
% L: dimensionless dipole L value
% and y = sin(alpha_equatorial)

global odc_constants
I = L.*SL_Y(y);

if nargout > 1,
    [gamma,v,m] = MeVtogamma(MeV,species);
    p = m.*v/odc_constants.mks.MeV*odc_constants.mks.c; % MeV/c
    J = 2.*p.*I; % RE*MeV/c
end


function MeV = MBtoMeV(M,B,alpha,species)
% given M in MeV/G, and B in nT, alpha in degrees, and species
% returns energy, in MeV

global odc_constants
mks = odc_constants.mks;

Bm = B./sind(alpha).^2; % Bmirror
species = SelectSpecies(species);
m0 = mks.(species).m0;

% M = p^2/2/m0/Bm
p2 = 2*m0*Bm.*M; % (kg * nT * MeV/G) = (1e5)*(kg MeV)
p2 = 1e-5*p2*mks.MeV; % kg J = kg^2 m^2 / s^2
gamma = sqrt(1+p2/m0^2/mks.c^2);
E = (gamma-1)*m0*mks.c^2; % J
MeV = E/mks.MeV;

function R = rigidity(MeV,species) 
% returns rigidity (p/q) in GV/c (/c is usually dropped)
global odc_constants
mks = odc_constants.mks;

[~,v,m] = MeVtogamma(MeV,species);
% v in m/s
% m in kg

R = m.*v/mks.e*mks.c; % kg m^2 /s^2 / C = V 
R = R*1e-9; % R in GV


function psd = flux2psd(flux,energy,species,energy_unit)
% psd = flux2psd(flux,energy,species,energy_unit)
% psd = flux2psd(flux,energy,'e','MeV');
% species is:
% electrons: 'electron','e','e-'
% protons: 'p','H+','p+'
% flux is expected to be a matrix Nt x NE
% flux is in units of #/cm^2/s/sr/(energy_unit)
% energy is expected to be a vector of length NE
% psd is in units of (MeV s)^(-3)

global odc_constants

c_cm = odc_constants.mks.c*100; % speed of light in cm/s

species = SelectSpecies(species);
m0c2 = odc_constants.mks.(species).m0*odc_constants.mks.c^2/odc_constants.mks.MeV;

inMeV = EnergyUnitInMeV(energy_unit);

W = energy*inMeV; % energy, MeV
gamma = W/m0c2+1; % relativistic factor
p2 = (gamma.^2-1)*m0c2^2/c_cm^2; % p^2 = (gamma^2-1)*m0^2*c^2; units of (MeV/cm*s)^2
psd = flux./inMeV./p2; % psd = flux/p^2
% #/cm^2/s/sr/MeV / (MeV/cm*s)^2
% cm^2 / (cm^2 s sr MeV MeV^2 s^2)
% # / (MeV^3 s^3)



function inMeV = EnergyUnitInMeV(energy_unit)
% inMeV = EnergyUnitInMeV(energy_unit)
% returns the MeV equivalent of 1 energy_unit
% e.g., EnergyUnitInMeV('GeV') = 1000
switch(energy_unit),
    case 'eV',
        inMeV = 1e-6;
    case 'keV',
        inMeV = 1e-3;
    case 'MeV',
        inMeV = 1;
    case 'GeV',
        inMeV = 1e3;
    otherwise
        error('Unknown energy unit %s',energy_unit);
end

function r = GyroRadius(Species,Energy,B)
% r = GyroRadius(Species,Energy,B)
% Species - 'e', 'p', etc
% Energy - particle kinetic energy MeV
% B - local magnetic field strength, nT
% r = gyroradius in m
global odc_constants
species = SelectSpecies(Species);
Bsi = B/1e9; % B in Tesla
[gamma,vmag,m] = MeVtogamma(Energy,species); % gamma, speed (m/s), m (kg)
q = abs(odc_constants.mks.(species).q); % charge, C
r = m.*vmag./q./Bsi; % kg * m/s / C / T = m

function M = EBalpha2M(E,B,alpha,species)
% for E in MeV, alpha in degrees
% returns M in MeV/G or MeV/nT, using same B unit as input
global odc_constants
mks = odc_constants.mks;
sp = SelectSpecies(species);
m0 = mks.(sp).m0; % kg
c = mks.c; % m/s
EJ = E*mks.MeV; % J = kg (m/s)^2
p2 = (EJ.^2+2*EJ*m0*c^2)/c^2;  % (kg m/s)^2
M = p2.*sind(alpha).^2./(2*m0*B)/mks.MeV; % J/G -> MeV/G

function [M,K] = Ealpha2MK(E,alpha,L,species,Bunit)
% for E in MeV, alpha in degrees
% Bunit is 'nT' or 'G' (default is 'G')
% returns dipole M,K
% M in MeV/G (or MeV/nT if Bunit='nT')
% K in RE*sqrt(G)(or RE*sqrt(nT) if Bunit='nT')
global odc_constants
mks = odc_constants.mks;
if nargin < 5,
    Bunit = 'G';
end
sp = SelectSpecies(species);
K = alphaL2K(alpha,L,Bunit);
BnT = dipoleB(L,0,0);
if isequal(Bunit,'G'),
    Beq = BnT/1e5; % 1 nT = 1E-5 G
else
    Beq = BnT;
end
m0 = mks.(sp).m0; % kg
c = mks.c; % m/s
EJ = E*mks.MeV; % J = kg (m/s)^2
p2 = (EJ.^2+2*EJ*m0*c^2)/c^2;  % (kg m/s)^2
M = p2.*sind(alpha).^2./(2*m0*Beq)/mks.MeV; % J/G -> MeV/G


function [E,alpha] = MK2Ealpha(M,K,L,species,Bunit)
% for M in MeV/G and K in RE*sqrt(G)
%  (or MeV/nT and RE*sqrt(nT) if Bunit = 'nT')
% returns E in MeV, alpha in degrees
% using dipole field
% Bunit = 'G' by default
if nargin < 5,
    Bunit = 'G';
end

alpha = KL2alpha(K,L,Bunit);
B = dipoleB(L,0,0); % nT
if isequal(Bunit,'G'),
    MG = M;
else
    MG = M*1e5; % MeV/nT -> MeV/G
end

E = MBtoMeV(MG,B,alpha,species);

function dfdL = dbydL(f,E,alpha,L,species,dL)
% returns df/dL at constant M,K
% f is a function handle with 3 arguments: E,alpha,L
% performs numerical derivative with dL=0.001 or specified by user
% E in MeV, alpha in degrees
% f is evaluated at L and at L+dL, 
%  at L+dL, E and alpha are adjusted to preserve M,K
% assumes dipole field for E,alpha <-> M,K conversions
if nargin < 6,
    dL = 0.001;
end

f1 = f(E,alpha,L);
[M1,K1] = Ealpha2MK(E,alpha,L,species,'G');
L2 = L+dL;
[E2,alpha2] = MK2Ealpha(M1,K1,L2,species,'G');
f2 = f(E2,alpha2,L2);
dfdL = (f2-f1)/(L2-L);


