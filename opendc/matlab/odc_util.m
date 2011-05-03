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
% [gamma,v,m] = util.MeVtogamma(MeV,species)
% convert Energy in MeV to gamma factor
% and velocity in  m / s, if requested
% and relativistic mass in kg, if requested
% species is 'e' for electrons, 'p' for protons
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
% mks.RE =  6378.136e3; % Earth Radius, m

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
    mks.RE =  6378.136e3; % Earth Radius, m
    
    
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
util.dipoleB = @dipoleB;

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
% note, below fixes error in Shprits thesis, eqn  F13
% which has sina0^2. Instead use sina0^4 from Shprits 2006 eqn 10.
Px = [1, 0, 0, 0, 0, 3*sina0^4, -4*sina0^4];
xroots = roots(Px);
xroot = xroots((imag(xroots)==0) & (xroots>0));
mirror_lat = acos(sqrt(xroot))/torad; % mirror latitude

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
    maglat_table.deg = (0:0.1:90)';
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

global odc_constants
mks = odc_constants.mks;

species = SelectSpecies(species);

m0 = mks.(species).m0;

W = MeV*mks.MeV; % Energy, Joules

gamma = 1+W/(m0*mks.c^2); % relativistic factor
if nargout >= 2,
    v = sqrt((gamma.^2-1))*mks.c./gamma; % relativistic velocity
end
if nargout >= 3,
    m = m0*gamma;
end

