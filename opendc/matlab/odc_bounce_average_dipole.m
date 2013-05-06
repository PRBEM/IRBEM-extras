function [ba,denom] = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,varargin)
% ba = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,varargin)
% bounce average "local" over dipole field line
% L - L value of field line
% MLT - Magnetic local time of field line
% alpha0_deg - equatorial pitch angle, degrees
% local - what to average
%  @function(XYZ,Blocal,Bm,maglat,sign_cospa) - evaluated at each point on field line
%    XYZ - RE in dipole-local time coordinates: z is along dipole, x is along MLT=0, y is along MLT=0600
%    Bm is mirror field strength in nT
%    maglat is signed magnetic latitude, degrees
%    can return scalar or 1 x M vector
%    sign_cospa is sign of cos(local pitch angle) (will do for -1 and 1 to get
%    full bounce orbit)
% options:
% [bint,denom] = odc_bounce_average_dipole(...,'no_avg') - returns bounce integral
%  and denominator (in RE), does not average (bav = bint/denom)
%  NOTE: the true bounce integral is bint*v, where v is the velocity
%  and denom*v is the true full bounce period
%
% [...] = odc_bounce_average_dipole(...,'method',method) - specifies
%  integration method: 'quad', 'quadl', 'quadgk', 'quadv', 'trace'
%  method defaults to 'quad', unless "vectorized" in which case
%  it defaults to quadv. Note that quad and quadv aren't good
%  for functions (like diffusion coefficients) that are zero over
%  large parts of the field line. Use 'trace' instead
% [...] = odc_bounce_average_dipole(...,'tol',tol) - specifies absolute
%  tolerance for numerical integrals (quad, quadl, quadv, quadgk), default
%  is 1e-6. Ignored for method "trace"
% [...] = odc_bounce_average_dipole(...,'reltol',reltol) - specifies relative
%  tolerance for numerical integrals by quadgk only default 1e-6.
% [...] = odc_bounce_average_dipole(...,'Nlats',Nlats) - specifies number of
% latitudes to use in trace integral (method='trace' only). Default 201 (an
% odd number is probably better)
% [...] = odc_bounce_average_dipole(...,'symmetric') - local does not depend
%  on sign of the pitch angle (i.e., northward and southward legs are
%  identical)
% [...] = odc_bounce_average_dipole(...,'hemi-symmetric') - local does not depend
%  on sign of magnetic latitude
%  identical)
% [...] = odc_bounce_average_dipole(...,'Beq',Beq) - provide equatorial field
%    strength at L (otherwise it's 31e3/L^3 (from util.SL.B0)
% [...] = odc_bounce_average_dipole(...,'vectorized',true) - local accepts multiple input points
%
% Implements psi transform from Orlova and Shprits, 2011, Phys. Plasmas

% O&S denotes equation or info from Orlova and Shprits, 2011

util = odc_util; % load utility functions and constants

method = '';
symmetric = false; % user tells us it's symmetric?
hemi_symmetric = false; % user tells us it's invariant to sign of latitude?
vectorized = false;
do_avg = true; % do bounce average? vs integral?
Beq = util.SL.B0/L^3*1e9; % nT
i = 1;
Nlats = 201;
tol = 1e-6; % absolute tolerance for numerical integrals
reltol = 1e-6; % relative tolerance for numerical integrals
while i <= length(varargin),
    switch(lower(varargin{i})),
        case 'no_avg',
            do_avg = false;
        case 'method',
            i = i+1;
            method = varargin{i};
        case 'symmetric',
            symmetric = true;
        case {'hemi-symmetric','hemi_symmetric'},
            hemi_symmetric = true;
        case 'vectorized',
            vectorized = true;
        case 'bm',
            i = i+1;
            Beq = varargin{i};
        case 'nlats',
            i = i+1;
            Nlats = varargin{i};
        case 'tol',
            i = i+1;
            tol = varargin{i};
        case 'reltol',
            i = i+1;
            reltol = varargin{i};
        otherwise
            error('Unknown argument "%s"',varargin{i});
    end
    i = i+1;
end

a0 = alpha0_deg*pi/180; % radians
sina0 = sin(a0);
Bm = Beq./sina0.^2; % mirror field strength
% find mirror lat:
lambdam = util.dipole_mirror_latitude(a0,'rad');

lam2 = lambdam; % northern mirror latitude
if hemi_symmetric, % only do northern hemisphere
    lam1 = 0;
    Tfactor = 1;
    psi1 = 0; % southern mirror point
    psi2 = pi/2; % equator
else
    lam1 = -lam2;
    Tfactor = 2;
    psi1 = 0; % southern mirror point
    psi2 = pi; % norhtern mirror point
end

if isempty(method),
    method = 'quadv';
end

if strcmpi(method,'trace'),
    options = {};
    if ~do_avg,
        options{end+1} = 'no_avg';
    end
    if symmetric,
        options{end+1} = 'symmetric';
    end
    
    maglat_deg = linspace(lam1*180/pi,lam2*180/pi,Nlats)';
    [Blocal,Bvec,XYZ] = util.dipoleB(L,maglat_deg,MLT*15);
    Blocal = Blocal*Beq/(util.SL.B0*1e9/L.^3);
    [ba,g] = odc_bounce_average_trace(XYZ,Blocal,local,sign(maglat_deg),options{:});
    
    if ~do_avg && hemi_symmetric, % add in southern half of field line
        ba = ba*2;
        g = g*2;
    end
    
    denom = g;
    
    return
end

switch(lower(method)),
    case 'quad',
        ifunc = @(func,x1,x2,tol,reltol)quad(func,x1,x2,tol);
    case 'quadl',
        ifunc = @(func,x1,x2,tol,reltol)quadl(func,x1,x2,tol);
    case 'quadv',
        ifunc = @(func,x1,x2,tol,reltol)quadv(func,x1,x2,tol);
    case 'quadgk',
        ifunc = @(func,x1,x2,tol,reltol)quadgk(func,x1,x2,'AbsTol',tol,'RelTol',reltol);
    otherwise
        error('Unrecognized integration method "%s"',method);
end

if symmetric,
    SIGN_COSPA = 1; % double a single half-bounce
else
    SIGN_COSPA = [-1,1]; % do both half-bounces
end

T = util.T(sina0);
% <f>_b = 1/T(a0)*int_0^lambdam [f(lam)*cos(lam)*sqrt(1+3*sin^2(lam))/cos(a)]dlam

ba = 0; % numerator for now, normalized after loop
g = 0;
for sign_cospa = SIGN_COSPA, % which half-cycle is this?
    ba = ba+ifunc(@(psi)wrapper(local,psi,L,MLT,Beq,Bm,sign_cospa,vectorized),psi1,psi2,tol,reltol);
    g = g+T*Tfactor;
end

% now normalize: to bounce average
% or to bounce int (without 1/v term)

if do_avg,
    ba = ba/g; % compute average
else
    if symmetric, % dupilcate for sing_cospa=-1 case
        ba = ba*2;
    end
    if hemi_symmetric, % dupilcate for southern hemisphere case
        ba = ba*2;
    end
end
denom = 4*L*T;


function out = wrapper(local,psi,L,MLT,Beq,Bm,sign_cospa,vectorized)
% note: psi is in rads, MLT in hours

util = odc_util;

if size(psi,2)==numel(psi),
    psi = psi'; % make column vector
    transpose = true;
else
    transpose = false;
end

N = length(psi);
a0_deg = asind(sqrt(Beq./Bm));
BmoverB0 = Bm./Beq;
BoverBm = 1-cosd(a0_deg).^2.*sin(psi).^2; % eq 3 from O&S
BoverB0 = max(1,BoverBm*BmoverB0); % sometimes these formulae give BoverB0<1 by a tiny bit at lambda=0
maglat = util.BB0toMagLat(BoverB0,'rad');
maglat = -maglat.*sign(cos(psi));
maglat_deg = maglat*180/pi;

[Blocal,Bvec,XYZ] = util.dipoleB(L,maglat_deg,MLT*15); % need Blocal, XYZ
Blocal = Beq*BoverB0; % replace, ignore Bvec

y = sqrt(Beq./Bm); % sin(a0)
smlat = sin(maglat);
cmlat = cos(maglat);
factor = abs(cos(psi)./smlat); % positive definite
factor(maglat==0) = 3/sqrt(2)*tand(a0_deg); % O&S deal with singularity, after eq 14
F = 2/3*sqrt(1-y.^2)./y.^2.*(1+3*smlat.^2).*cmlat.^8./(3+5.*smlat.^2).*factor; % O&S 13

if (N==1) || vectorized,
    out = local(XYZ,Blocal,Bm,maglat_deg,sign_cospa).*F;
else
    out = local(XYZ,Blocal,Bm,maglat_deg(1),sign_cospa).*F(1);
    out = repmat(out,size(maglat));
    
    for i = 2:N,
        out(i,:) = local(XYZ(i,:),Blocal(i),Bm,maglat_deg(i),sign_cospa).*F(i);
    end
end

if transpose,
    out = out';
end
