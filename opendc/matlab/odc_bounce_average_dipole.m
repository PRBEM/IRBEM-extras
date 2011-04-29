function [ba,denom] = bounce_average_dipole(L,MLT,alpha0_deg,local,varargin)
% ba = bounce_average_dipole(L,MLT,alpha0_deg,local,varargin)
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
% [bint,denom] = bounce_average_dipole(...,'no_avg') - returns bounce integral
%  and denominator (in RE), does not average (bav = bint/denom)
%  NOTE: the true bounce integral is bint*v, where v is the velocity
%  and denom*v is the true full bounce period
%
% [...] = bounce_average_dipole(...,'tol',tol) - absolute error tolerance
%  for numerical integral (quadv) (default is 1e-6)
% [...] = bounce_average_dipole(...,'symmetric') - local does not depend
%  on sign of the pitch angle (i.e., northward and southward legs are
%  identical)
% [...] = bounce_average_dipole(...,'hemi-symmetric') - local does not depend
%  on sign of magnetic latitude
%  identical)
% [...] = bounce_average_dipole(...,'Beq',Beq) - provide equatorial field
%    strength at L (otherwise it's 31e3/L^3 (from util.SL.B0)
% [...] = bounce_average_dipole(...,'vectorized',true) - local accepts multiple input points

util = odc_util; % load utility functions and constants

symmetric = false; % user tells us it's symmetric?
hemi_symmetric = false; % user tells us it's invariant to sign of latitude?
vectorized = false;
do_avg = true; % do bounce average? vs integral?
Beq = util.SL.B0/L^3*1e9; % nT
i = 1;
tol = 1e-6; % default
while i <= length(varargin),
    switch(lower(varargin{i})),
        case 'no_avg',
            do_avg = false;
        case 'symmetric',
            symmetric = true;
        case {'hemi-symmetric','hemi_symmetric'},
            hemi_symmetric = true;
        case 'vectorized',
            vectorized = true;
        case 'bm',
            i = i+1;
            Beq = varargin{i};
        case 'tol',
            i = i+1;
            tol = varargin{i};
        otherwise
            error('Unknown argument "%s"',varargin{i});
    end
    i = i+1;
end

if symmetric,
    SIGN_COSPA = 1; % double a single half-bounce
else
    SIGN_COSPA = [-1,1]; % do both half-bounces
end

a0 = alpha0_deg*pi/180; % radians
sina0 = sin(a0);
Bm = Beq./sina0.^2; % mirror field strength
% find mirror lat:
lambdam = util.dipole_mirror_latitude(a0,'rad');

tiny = 1e-6; % small distance to avoid numerical errors at mirror point
lam2 = lambdam*(1-tiny); % northern mirror latitude
if hemi_symmetric, % only do northern hemisphere
    lam1 = 0;
    Tfactor = 1;
else
    lam1 = -lam2;
    Tfactor = 2;
end

T = util.T(sina0);
% <f>_b = 1/T(a0)*int_0^lambdam [f(lam)*cos(lam)*sqrt(1+3*sin^2(lam))/cos(a)]dlam

ba = 0; % numerator for now, normalized after loop
g = 0;
for sign_cospa = SIGN_COSPA, % which half-cycle is this?
    ba = ba+quadv(@(lam)wrapper(local,lam,L,MLT,Beq,Bm,sign_cospa,vectorized),lam1,lam2,tol);
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
    if hemi_symmetric, % dupilcate for souther hemisphere case
        ba = ba*2;
    end
end
denom = 4*L*T;


function out = wrapper(local,maglat,L,MLT,Beq,Bm,sign_cospa,vectorized)
% note: maglat is in rads, MLT in hours
N = length(maglat);
MLTrads = MLT*pi/12;
cmlat = cos(maglat(:));
smlat = sin(maglat(:));
R = L*cmlat.^2;
XYZ = (R.*cmlat)*[cos(MLTrads) sin(MLTrads) 0]; % set X, Y
XYZ(:,end) = R.*smlat;
Blocal = Beq*(1+3*smlat.^2).^(1/2)./cmlat.^6;

maglat_deg = maglat*180/pi;

if (N==1) || vectorized,
    out = local(XYZ,Blocal,Bm,maglat_deg,sign_cospa);
else
    out = local(XYZ,Blocal,Bm,maglat_deg(1),sign_cospa);
    out = repmat(out,length(maglat));
    for i = 2:N,
        out(i,:) = local(XYZ(i,:),Blocal(i),Bm,maglat_deg(i),sign_cospa);
    end
end

cosa0 = sqrt(1-Blocal./Bm);

out = out.*repmat(cmlat.*sqrt(1+3*smlat.^2)./cosa0,1,size(out,2)); % apply dlat vs dt scaling
