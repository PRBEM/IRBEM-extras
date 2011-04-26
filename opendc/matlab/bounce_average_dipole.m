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
% [...] = bounce_average_dipole(...,'symmetric') - local does not depend
%  on sign of the pitch angle (i.e., northward and southward legs are
%  identical)
% [...] = bounce_average_dipole(...,'hemi-symmetric') - local does not depend
%  on sign of magnetic latitude
%  identical)
% [...] = bounce_average_dipole(...,'Beq',Beq) - provide equatorial field
%    strength at L (otherwise it's 30e3/L^3)
% [...] = bounce_average_dipole(...,'vectorized',true) - local accepts multiple input points

symmetric = false; % user tells us it's symmetric?
hemi_symmetric = false; % user tells us it's invariant to sign of latitude?
vectorized = false;
do_avg = true; % do bounce average? vs integral?
Beq = 30e3/L^3;
i = 1;
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
lambdam = dipole_mirror_latitude(a0,'rad');

tiny = 1e-6; % small distance to avoid numerical errors at mirror point
lam2 = lambdam*(1-tiny); % northern mirror latitude
if hemi_symmetric, % only do northern hemisphere
    lam1 = 0;
    Tfactor = 1;
else
    lam1 = -lam2;
    Tfactor = 2;
end

T0 = 1+1/(2*sqrt(3))*log(2+sqrt(3)); % S&L 1.28a
T1 = pi/6*sqrt(2); % S&L 1.28b

T = T0-0.5*(T0-T1)*(sin(a0)+sqrt(sin(a0))); % 1/4 bounce integral of 1
% <f>_b = 1/T(a0)*int_0^lambdam [f(lam)*cos(lam)*sqrt(1+3*sin^2(lam))/cos(a)]dlam

% sin(a) = sin(a0)*sqrt[sqrt(1+3*sin^2(lam))/cos^6(lam)]
% cos(a) = sqrt(1-sin^2(a))
sin_a_lam = @(lam)sin(a0).*sqrt(sqrt(1+3.*sin(lam).^2)./cos(lam).^6);
cos_a_lam = @(lam)sqrt(1-sin_a_lam(lam).^2);


ba = 0; % numerator for now, normalized after loop
g = 0;
for sign_cospa = SIGN_COSPA, % which half-cycle is this?
    ba = ba+quadv(@(lam)wrapper(local,lam,L,MLT,Beq,Bm,sign_cospa,vectorized).*cos(lam).*sqrt(1+3*sin(lam).^2)./cos_a_lam(lam),lam1,lam2);
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

if (N==1) || vectorized,
    out = local(XYZ,Blocal,Bm,maglat,sign_cospa);
else
    out = local(XYZ,Blocal,Bm,maglat(1),sign_cospa);
    out = repmat(out,length(maglat));
    for i = 2:N,
        out(i,:) = local(XYZ(i,:),Blocal(i),Bm,maglat(i),sign_cospa);
    end
end
