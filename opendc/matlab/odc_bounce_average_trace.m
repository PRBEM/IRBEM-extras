function [ba,denom] = odc_bounce_average_trace(XYZ,Blocal,local,hemi,varargin)
% ba = odc_bounce_average_trace(XYZ,Blocal,local,hemi,...)
% bounce average "local" over traced field line
% XYZ - Nx3 points along field line, any cartesian coordinates, RE
% Blocal - local magnetic field, nT, two options;
%  Nx1 field strength at each point on field line
%  Nx3 field vector at each point on field line, any cartesian coordinates
% local - what to average, two possibilities
%  NxM doubles - provides local @ XYZ, Blocal
%  @function(XYZ,Blocal,Bm,maglat,sign_cospa) - evaluated at each point on field line
%    Bm is mirror field strength in nT
%    maglat is signed magnetic latitude (from B/Bmin), degrees
%    can return scalar or 1 x M vector
%    sign_cospa is sign of cos(local pitch angle) (will do for -1 and 1 to get
%    full bounce orbit)
% hemi Nx1 - sign of hemisphere +1 for northern, -1 for southern
% options:
% [bint,denom] = bounce_average_trace(...,'no_avg') - returns bounce integral
%  and denominator (in RE), does not average (bav = bint/denom)
%  NOTE: the true bounce integral is bint*v, where v is the velocity
%  and denom*v is the true full bounce period
%
% [...] = bounce_average_trace(...,'symmetric') - local does not depend
%  on sign of the pitch angle (i.e., northward and southward legs are
%  identical)
% [...] = bounce_average_trace(...,'Bm',Bm) - provide Bmirror
%   (otherwise it's max(|B|)

% ba = [int_0^(2pi) local(...) ds/dB cos(psi) dpsi] /
%      [int_0^(2pi) ds/dB cos(psi) dpsi]
% bounce averaging involves a singularity
% we handle this singularity by representing B and local
% linearly between the fiducial points
% Implements psi transform from Orlova and Shprits, 2011, Phys. Plasmas

% O&S denotes equation or info from Orlova and Shprits, 2011

force_symmetric = false; % user tells us it's symmetric?
do_avg = true; % do bounce average? vs integral?
Bm = nan;
i = 1;
while i <= length(varargin),
    switch(lower(varargin{i})),
        case 'no_avg',
            do_avg = false;
        case 'symmetric',
            force_symmetric = true;
        case 'bm',
            i = i+1;
            Bm = varargin{i};
        otherwise
            error('Unknown argument "%s"',varargin{i});
    end
    i = i+1;
end

if size(Blocal,2)==3,
    Bmag = sqrt(sum(Blocal.^2,2));
else
    Bmag = Blocal;
end

N = length(Bmag);
if ~isfinite(Bm),
    Bm = max(Bmag); % mirror field strength
end
Beq = min(Bmag);
a0_deg = asind(sqrt(Beq./Bm));
util = odc_util; % load utility functions and constants
maglat = util.BB0toMagLat(Bmag./Beq).*hemi;

if isnumeric(local), % local is a table of numbers for one half-bounce
    loc_func = @(XYZ,Blocal,Bm,maglat,sign_cospa,i)local(i,:);
    symmetric = true;
else % local a function pointer, need to do both half-bounces
    loc_func = @(XYZ,Blocal,Bm,maglat,sign_cospa,i)local(XYZ,Blocal,Bm,maglat,sign_cospa);
    symmetric = false; % assume asymmetric, but may be overriden by force_symmetric
end

if force_symmetric,
    symmetric = force_symmetric;
end

if symmetric,
    SIGN_COSPA = 1; % double a single half-bounce
    double_result = true;
else
    SIGN_COSPA = [-1,1]; % do both half-bounces
    double_result = false;
end

% next compute ba' = pathint local/sqrt(Bm-B),
% and g' = pathint 1/sqrt(Bm-B)
% (we want pathint local/sqrt(1-B/Bm) but we do it this way for simpler
% numerics)
% bint = ba'/sqrt(Bm)
% denom = g'/sqrt(Bm)



loc1 = loc_func(XYZ(1,:),Blocal(1,:),Bm,maglat(1),-1,1);
loc = nan(N,size(loc1,2));
ba = zeros(1,size(loc1,2)); % numerator for now, normalized after loop
psi = asin(min(1,sqrt(1-Bmag/Bm)/cosd(a0_deg))); % min(1,...) handles round-off error
% psi is presently only on [0,pi/2], doubles back in North hemi
iN = hemi>=0; % northern hemisphere points must be put on interval [pi/2,pi]
psi(iN) = pi-psi(iN); % handle psi quadrant
% compute dsdB*cos(psi)
dsdBcospsi = nan(N,1);
if N==1,
    dsdBcospsi = 1;
else
    dsdBcospsi(1) = cos(psi(1))*sqrt(sum((XYZ(2,:)-XYZ(1,:)).^2))./(Bmag(2)-Bmag(1));
    for i = 2:(N-1),
        if Bmag(i+1)==Bmag(i-1),
            dsdBcospsi(i) = -tand(a0_deg)^2/2/Beq*sqrt(sum((XYZ(i+1,:)-XYZ(i-1,:)).^2))./(psi(i+1)-psi(i-1)); % O&S 7
        else
            dsdBcospsi(i) = cos(psi(i))*sqrt(sum((XYZ(i+1,:)-XYZ(i-1,:)).^2))./(Bmag(i+1)-Bmag(i-1));
        end
    end
    dsdBcospsi(N) = cos(psi(N))*sqrt(sum((XYZ(N,:)-XYZ(N-1,:)).^2))./(Bmag(N)-Bmag(N-1));
end

dsdBcospsi = abs(dsdBcospsi); % positive definite

for sign_cospa = SIGN_COSPA, % which half-cycle is this?
    for i = 1:N,
        loc(i,:) = loc_func(XYZ(i,:),Blocal(i,:),Bm,maglat(i),sign_cospa,i);
    end
    for i = 1:size(loc,2),
        ba(i) = ba(i) + trapz(psi,dsdBcospsi.*loc(:,i));
    end
end
g = trapz(psi,dsdBcospsi)*length(SIGN_COSPA); % compute denominator

if double_result,
    ba = ba*2;
    g = g*2;
end

% now normalize: to bounce average
% or to bounce int (without 1/v term)
prefactor = 2*cosd(a0_deg)*Bm;

if do_avg,
    ba = ba/g; % normalize: ba' to ba
else
    ba = ba*prefactor; % ba' to bint
end

denom = g*prefactor;

