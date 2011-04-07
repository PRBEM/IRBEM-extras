function [ba,denom] = bounce_average_trace(XYZ,Blocal,local,varargin)
% ba = bounce_average_trace(XYZ,Blocal,local,...)
% bounce average "local" over traced field line
% XYZ - Nx3 points along field line, any cartesian coordinates, RE
% Blocal - local magnetic field, nT, two options;
%  Nx1 field strength at each point on field line
%  Nx3 field vector at each point on field line, any cartesian coordinates
% local - what to average, two possibilities
%  NxM doubles - provides local @ XYZ, Blocal
%  @function(XYZ,Blocal,Bm,sign_pa) - evaluated at each point on field line
%    Bm is mirror field strength in nT
%    can return scalar or 1 x M vector
%    sign_pa is sign of cos(local pitch angle) (will do for -1 and 1 to get
%    full bounce orbit)
% options:
% [bint,denom] = bounce_average_trace(...,'no_avg') - returns bounce integral
%  and denominator, does not naverage (bav = bint/denom)
%  NOTE: the true bounce integral is bint*v, where v is the velocity
%  and denom*v is the true full bounce period
%
% [...] = bounce_average_trace(...,'symmetric') - local does not depend
%  on sign of the pitch angle (i.e., northward and southward legs are
%  identical)

% ba = [int_s1^s2 local(XYZ,B) *ds / sqrt(Bmirror - B)] /
%      [int_s1^s2 ds / sqrt(Bmirror - B)]
% bounce averaging involves a singularity
% we handle this singularity by representing B and local
% linearly between the fiducial points

force_symmetric = false; % user tells us it's symmetric?
do_avg = true; % do bounce average? vs integral?
i = 1;
while i <= length(varargin),
    switch(lower(varargin{i})),
        case 'no_avg',
            do_avg = false;
        case 'symmetric',
            force_symmetric = true;
        otherwise
            error('Unknown argument "%s"',varargin{i});
    end
end

if size(Blocal,2)==3,
    Bmag = sqrt(sum(Blocal.^2,2));
else
    Bmag = Blocal;
end

N = length(Bmag);
Bm = max(Bmag); % mirror field strength

if isnumeric(local), % local is a table of numbers for one half-bounce
    loc_func = @(XYZ,Blocal,Bm,sign_pa,i)local(i,:);
    symmetric = true;
else % local a function pointer, need to do both half-bounces
    loc_func = @(XYZ,Blocal,Bm,sign_pa,i)local(XYZ,Blocal,Bm,sign_pa);
    symmetric = false; % assume asymmetric, but may be overriden by force_symmetric
end

if force_symmetric,
    symmetric = force_symmetric;
end

if symmetric,
    SIGN_PA = 1; % double a single half-bounce
    double_result = true;
else
    SIGN_PA = [-1,1]; % do both half-bounces
    double_result = false;
end

% next compute ba' = pathint local/sqrt(Bm-B), 
% and g' = pathint 1/sqrt(Bm-B)
% (we want pathint local/sqrt(1-B/Bm) but we do it this way for simpler
% numerics)
% bint = ba'/sqrt(Bm)
% denom = g'/sqrt(Bm)


ba = 0; % numerator for now, normalized after loop
g = 0; % denominator
for sign_pa = SIGN_PA, % which half-cycle is this?
    % cache i=1 point
    B2  = Bmag(1);
    loc2 = loc_func(XYZ(1,:),Blocal(1,:),Bm,1);
    for i = 2:N,
        B1 = B2;
        loc1 = loc2;
        B2  = Bmag(i);
        loc2 = loc_func(XYZ(i,:),Blocal(i,:),Bm,sign_pa,i);
        ds = sqrt(sum((XYZ(i,:)-XYZ(i-1,:)).^2));
        
        % want int(local/sqrt(Bm-B)*ds,0,s)
        % represent as int((c+d*s)/sqrt(a-b*s),s,0,ds)
        % B = B1+(B2-B1)/ds*s
        % Bm-B = a-b*s = Bm - B1 - (B2-B1)/ds*s
        % local = loc1 + (loc2-loc1)/ds*s = c +d*s
        a = Bm - B1;
        b = (B2-B1)/ds;
        c = loc1;
        d = (loc2-loc1)/ds;
        
        % int((c+d*s)/sqrt(a-b*s),s,0,ds)
        % piecewise([a = 1 and b = 1 and ds = 1, 2*c + (4*d)/3],
        %  [a <> 0 and a = b*ds or a = 0 and b = 0, (2*a^(1/2)*(2*a*d + 3*b*c))/(3*b^2)],
        %  [a <> 0 and a <> b*ds or a = 0 and b <> 0, - (d*((4*a*(a - b*ds)^(1/2))/3 - (4*a^(3/2))/3 + (2*b*ds*(a - b*ds)^(1/2))/3))/b^2 - (c*(2*(a - b*ds)^(1/2) - 2*a^(1/2)))/b],
        %  [(a <> 0 and a = b*ds or a = 0 and b = 0) and (a <> 0 and a <> b*ds or a = 0 and b <> 0), (2*a^(1/2)*c)/b - (d*((4*a*(a - b*ds)^(1/2))/3 - (4*a^(3/2))/3 + (2*b*ds*(a - b*ds)^(1/2))/3))/b^2])

        if b==0, % two points with same B, e.g., symmetric about equator
            % int((c+d*s)/sqrt(a),s,0,ds) = (c*ds+d*ds^2/2)/sqrt(a)
            dba = (c*ds+d*ds^2/2)/sqrt(a);
            dg = c*ds/sqrt(a);            
        elseif (a==1) && (b==1) && (ds==1),
            dba = 2*c+(4*d)/3;
            dg = 2;
        elseif (a ~= 0) && (a == b*ds) || (a == 0) && (b == 0),
            dba = (2*a^(1/2)*(2*a*d + 3*b*c))/(3*b^2);
            dg = (2*a^(1/2))/b; % c=1,d=0
        elseif (a ~= 0) && (a ~= b*ds) || (a == 0) && (b ~= 0),
            dba = - (d*((4*a*(a - b*ds)^(1/2))/3 - (4*a^(3/2))/3 + (2*b*ds*(a - b*ds)^(1/2))/3))/b^2 - (c*(2*(a - b*ds)^(1/2) - 2*a^(1/2)))/b;
            dg = - ((2*(a - b*ds)^(1/2) - 2*a^(1/2)))/b; % c=1,d=0
        else % (a <> 0 and a = b*ds or a = 0 and b = 0) and (a <> 0 and a <> b*ds or a = 0 and b <> 0)
            dba = (2*a^(1/2)*c)/b - (d*((4*a*(a - b*ds)^(1/2))/3 - (4*a^(3/2))/3 + (2*b*ds*(a - b*ds)^(1/2))/3))/b^2;
            dg = (2*a^(1/2))/b; % c=1,d=0
        end
        ba = ba+dba;
        g = g+dg;
    end
end

if double_result,
    ba = ba*2;
    g = g*2;
end

% now normalize: to bounce average
% or to bounce int (without 1/v term)

if do_avg,
    ba = ba/g; % normalize: ba' to ba
else
    ba = ba/sqrt(Bm); % ba' to bint
end

denom = g/sqrt(Bm);
