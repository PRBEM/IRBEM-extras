function [Daa,Dap,Dpp] = Daa_FA_local(species,E,alpha,L,MLT,B,Beq,hemi,wave_model,varargin)
% [Daa,Dap,Dpp] = Daa_FA_local(species,E,alpha,L,MLT,B,Beq,hemi,wave_model)
% compute local diffusion coefficients
% for field-aligned waves. Coefficients are momentum-normalized
% Daa is rad^2/s
% Dap is rad/s, actually Dap/p
% Dpp is 1/s, actually Dpp/p^2
%
% After Summers, 2005, JGR, Quasi-linear diffusion coefficients...
% rewritten with a lot of comparisons to daa-fa-ucla provided
% by Y. Shprits.
% See also Shprits et al, 2006, GRL, Bounce-averaged diffusion
%    coefficients...
% species - particle species, string or struct
%   as string: 'e' (electron), 'p' (proton)
%   as struct:
%      .m0c2 - rest mass, MeV
%      .q - charge in units of proton charge
%      .lambda - species-specific lambda (-1 for e-, me/mp for protons)
% Energy - particle energy in MeV
% alpha - local pitch angle, degrees
% L - L shell (McIlwain) in RE
% MLT - Magnetic Local Time in Hours
% B - local magnetic field strength, nT
% Beq - equatorial magnetic field strength, nT
% hemi - hemisphere +1 for northern, -1 for southern
% wave_model - structure
%  .mode - 'R' or 'L' - wave mode (polarization)
%  .composition' - [H+ fraction,He+ fraction,O+ fraction]
%   (provides ion composition, must sum to 1. Default = [1 0 0])
%  .normalization - string
%    changes meaning of omega_m, domega, omega1, omega2
%    to be in terms of
%    'Omega_e' local cold electron gyro
%    'Omega_e_eq' equatorial cold electron gyro
%    'Omega_p' local cold proton gyro
%    'Omega_p_eq' equatorial cold proton gyro
%    'Omega_He+', 'Omega_He+_eq', 'Omega_O+', 'Omega_O+_eq' are also supported
%    'Omega_s' local cold species gyro
%    'Omega_s_eq' equatorial cold species gyro
%  .omega_pe_normalization - string
%     changes meaning of omega_pe. Uses same values as .normalization
%  [each of the following can be a scalar or
%   a function handle that takes (L,MLT,maglat), maglat in degrees]
%  .directions -
%    'F' forward (k || to B)
%    'B' backward (k || to -B)
%    'FB' both
%    'P' poleward (k || to hemi*B)
%    'E' equatorward (k || to -hemi*B)
%  .dB - function provides peak wave amplitude, nT
%  .omega_m - angular frequency of peak wave power, rad/sec
%  .domega - angular frequency width of wave power, rad/sec
%  .omega1, .omega2 - lower and upper bounds of wave power Gaussian, rad/sec
%   ** can supply .sigma instead, where sigma = (omega2-omega1) / (2 domega)
%  .N0 = number density, #/cm^3
%   ** can supply .alpha_star = (Omega_e/omega_pe)^2 (local)
%   ** can supply .omega_pe (rad/s, unless omega_pe_normalization provided)
%
% options:
% [...] = Daa_FA_local(...,'method','Summers2007'); use Summers 2007 method
%   (multispecies plasma)
% [...] = Daa_FA_local(...,'method','Summers2005'); use Summers 2005 method
%   (Hydrogen plasma only, maybe some errors)
% [...] = Daa_FA_local(...,'method','auto'); select Summers 2005 for
% Hydrogen plasma, Summers 2007 for multispecies (default)
% [...] = Daa_FA_local(...,'usenu'); use nu versus rho
%   (nu is 2005 normalization and appears to be incorrect, rho is 2007
%   normalization)
% D = Daa_FA_local(...,'join_outputs'): D = [Daa,Dap,Dpp];

method = 'auto';
join_outputs = false;
last_root_only = false; % this is for debugging only
usenu = false;
i = 1;
while i <= length(varargin),
    switch(lower(varargin{i})),
        case 'method',
            i = i+1;
            method = lower(varargin{i});
        case 'join_outputs',
            join_outputs = true;
        case 'last_root_only', % for debugging only
            last_root_only = true;
        case 'usenu', % for debugging only
            usenu = true;
        otherwise
            error('Unknown argument "%s"',varargin{i});
    end
    i = i+1;
end

if isfield(wave_model,'composition'),
    composition = wave_model.composition;
else    
    composition = [1 0 0]; % H+, He+, O+, default is all Hydrogen
end


if strcmpi(method,'auto'),
    if composition(1) == 1, % Hydrogen plasma
        method = 'Summers2005';
    else
        method = 'Summers2007'; % multispecies
    end
end

maglat = BBeq_to_maglat(B/Beq)*hemi;

directions = upper(evaluate_scalar_or_handle(wave_model.directions,L,MLT,maglat));
switch(directions),
    case 'F', % forward (k || to B)
        yrange = [0 inf];
    case 'B', % backward (k || to -B)
        yrange = [-inf 0];
    case 'P', % poleward (k || to hemi*B)
        yrange = sort([0 hemi*inf]);
    case 'E', % equatorward (k || to -hemi*B)
        yrange = sort([0 -hemi*inf]);
    case {'FB','PE','BF','EP'}, % both
        yrange = [-inf inf];
end

Daa = 0;
Dap = 0;
Dpp = 0;

% test for nonzero, finite wave field
dB = evaluate_scalar_or_handle(wave_model.dB,L,MLT,maglat);

if dB==0,
    return; % zeros
elseif ~isfinite(dB),
    if join_outputs,
        Daa = nan(1,3);
    else
        Daa = nan;
        Dap = nan;
        Dpp = nan;
    end
    return
end

species = species2struct(species);
electron = species2struct('e');
proton = species2struct('p');

% standard frequencies
Omega_e = abs(electron.qC) * B*1e-9 / electron.m0kg; % local cold electron gyro
Omega_e_eq = abs(electron.qC) * Beq * 1e-9 / electron.m0kg; % equatorial cold electron gyro
Omega_p = proton.qC * B*1e-9 / proton.m0kg; % local cold proton gyro
Omega_p_eq = proton.qC * Beq * 1e-9 / proton.m0kg; % equatorial cold proton gyro
Omega_s = abs(species.qC) * B*1e-9 / species.m0kg; % local cold species gyro
Omega_s_eq = abs(species.qC) * Beq * 1e-9 / species.m0kg; % equatorial cold species gyro


if isfield(wave_model,'normalization'),
    freq_norm = interpret_normalization(wave_model.normalization,Omega_e,Omega_e_eq,Omega_p,Omega_p_eq,Omega_s,Omega_s_eq);
else
    freq_norm = 1;
end

if isfield(wave_model,'omega_pe_normalization'),
    omega_pe_freq_norm = interpret_normalization(wave_model.omega_pe_normalization,Omega_e,Omega_e_eq,Omega_p,Omega_p_eq,Omega_s,Omega_s_eq);
else
    omega_pe_freq_norm = 1;
end


omega_m = evaluate_scalar_or_handle(wave_model.omega_m,L,MLT,maglat)*freq_norm;
domega = evaluate_scalar_or_handle(wave_model.domega,L,MLT,maglat)*freq_norm;

% use sigma or omega1,omega2 based on Summers before (30)
if isfield(wave_model,'sigma'),
    sigma = evaluate_scalar_or_handle(wave_model.sigma,L,MLT,maglat);
    omega1 = omega_m - domega*sigma;
    omega2 = omega_m + domega*sigma;
else
    omega1 = evaluate_scalar_or_handle(wave_model.omega1,L,MLT,maglat)*freq_norm;
    omega2 = evaluate_scalar_or_handle(wave_model.omega2,L,MLT,maglat)*freq_norm;
    sigma = (omega2-omega1)/ (2*domega);
end
R = (dB/B)^2; % after Summers (35)
Erel = E/species.m0c2; % E relative to rest mass

epsilon = electron.m0c2 / proton.m0c2; % ratio of electron to proton mass, after Summers (22)

% set s, Summers after (19)
switch(lower(wave_model.mode)),
    case 'r',
        s = 1;
        xrange = [0 1]; % paragraph 42
    case 'l',
        s = -1;
        xrange = [0 epsilon];  % paragraph 42
    otherwise
        error('Unknown wave mode "%s", expect R or L',wave_model.mode);
end


% compute local electron plasma frequency omega_pe, rad/s
if isfield(wave_model,'alpha_star'), % use alpha_star
    alpha_star = evaluate_scalar_or_handle(wave_model.alpha_star,L,MLT,maglat);
    omega_pe = sqrt(Omega_e^2 / alpha_star); % Summers (22)
elseif isfield(wave_model,'N0'),
    epsilon0 = 8.854187817e-12; % F/m = C^2 s^2  / kg / m^3 - permitivity of free space
    N0 = evaluate_scalar_or_handle(wave_model.N0,L,MLT,maglat);
    omega_pe = sqrt((N0*1e2^3)*electron.qC^2 / epsilon0 / electron.m0kg);
elseif isfield(wave_model','omega_pe'),
    omega_pe = evaluate_scalar_or_handle(wave_model.omega_pe,L,MLT,maglat)*omega_pe_freq_norm;
else
    error('No specification of omega_pe, alpha_star, or N0 to provide omega_pe');
end

beta = sqrt(Erel*(Erel+2))/(Erel+1); % v/c, Summers after A1
gamma = Erel+1; % E/m0c2 + 1, Summers after A1
mu = cosd(alpha); % alpha is in degrees, Summers after (6)

% for a Hydrogen plasma
alpha_star = (Omega_e / omega_pe)^2;  % Summers (22)
a = s*species.lambda/gamma; % Summers after (24)
x_m = omega_m / Omega_e; % Summers after (35)
dx = domega / Omega_e; % Summers after (35)
x1 = omega1 / Omega_e;
x2 = omega2 / Omega_e;
if usenu,
    nu = sqrt(pi)*erf(sigma); % after Summers (32)
    rho = nu;
else
    rho = sqrt(pi)/2*(erf((omega_m-omega1)/domega)+erf((omega2-omega_m)/domega)); % Summers 2007, (3)
    % rho is the more general normalization for 2005 or 2007 method. nu appears to be
    % right only when omega1 and omega2 are symmetric around omega_m    
end


switch(lower(method)),
    case {'summers2005','2005',2005},        
        if composition(1) ~= 1,
            warning('Summers2005 method ignores composition');
        end
        b = (1+epsilon) / alpha_star; % Summers after (25)
        F_func = @(x,y)F2005(x,y,s,epsilon,b);
        roots_func = @(x1,x2)roots2005(x1,x2,s,beta,mu,epsilon,a,b);
        y_func = @(x)y2005(x,b,s,epsilon);
    case {'summers2007','2007',2007},
        masses = [1;4;16]; % H, He, O
        eta = composition(:)/sum(composition); % ensure it sums to 1
        F_func = @(x,y)F2007(x,y,s,epsilon,alpha_star,eta,masses);
        roots_func = @(x1,x2)roots2007(x1,x2,s,epsilon,alpha_star,mu,beta,a,eta,masses);
        y_func = @(x)y2007(x,alpha_star,s,epsilon,eta,masses);
    otherwise
        error('Unknown method "%s"',method);
end        

if mu ~= 0,
    % Evaluate Summers 2005 (33)-(35), or equivalent Summers 2007 (5)-(7)
    rr = roots_func(x1,x2);
    rr = rr((rr>=xrange(1)) & (rr < xrange(end)));
    if last_root_only && ~isempty(rr),
        rr = rr(end);
    end
    for iroot = 1:length(rr),
        x = rr(iroot); % omega(i) / Omega_e, Summers after (35)
        y = (x+a)/(beta*mu); % Summers (24)
        if (y>=yrange(1)) && (y <= yrange(end)),
            F_tmp = F_func(x,y);
            common_factor = R * abs(F_tmp) / ...
                (dx*abs(beta*mu-F_tmp)) * exp(-((x-x_m)/dx)^2);
            a_factor = (1-x*mu/(y*beta));
            p_factor = (x/y);
            Daa =  Daa + common_factor*a_factor^2;
            Dap =  Dap + common_factor*a_factor*p_factor;
            Dpp =  Dpp + common_factor*p_factor^2;
        end
    end
    common_prefactor = (pi/2)*(1/rho)*Omega_s^2/Omega_e/(Erel+1)^2;
    p_prefactor = sind(alpha)/beta; % alpha in deg
    Daa = common_prefactor*Daa; % Summers (33)
    Dap = -common_prefactor*Dap*p_prefactor; % Summers (34)
    Dpp = common_prefactor*Dpp*p_prefactor^2; % Summers (35)
else
    % special case for alpha=90, given in Summers paragraph [45]
    % note, resonates with forward and backward going waves, so
    % sum over +y0 and -y0
    if ((s==1) && (species.lambda==-1)) ... % Summers case (i), A2 (electrons with R-mode waves)
            || ((s==-1) && (species.lambda==epsilon)), % Summers case (ii), A3 (protons with L-mode waves)
        x0 = -a; % only finite y in (24) for mu=0 is x = -a
        y0 = y_func(x0);
    else % other cases require particle to overtake wave, which can't happen for mu=0, alpha=90, v|| = 0
        x0 = nan;
    end
    if (x0 > x1) && (x0 < x2),
        Daa_tmp = pi/2/rho*Omega_s^2/Omega_e/(Erel+1)^2*R/dx*exp(-((x0-x_m)/dx)^2); % Summers (36)
        for signy = [-1,1],
            y = signy*y0;
            if (y>=yrange(1)) && (y <= yrange(end)),
                Daa = Daa+Daa_tmp;
                Dap = Dap+Daa_tmp/beta*(x0/y);  % Summers (37)
                Dpp = Dpp+Daa_tmp/beta^2*(x0/y)^2; % Summers (38)
            end
        end
    end
end
        
if join_outputs,
    Daa = [Daa,Dap,Dpp];
end

function y = y2005(x,b,s,epsilon)
y = abs(x)*sqrt(1-b/((x-s)*(x+s*epsilon))); % (25), 2005 (also A2, A3)

function y = y2007(x,alpha_star,s,epsilon,eta,masses)
y = sqrt(x^2+x/alpha_star*(1/(s-x)-epsilon*sum(eta./(masses*x+s*epsilon)))); % (11), 2007

function y = evaluate_scalar_or_handle(x,varargin)
% y = evaluate_scalar_or_hande(x,...)
% if x is a scalar, returns it
% otherwise, evaluates it as a function y = x(...)

if isa(x,'function_handle'),
    y = x(varargin{:});
else
    y = x;
end

function maglat = BBeq_to_maglat(BB0)
% convert B/Bequator to magnetic latitude (degrees)
persistent table
if isempty(table),
    table.maglat_deg = (0:0.1:89.9)';
    table.maglat_bb0 = (1+3*sind(table.maglat_deg).^2).^(1/2)./cosd(table.maglat_deg).^6;
end
maglat = interp1(table.maglat_bb0,table.maglat_deg,BB0,'linear');

function species = species2struct(species)
% populate species structure
% m0c2 - rest mass, MeV
% m0kg - rest mass, kg
% q - charge, in units of proton charge
% qC - charge in C
epsilon = 0.510998910 / 938.272013 ; % ratio of electron to proton mass, after Summers (22)
% lam is defined after Summers (24)
if ischar(species),
    switch(lower(species)),
        case {'e','e-','beta','beta-','electron','ele'},
            species = struct('m0c2',0.510998910,'q',-1,'lambda',-1);
        case {'e+','beta+','positron'},
            species = struct('m0c2',0.510998910,'q',+1);
        case {'p','proton','p+','h+'},
            species = struct('m0c2',938.272013,'q',+1,'lambda',epsilon);
        case {'p-','antiproton'},
            species = struct('m0c2',938.272013,'q',-1);
        otherwise
            error('Unknown species "%s"',species);
    end
end
species.qC = species.q*1.602176487e-19;
species.m0kg = species.m0c2*1.782661758877380e-030;


function rr = roots2005(x1,x2,s,beta,mu,epsilon,a,b)
% real roots of 4-th order polynomial- Summers, (A1)
% x1, x2 - lower and upper bounds of frequency range
%   in units of the *local* electron gyro frequency
% Erel = energy in units of rest mass
% s = sign of particle charge
% lam_sp = -1 for electrons, epsilon for protons
% beta - v/c
% mu - cos local pitch angle
% epsilon - ratio of proton to electron mass
% a - s * lambda / gamma
% b - (1 + epsilon) / alpha_star

% turn common divisor a1-a4 into a0 =  (1-(beta*mu)^2)
a0 = 1-(beta*mu)^2;
a1 = 2*a+ s*(epsilon-1)*(1 - (beta*mu)^2);
a2 = a^2 + 2*a*s*(epsilon-1) - epsilon + (beta*mu)^2*(b+epsilon);
a3 = a^2*s*(epsilon-1)-2*a*epsilon;
a4 = -a^2*epsilon;

% fprintf('x1=%g, x2=%g\n',x1,x2);
% fprintf('s=%g,beta=%g,mu=%g,epsilon=%g,a=%g,b=%g\n',s,beta,mu,epsilon,a,b);
% poly = [a0 a1 a2 a3 a4]  / a0 ; % Summers's format
% a3 / a0
% a4 / a0

Croots = roots([a0 a1 a2 a3 a4]);
rr = Croots((Croots > x1) & (Croots < x2) & ~imag(Croots));
% validated up to this point against Shprits

function rr = roots2007(x1,x2,s,epsilon,alpha_star,mu,beta,a,eta,masses)
% implements Summers 2007 (18)
% returns unique real roots in x1,x2
xi = beta*mu;
A = 64+epsilon*(64*eta(1)+16*eta(2)+4*eta(3));
B = s*epsilon*(84-((64-20*epsilon)*eta(1)+(16-17*epsilon)*eta(2)+(4-5*epsilon)*eta(3)));
C = epsilon^2*(21-((20-epsilon)*eta(1)+(17-epsilon)*eta(2)+(5-epsilon)*eta(3)));
A0 = 64*(1-xi^2);
A1 = s*(84*epsilon-64)+128*a+xi^2*s*(64-84*epsilon);
A2 = epsilon*(21*epsilon-84)+a*s*(168*epsilon-128)+64*a^2 ...
    + xi^2*epsilon*(84-21*epsilon)+A*xi^2/alpha_star;
A3 = s*epsilon^2*(epsilon-21)+a*epsilon*(42*epsilon-168)+a^2*s*(84*epsilon-64) ...
    +xi^2*epsilon^2*s*(21-epsilon)+B*xi^2/alpha_star;
A4 = s*a*epsilon^2*(2*epsilon-42)+a^2*epsilon*(21*epsilon-84)-epsilon^3*(1-xi^2) ...
    +C*xi^2/alpha_star;
A5 = a*epsilon^2*(a*s*(epsilon-21)-2*epsilon);
A6 = -a^2*epsilon^3;

P = [A0 A1 A2 A3 A4 A5 A6];

Croots = roots(P);
rr = Croots((Croots > x1) & (Croots < x2) & ~imag(Croots));

% remove noise roots that don't solve y9^2-y11^2
for i = 1:length(rr),
    x = rr(i);
    y9 = (x+a)/(beta*mu); % (9), 2007
    y11 = y2007(x,alpha_star,s,epsilon,eta,masses); % (11), 2007
    if abs(y9^2-y11^2)>1e-6,
        rr(i) = nan;
    end
end
rr = rr(isfinite(rr));


function out = F2005(x,y,s,epsilon,b)
% Summers (C1)
% dispersion relation

c1 = 2*s*(epsilon-1);
c2 = 1-4*epsilon+epsilon^2;
c3 = s*(1-epsilon)*(b+4*epsilon)/2;
c4 = epsilon*(b+epsilon);

g = polyval([1 c1 c2 c3 c4],x);

out = y*(x-s)^2*(x+s*epsilon)^2 / (x * g);
% validated against Shprits

function out = F2007(x,y,s,epsilon,alpha_star,eta,masses)
% dispersion relation for Summers 2007, (11)
% dxdy

dydx = 2*x*alpha_star + 1/(s-x) + x/(s-x)^2; % constant and electron terms
dydx = dydx - epsilon*sum(eta./(masses*x+s*epsilon)) + x*epsilon*sum(eta.*masses./(masses*x+s*epsilon).^2); % ion terms
dydx = dydx/(2*y*alpha_star);

out = 1/dydx; % dx/dy


function freq_norm = interpret_normalization(normalization,Omega_e,Omega_e_eq,Omega_p,Omega_p_eq,Omega_s,Omega_s_eq)
switch(normalization),
    case 'Omega_e',
        freq_norm = Omega_e;
    case 'Omega_e_eq',
        freq_norm = Omega_e_eq;
    case 'Omega_p',
        freq_norm = Omega_p;
    case 'Omega_p_eq',
        freq_norm = Omega_p_eq;
    case 'Omega_He+',
        freq_norm = Omega_p/4;
    case 'Omega_He+_eq',
        freq_norm = Omega_p_eq/4;
    case 'Omega_O+',
        freq_norm = Omega_p/16;
    case 'Omega_O+_eq',
        freq_norm = Omega_p_eq/16;
    case 'Omega_s',
        freq_norm = Omega_s;
    case 'Omega_s_eq',
        freq_norm = Omega_s_eq;
    otherwise
        error('Unknown normalization "%s"',wave_model.normalization);
end
