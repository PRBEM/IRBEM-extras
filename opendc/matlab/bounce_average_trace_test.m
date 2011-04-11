% validate bounce_average_trace

% use expressions in Shpris' thesis, appendix F

alpha0_deg = 89.9*rand(1); % equatorial pitch angle of particle

a0 = alpha0_deg*pi/180;

% s = @(a0)1.3-0.56*sin(a0); % F11 (good)
s = @(a0)1.38-0.32*(sin(a0)+sqrt(sin(a0))); % F12 (better)
% <f>_b = 1/s(a0)*int_0^lambdam [f(lam)*cos(lam)*sqrt(1+3*sin^2(lam))/cos(a)]dlam
% x = cos^2(lam)
% x^6+(3*sin^2(a0))*x - 4*sin^2(a0) = 0

% sin(a) = sin(a0)*sqrt[sqrt(1+3*sin^2(lam))/cos^6(lam)]
% cos(a) = sqrt(1-sin^2(a))

sin_a_lam = @(lam)sin(a0).*sqrt(sqrt(1+3.*sin(lam).^2)./cos(lam).^6);
cos_a_lam = @(lam)sqrt(1-sin_a_lam(lam).^2);

% find mirror lat:
sina0 = sin(a0);
lambdam = dipole_mirror_latitude(a0,'rad');

%f = @(lambda)ones(size(lambda)); % trivial function to integrate
%f = @(lambda)cos_a_lam(lambda); % cosine of local pitch angle
P = randn(5,1);
f = @(lambda)polyval(P,lambda); % random polynomial
f_bav = quadl(@(lam)f(lam).*cos(lam).*sqrt(1+3*sin(lam).^2)./cos_a_lam(lam),0,lambdam*(1-1e-6))/s(a0);
if abs(imag(f_bav)) < abs(real(f_bav))/1e4,
    f_bav = real(f_bav);
else
    error('a0 = %.1f, <f>_b = %g + %gi\n',alpha0_deg,real(f_bav),imag(f_bav));
end

% compute B along field line
RE = 6378;
lambda = linspace(-lambdam,lambdam,200)'; % radians
L = 4.5;
B0 = 30e3; % equatorial field strength at L=1
phi = 37*pi/180; % arbitrary azimuth
r = L*cos(lambda).^2;
B = B0./r.^3.*sqrt(1+3*sin(lambda).^2);
Beq = B0./L^3;
X = r.*cos(lambda).*cos(phi);
Y = r.*cos(lambda).*sin(phi);
Z = r.*sin(lambda);

local = @(XYZ,Blocal,Bm,maglat,sign_cospa)f(BB0toMagLat(Blocal./Beq,'rad'));

f_bavt = bounce_average_trace([X,Y,Z],B,local,sign(lambda));
[f_bav f_bavt]


