% validate bounce_average_trace

alpha0_deg = 89.9*rand(1); % equatorial pitch angle of particle
a0 = alpha0_deg*pi/180;

%f = @(lambda)ones(size(lambda)); % trivial function to integrate
%f = @(lambda)cos_a_lam(lambda); % cosine of local pitch angle
P = rand(2,1)*2-1;
f = @(lambda)polyval(P,lambda); % random polynomial


% find mirror lat:
sina0 = sin(a0);
lambdam = dipole_mirror_latitude(a0,'rad');



% compute field line
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

local = @(XYZ,Blocal,Bm,maglat,sign_cospa)f(maglat);

f_bavd = bounce_average_dipole(L,phi*12/pi,alpha0_deg,local);
f_bavd1 = bounce_average_dipole(L,phi*12/pi,alpha0_deg,local,'symmetric');
f_bavd2 = bounce_average_dipole(L,phi*12/pi,alpha0_deg,local,'hemi_symmetric');
f_bavd3 = bounce_average_dipole(L,phi*12/pi,alpha0_deg,local,'symmetric','hemi_symmetric');

f_bavt = bounce_average_trace([X,Y,Z],B,local,sign(lambda));
f_bavt1 = bounce_average_trace([X,Y,Z],B,local,sign(lambda),'symmetric');
[f_bavd f_bavt ; f_bavd1 f_bavt1 ; f_bavd2 f_bavd3]
