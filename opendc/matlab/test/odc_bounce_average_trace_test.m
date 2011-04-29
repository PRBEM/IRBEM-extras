% validate bounce_average_trace

alpha0_deg = 89.9*rand(1); % equatorial pitch angle of particle
a0 = alpha0_deg*pi/180;

%f = @(lambda)ones(size(lambda)); % trivial function to integrate
P = rand(6,1)*2-1;
f = @(lambda)polyval(P,lambda/90); % random polynomial

util = odc_util; % load utility functions and constants

% find mirror lat:
sina0 = sin(a0);
lambdam = util.dipole_mirror_latitude(a0,'rad');

% compute field line
lambda = linspace(-lambdam,lambdam,200)'; % radians
L = 4.5;
B0 = 30e3; % equatorial field strength at L=1
phi = 37*pi/180; % arbitrary azimuth
MLT = phi*12/pi;
r = L*cos(lambda).^2;
B = B0./r.^3.*sqrt(1+3*sin(lambda).^2);
Beq = B0./L^3;
X = r.*cos(lambda).*cos(phi);
Y = r.*cos(lambda).*sin(phi);
Z = r.*sin(lambda);

fprintf('Within groups, these should all agree to 2-3 digits:\n');


local = @(XYZ,Blocal,Bm,maglat,sign_cospa)f(maglat);

f_bavd = odc_bounce_average_dipole(L,MLT,alpha0_deg,local);
f_bavd1 = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,'symmetric');

f_bavt = odc_bounce_average_trace([X,Y,Z],B,local,sign(lambda));
f_bavt1 = odc_bounce_average_trace([X,Y,Z],B,local,sign(lambda),'symmetric');

fprintf('Results for non-hemi-symmetric tests\n');
[f_bavd f_bavt ; f_bavd1 f_bavt1]

% now make function that's north-south symmetric
local = @(XYZ,Blocal,Bm,maglat,sign_cospa)f(abs(maglat));

f_bavd = odc_bounce_average_dipole(L,MLT,alpha0_deg,local);
f_bavd1 = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,'symmetric');
f_bavd2 = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,'hemi_symmetric');
f_bavd3 = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,'symmetric','hemi_symmetric');

f_bavt = odc_bounce_average_trace([X,Y,Z],B,local,sign(lambda));
f_bavt1 = odc_bounce_average_trace([X,Y,Z],B,local,sign(lambda),'symmetric');

fprintf('Results for hemi-symmetric tests\n');
[f_bavd f_bavt ; f_bavd1 f_bavt1 ; f_bavd2 f_bavd3]

