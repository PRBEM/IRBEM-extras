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
phi_deg = 37; % arbitrary azimuth
MLT = phi_deg/15; % hours
[B,Bvec,XYZ] = util.dipoleB(L,lambda*180/pi,phi_deg);


fprintf('Within groups, these should all agree to 2-3 digits:\n');


local = @(XYZ,Blocal,Bm,maglat,sign_cospa)f(maglat);

f_bavd = odc_bounce_average_dipole(L,MLT,alpha0_deg,local);
f_bavd1 = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,'symmetric');

f_bavt = odc_bounce_average_trace(XYZ,B,local,sign(lambda));
f_bavt1 = odc_bounce_average_trace(XYZ,B,local,sign(lambda),'symmetric');

fprintf('Results for non-hemi-symmetric tests\n');
[f_bavd f_bavt ; f_bavd1 f_bavt1]

methods = {'quad','quadv','quadl','quadgk','trace'};
for i = 1:length(methods),
    f_bavd = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,'method',methods{i});
    f_bavd1 = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,'symmetric','method',methods{i});
    fprintf('Results for hemi-symmetric tests (method=%s)\n',methods{i});
    [f_bavd f_bavd1]
end

fprintf('----------------------------------------------\n');

% now make function that's north-south symmetric
local = @(XYZ,Blocal,Bm,maglat,sign_cospa)f(abs(maglat));

f_bavd = odc_bounce_average_dipole(L,MLT,alpha0_deg,local);
f_bavd1 = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,'symmetric');
f_bavd2 = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,'hemi_symmetric');
f_bavd3 = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,'symmetric','hemi_symmetric');

f_bavt = odc_bounce_average_trace(XYZ,B,local,sign(lambda));
f_bavt1 = odc_bounce_average_trace(XYZ,B,local,sign(lambda),'symmetric');

fprintf('Results for hemi-symmetric tests\n');
[f_bavd f_bavt ; f_bavd1 f_bavt1 ; f_bavd2 f_bavd3]

methods = {'quad','quadv','quadl','quadgk','trace'};
for i = 1:length(methods),
    f_bavd = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,'method',methods{i});
    f_bavd1 = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,'symmetric','method',methods{i});
    f_bavd2 = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,'hemi_symmetric','method',methods{i});
    f_bavd3 = odc_bounce_average_dipole(L,MLT,alpha0_deg,local,'symmetric','hemi_symmetric','method',methods{i});    
    fprintf('Results for hemi-symmetric tests (method=%s)\n',methods{i});
    [f_bavd f_bavd1 f_bavd2 f_bavd3]
end
