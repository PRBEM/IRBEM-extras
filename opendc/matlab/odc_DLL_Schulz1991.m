function [DLLM,DLLE] = odc_DLL_Schulz1991(L,alpha0_deg,MeV,tau)
% [DLLM,DLLE] = odc_DLL_Schulz1991(L,alpha0_deg,MeV)
% [DLLM,DLLE] = odc_DLL_Schulz1991(L,alpha0_deg,MeV,tau)
% returns electrostatic (DLLE) and electromagnetic (DLLM)
% diffusion coefficients from Schulze 1991 (Geomagnetism)
% in units of 1/day, for electrons
% L - L shell
% alpha0_deg - equatorial pitch angle, degrees
% MeV - particle energy, MeV
% tau - decay time of electrostatic impulse, seconds (default is 1200)
% (note sizes must be compatible so that L.*alpha0_deg.*MeV is valid)

if nargin < 4,
    tau = 1200; % seconds, decay time of electrostatic impulse
end

util = odc_util;

y = sind(alpha0_deg);
Q = util.Q(y);
D = util.D(y);
T = util.T(y);

% equation 203
DLLM = 7E-9.*(Q./D./180).^2.*L.^10; 

Td = util.DriftPeriod('e-',MeV,alpha0_deg,L);
Omega3tau = 2*pi./Td.*tau;

gamma = util.MeVtogamma(MeV,'e-');

B0nT = util.dipoleB(1,0,0);
B0G = B0nT/1e5;

E0 = util.mks.electron.m0.*util.mks.c.^2/util.mks.MeV;
M = y.^2.*L.^3.*MeV.*(MeV+2*E0)./(2.*E0.*B0G);

% equation 204 with Z = 1, M0 = 1 GeV/G
M0 = 1000; % MeV/G
DLLE = 1E-10*L.^10.*(gamma.*y.^2.*M0./M).^2.*(T./2./D).^2./(1+Omega3tau.^-2);

