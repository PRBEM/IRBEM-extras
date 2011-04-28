function [DLLM,DLLE] = odc_DLL_BA2000(L,Kp,alpha0_deg,MeV)
% DLLM = odc_DLL_BA2000(L,Kp)
% DLLM = odc_DLL_BA2000(L,Kp,alpha0_deg)
% [DLLM,DLLE] = odc_DLL_BA2000(L,Kp,alpha0_deg,MeV)
% [DLLM,DLLE] = odc_DLL_BA2000(L,Kp,alpha0_deg,-Td)
% compute radial diffusion coefficient
% according to Brautigam and Albert, 2000
% with (optional) Schulz & Lanzerotti angular factor for DLLM *** check
% L - dipole L shell
% Kp - Kp magnetic index
% alpha0_deg - equatorial pitch angle, degrees (default is 90)
% MeV - particle energy, MeV
% Td - drift period, seconds
% 
% DLLM - magnetic radial diffusion coefficient (1/day)
% DLLE - electric radial diffusion coefficient (1/day)
%

Kp = min(max(Kp,1),6); % B&A applies only for Kp = 1 to 6

% B&A is for equatorially mirroring particles
% we'll apply the appropriate angular correction per S&L
% DLL_Kp = @(alpha0,E,L,Kp)DLL_BA_90(Kp,L).*QD2(alpha0);
if nargin >= 3,
    util = odc_util;
    y = sind(alpha0_deg);
    %(Q(y)/D(y))^2,  % angular dependence of DLL M, S&L, 3.13 
    % *** check that applies to DLL M only according to S&L
    QD2 = (util.Q(y)./util.D(y)).^2; 
else
    alpha0_deg = 90;
    QD2 = 1; % no angular correction    
end
DLLM = 10.^(0.506*Kp-9.325).*L.^10.*QD2; % D^M_LL Brautigam and Albert, 2000, eq. 6

if nargout >= 2,
    warning('B&A DLLE is not working -- it does not match published case');
    % B&A 2000, eqn 4, 5
    util = odc_util;
    if MeV>0, % MeV style call
        Td = util.DriftPeriod('e-',MeV,alpha0_deg,L); % seconds
    else % -Td style call
        Td = -MeV;
    end
    omega_d = 2*pi./Td*24*60*60; % convert from seconds to rad/day    
    T = 0.75/24; % decay time, days from B&A
    B0 = 0.311; % Earth Dipole Moment, Gauss, from B&A
    Erms = 0.26*(Kp-1)+0.1; % rms E field mV/m
    
    % factor of 10 accounts for G and mV/m to T and V/m
    DLLE = 1/4*(10*Erms/B0/util.mks.c).^2.*(T./(1+(omega_d.*T/2).^2)).*L.^6; % from B&A
    % there's still some missing factor, because these results are 10^14
    % too high or so...?
else
    DLLE = nan(size(DLLM));
end
