function DLLM = odc_DLL_Lejosne2013(L,Kp,alpha0_deg,MeV)
% ********************************************
% These coefficients were published in Lejosne et al, 2013, JGR
% this routine is now a wrapper for that one
% ********************************************
% DLLM = odc_DLL_Lejosne2013(L,Kp)
% DLLM = odc_DLL_Lejosne2013(L,Kp,alpha0_deg)
% DLLM = odc_DLL_Lejosne2013(L,Kp,alpha0_deg,MeV)
% DLLM = odc_DLL_Lejosne2013(L,Kp,alpha0_deg,MeV,QDoption)
% compute radial diffusion coefficient
% according to Lejosne et al., 2013
% L - dipole L shell
% Kp : 0-4
% alpha0_deg - equatorial pitch angle, degrees (default is 90)
% MeV - particle energy, MeV
% DLLM - radial diffusion coefficient (1/day)
% includes angular dependence for dipole
%
% with no arguments, makes figure 8 from paper
%  requires odc_DLL_BA2000

if nargin == 0,
    DLM = [];
    test;
    return;
end

if nargin < 3,
    alpha0_deg = 90;
    QD2 = 1; % no angular correction
else
    y = sind(alpha0_deg);
    util = odc_util;
    y = sind(alpha0_deg);
    %(Q(y)/D(y))^2,  % angular dependence of DLL M, S&L, 3.13
    QD2 = (util.Q(y)./util.D(y)).^2/(util.Q(1)/util.D(1))^2;
end

util = odc_util;
Td = util.DriftPeriod('e-',MeV,alpha0_deg,L);
Omega = 2*pi./Td*60; % rad/minute

Kp = round(Kp);
if (Kp > 4),
    error('Model not defined for Kp>4');
end
iKp = 1+Kp;
% table 2
C1 = [1.3e-2,4.26e-2,8.48e-2,2.45e-1,6.78e-1]; % nT^2 RE^-2 min^-2
C2 = [2.28e-4,3.83e-4,8.31e-4,1.31e-3,1.36e-3];
lambda1 = 2; % 1/minutes (apparent typo in published version has min^-4)
lambda2 = 1.43e-2;           % 1/minutes
% equation 23
RE = 1;
k0 = util.dipoleB(1,0,0); % nT
DLLM = L.^10*RE^8/(2*k0)^2*(5/7)^2.*(lambda1*C1(iKp)./(lambda1^2+Omega.^2)+lambda2*C2(iKp)./(lambda2^2+Omega.^2)); % 1/minute
DLLM = DLLM*24*60; % 1/day
DLLM = DLLM.*QD2; % apply angular correction

function test
keV = 10.^(0:0.01:4);
MeV = keV/1e3;
L = 6.6;
alpha0_deg = 90;
figure;
colors = {'k','r','g','b','m'};
hba = nan(5,1);
hlj = nan(5,1);
for Kp=0:4,
    DLLM_BA = odc_DLL_BA2000(L,Kp);
    DLLM_BA = repmat(DLLM_BA,size(MeV));
    DLLM_LJ = odc_DLL_Lejosne2013(L,Kp,alpha0_deg,MeV);
    hba(Kp+1) = loglog(keV,DLLM_BA,'-','color',colors{1+Kp},'linew',2);
    hold on;
    hlj(Kp+1) = loglog(keV,DLLM_LJ,'--','color',colors{1+Kp},'linew',2);
end
grid on;
xlabel('keV');
ylabel('D_{LL} 1/day');
legend([hba;hlj(1)],'Kp=0 BA2K','Kp=1','Kp=2','Kp=3','Kp=4','Kp=0 L13','location','eo');
