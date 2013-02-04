function [DLLM,DLLE] = odc_DLL_Ozeke2013(L,param,value,alpha0_deg,MeV)
% ********************************************
% WARNING! This code uses prepublication values! WARNING!
% ********************************************
% DLLM = odc_DLL_Ozeke2013(L,param,value)
% DLLM = odc_DLL_Ozeke2013(L,param,value,alpha0_deg)
% [DLLM,DLLE] = odc_DLL_Ozeke2012(L,param,value,alpha0_deg,MeV)
% compute radial diffusion coefficient
% according to Ozeke et al. 2013 (prepublication)
% L - dipole L shell
% param - 'Kp' (at a later date, this might include an option for 'Vsw')
% value - Kp
% alpha0_deg - equatorial pitch angle, degrees (default is 90)
% MeV - particle energy, MeV
%
% DLLM - magnetic radial diffusion coefficient (1/day)
% DLLE - electric radial diffusion coefficient (1/day)
%
%
% run with no input arguments to generate test figures
%
% test requires odc_DLL_Ozeke2012

persistent warned
if isempty(warned) || ~warned,
    warning('%s uses pre-publication values of DLLs',mfilename);
    warned = true;
end

if nargin == 0,
    test;
    DLLM = [];
    DLLE = [];
    return
end

if nargin < 4,
    alpha0_deg = 90;
    QD2 = 1; % no angular correction
else
    y = sind(alpha0_deg);    
    % formula is for equatorially mirroring particles
    % we'll apply the appropriate angular correction per S&L
    % DLL_Kp = @(alpha0,E,L,Kp)DLL_90(Kp,L).*QD2(alpha0);
    util = odc_util;
    y = sind(alpha0_deg);
    %(Q(y)/D(y))^2,  % angular dependence of DLL M, S&L, 3.13
    QD2 = (util.Q(y)./util.D(y)).^2/(util.Q(1)/util.D(1))^2;
end

switch(lower(param)),
    case 'kp',
        Kp = value;
        DLLM = 1.94e-6*L.^8.*10.^(L.^2*0.112 + L.*Kp.*0.0824 - L.*1.32+Kp.^2.*0.0291-Kp.*0.270); % eqn 18
        DLLE = 5.75e-9*L.^6.*10.^(0.208.*L+0.457.*Kp); % eqn 21
    otherwise,
        error('Unknown parameter "%s"',param);
end

DLLM = DLLM.*QD2; % apply angular correction

function test
Ms = [500, 2000, 3500, 5000]; % MeV/G
Kps = [5,3,1];
L = 3:7;
alpha0_deg=90;
util = odc_util;

% DLL M plot
f = figure;
pos = get(f,'pos');
set(f,'pos',[pos(1)-pos(3)/2,pos(2)-pos(4),pos(3)*2,pos(4)*2]);
for iM = 1:4,
    M = Ms(iM);
    B = util.SL.B0*1e9./L.^3;
    MeV = util.MBtoMeV(M,B,alpha0_deg,'e-');
    for iKp = 1:3,
        Kp = Kps(iKp);
        subplot(3,4,4*(iKp-1)+iM);
        DLLMm1 = nan*L;
        DLLEm1 = nan*L;
        DLLMm10 = nan*L;
        DLLEm10 = nan*L;
        DLLM = nan*L;
        DLLE = nan*L;
        for iL = 1:length(L),
            [DLLMm1(iL),DLLEm1(iL)] = odc_DLL_Ozeke2012(L(iL),'Kp',Kp,alpha0_deg,MeV(iL),1);
            [DLLMm10(iL),DLLEm10(iL)] = odc_DLL_Ozeke2012(L(iL),'Kp',Kp,alpha0_deg,MeV(iL),10);
            [DLLM(iL),DLLE(iL)] = odc_DLL_Ozeke2013(L(iL),'Kp',Kp,alpha0_deg,MeV(iL));
        end
        loglog(L,DLLMm1,'s',L,DLLMm10,'o',L,DLLM,'-');
        axis([2 7 2e-6 2e-1]);
        set(gca,'xtick',3:7,'xticklabel',{'3','4','5','6','7'},'ytick',10.^(-5:-1));
        if (iKp==1) && (iM == 1),
            legend('m=1 ''12','m=10 ''12','2013','location','nw');
        end
        if iKp==1,
            title(sprintf('%g MeV/G',M));
        end
        if iM == 1,
            ylabel(sprintf('Kp=%g B_z D_{LL} /days',Kp));
        end
    end
end

% DLL E plot
f = figure;
pos = get(f,'pos');
set(f,'pos',[pos(1)-pos(3)/2,pos(2)-pos(4),pos(3)*2,pos(4)*2]);
for iM = 1:4,
    M = Ms(iM);
    B = util.SL.B0*1e9./L.^3;
    MeV = util.MBtoMeV(M,B,alpha0_deg,'e-');
    for iKp = 1:3,
        Kp = Kps(iKp);
        subplot(3,4,4*(iKp-1)+iM);
        DLLMm1 = nan*L;
        DLLEm1 = nan*L;
        DLLMm10 = nan*L;
        DLLEm10 = nan*L;
        DLLM = nan*L;
        DLLE = nan*L;
        for iL = 1:length(L),
            [DLLMm1(iL),DLLEm1(iL)] = odc_DLL_Ozeke2012(L(iL),'Kp',Kp,alpha0_deg,MeV(iL),1);
            [DLLMm10(iL),DLLEm10(iL)] = odc_DLL_Ozeke2012(L(iL),'Kp',Kp,alpha0_deg,MeV(iL),10);
            [DLLM(iL),DLLE(iL)] = odc_DLL_Ozeke2013(L(iL),'Kp',Kp,alpha0_deg,MeV(iL));
        end
        loglog(L,DLLEm1,'s',L,DLLEm10,'o',L,DLLE,'-');
        axis([2 7 2e-6 2e+1]);
        set(gca,'xtick',3:7,'xticklabel',{'3','4','5','6','7'},'ytick',10.^(-5:+1));
        if (iKp==1) && (iM == 1),
            legend('m=1 ''12','m=10 ''12','2013','location','nw');
        end
        if iKp==1,
            title(sprintf('%g MeV/G',M));
        end
        if iM == 1,
            ylabel(sprintf('Kp=%g E_\\phi D_{LL} /days',Kp));
        end
    end
end
