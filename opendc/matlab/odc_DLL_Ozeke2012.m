function [DLLB,DLLE] = odc_DLL_Ozeke2012(L,param,value,alpha0_deg,MeV,m)
% DLLB = odc_DLL_Ozeke2012(L,param,value)
% DLLB = odc_DLL_Ozeke2012(L,param,value,alpha0_deg)
% [DLLB,DLLE] = odc_DLL_Ozeke2012(L,param,value,alpha0_deg,MeV,m)
% compute radial diffusion coefficient
% according to Ozeke et al. 2012
% L - dipole L shell
% param - 'Kp' or 'Vsw'
% value - Kp or Vsw (in km/s)
% alpha0_deg - equatorial pitch angle, degrees (default is 90)
% MeV - particle energy, MeV
% m - assumed m number (default is 1)
%
% DLLB - magnetic radial diffusion coefficient (1/day)
% DLLE - electric radial diffusion coefficient (1/day)
%
%
% Ozeke, L. G., I. R. Mann, K. R. Murphy, I. J. Rae,
% D. K. Milling, S. R. Elkington, A. A. A. Chan, and H. J. Singer (2012),
% ULF Wave Derived Radiation Belt Radial Diffusion Coefficients,
% J. Geophys. Res., doi:10.1029/2011JA017463, in press.
% (accepted 23 February 2012)
%
% run with no input arguments to generate test figures from paper
%
% requires odc_PSD_Ozeke2012 and its data tables
% test also requires odc_DLL_BA2000

if nargin == 0,
    test;
    DLLB = [];
    DLLE = [];
    return
end

util = odc_util;

if nargin < 4,
    alpha0_deg = 90;
end
y = sind(alpha0_deg);

if nargin < 6,
    m = 1; % default m value
end

PSD = cell(1,nargout);

Td = util.DriftPeriod('e-',MeV,alpha0_deg,L); % seconds
BE = util.SL.B0*1e9; % equatorial field at Earth, nT
[ele.gamma,ele.v,ele.m] = util.MeVtogamma(MeV,'e-');
f_d = 1e3./Td; % mHz
f = m*f_d;
[PSD{:}] = odc_PSD_Ozeke2012(L,param,value,f);
p_perp = ele.m.*ele.v.*y; % kg m/s
M = p_perp.^2.*L.^3./(2.*util.mks.electron.m0.*BE); % equation (4)
% M ~ (kg m /s )^2 / (kg * nT) = kg m^2 / s^2 / nT

PSDM = PSD{1};
DLLB = M.^2./(8.*util.mks.electron.q.^2.*ele.gamma.^2.*BE^2.*util.mks.R_E^4).*L.^4.*m.^2.*PSDM; % equation (3)
% DLLB ~ (kg m^2 / s^2 / nT)^2 / (C^2 nT^2 m^4) * (nT^2/mHz)
% DLLB ~  nT^2 kg ^2 m^4 /s^4 / nT^2 / C^2 / nT^2 / m^4 / mHz
% DLLB ~  kg ^2 /s^4 / nT^2 / C^ 2 / mHz
% T = kg / C / s
% DLLB ~  (1/1e-21) s kg ^2 /s^4 / (kg/C/s)^2 / C^2
% DLLB ~  (1/1e-21) s C^2 s^2 kg^2 /s^4 / kg^2 / C^2
% DLLB ~  (1/1e-21) / s
DLLB = DLLB * 1e+21; % /s
DLLB = DLLB * 24*60*60; % /day

if nargout >=2,
    PSDE = PSD{2};
    DLLE = 1./(8.*BE^2.*util.mks.R_E^2).*L.^6.*PSDE; % equation (2)
    % DLLE ~ 1/(nT^2 m^2) * (mV/m)^2/mHz =
    % DLLE ~ (1e-3 V)^2 / m^2 / (1e-3 / s) / m^2 * m^4 / (1e-9 V)^2 / s ^2
    % DLLE ~ (1e-6 / 1e-3 / 1e-18) / s
    DLLE = DLLE * 1e15; % /s
    DLLE = DLLE * 24*60*60; % /day
else
    DLLE = nan(size(DLLB));
end


function test
% make selected figures from paper
util = odc_util;

% figure 9
figure(9);
Ls = [3:6,6.6]';
Kps = [0,2,6];
ms = [1 10];
Ms = [500 5000]; % MeV/G
styles = {'ro','bs'};
pl_styles = {'r-','b--'};
leg = cell(size(ms));
h = nan(length(ms),1); % scatter handle
alpha0_deg = 90;
for iM = 1:length(Ms),
    M = Ms(iM);
    for iKp = 1:length(Kps),
        Kp = Kps(iKp);
        DLLB = nan(length(Ls),length(ms));
        for im = 1:length(ms),
            m = ms(im);
            for iL = 1:length(Ls),
                L = Ls(iL);
                % compute MeV from M and B
                B = util.SL.B0*1e9/L^3;
                MeV = util.MBtoMeV(M,B,alpha0_deg,'e-');
                DLLB(iL,im) = odc_DLL_Ozeke2012(L,'Kp',Kp,alpha0_deg,MeV,m);
            end
            
            subplot(length(Kps),length(Ms),(iKp-1)*length(Ms)+iM);
            h(im) = loglog(Ls,DLLB(:,im),styles{im});
            hold on;
            P = polyfit(log(Ls),log(DLLB(:,im)),1);
            loglog(Ls,exp(polyval(P,log(Ls))),pl_styles{im});
            
            leg{im} = sprintf('m=%g',m);
        end
        if iKp==1,
            title(sprintf('M=%g MeV/g',M));
        end
        if iKp==length(Kps),
            xlabel('L Shell');
        end
        if (iKp==1) && (iM==1),
            legend(h,leg{:},'location','nw');
        end        
        if iM==1,
            ylabel(sprintf('D_{LL}^{M}, 1/day, Kp=%g',Kp));
        end
        axis([2 8 1e-7 1e1]);
        set(gca,'ytick',10.^(-6:0),'xtick',2:8);
    end
end

% figure 10
figure(10);
Ls = [2.55,2.98,4.21,4.26,5.4,6.51,7.94]';
Kps = [0,2,6];
ms = [1 10];
Ms = [500 5000]; % MeV/G
styles = {'ro','bs'};
pl_styles = {'r-','b--'};
leg = cell(size(ms));
h = nan(length(ms),1); % scatter handles
alpha0_deg = 90;
for iM = 1:length(Ms),
    M = Ms(iM);
    for iKp = 1:length(Kps),
        Kp = Kps(iKp);
        DLLB = nan(length(Ls),length(ms));
        DLLE = nan(length(Ls),length(ms));
        for im = 1:length(ms),
            m = ms(im);
            for iL = 1:length(Ls),
                L = Ls(iL);
                % compute MeV from M and B
                B = util.SL.B0*1e9/L^3;
                MeV = util.MBtoMeV(M,B,alpha0_deg,'e-');
                [DLLB(iL,im),DLLE(iL,im)] = odc_DLL_Ozeke2012(L,'Kp',Kp,alpha0_deg,MeV,m);
            end
            
            figure(10);
            subplot(length(Kps),length(Ms),(iKp-1)*length(Ms)+iM);
            h(im) = loglog(Ls,DLLE(:,im),styles{im});
            hold on;
            P = polyfit(log(Ls),log(DLLE(:,im)),1);
            loglog(Ls,exp(polyval(P,log(Ls))),pl_styles{im});
            
            leg{im} = sprintf('m=%g',m);
        end
        if iKp==1,
            title(sprintf('M=%g MeV/g',M));
        end
        if iKp==length(Kps),
            xlabel('L Shell');
        end
        if (iKp==1) && (iM==1),
            legend(h,leg{:},'location','nw');
        end        
        if iM==1,
            ylabel(sprintf('D_{LL}^{E}, 1/day, Kp=%g',Kp));
        end
        axis([2 8 1e-6 1e2]);
        set(gca,'ytick',10.^(-5:1),'xtick',2:8);
    end
end

% figure 11
figure(11);
Ms = [100 500 1000 5000];
Ls = 3:0.1:7;
Kps = [6 1];
m=1;
for iM = 1:length(Ms),
    M = Ms(iM);
    for iKp = 1:length(Kps),
        Kp = Kps(iKp);
        
        DLLE_BA = nan(length(Ls),1);
        DLLM_BA = nan(length(Ls),1);
        DLLE = nan(length(Ls),1);
        DLLB = nan(length(Ls),1);
        
        for iL = 1:length(Ls),
            L = Ls(iL);
            % compute MeV from M and B
            B = util.SL.B0*1e9/L^3;
            MeV = util.MBtoMeV(M,B,alpha0_deg,'e-');
            [DLLM_BA(iL),DLLE_BA(iL)] = odc_DLL_BA2000(L,Kp,alpha0_deg,MeV);
            [DLLB(iL),DLLE(iL)] = odc_DLL_Ozeke2012(L,'Kp',Kp,alpha0_deg,MeV,m);
        end
        
        subplot(length(Kps),length(Ms),(iKp-1)*length(Ms)+iM);
        loglog(Ls,DLLE_BA,'k--');
        hold on;
        loglog(Ls,DLLM_BA,'k:');
        loglog(Ls,DLLE,'rs-');
        loglog(Ls,DLLB,'g-.','linew',2);
        if iKp==1,
            title(sprintf('M=%g MeV/g',M));
        end
        if iKp==length(Kps),
            xlabel('L Shell');
        end
        if (iKp==1) && (iM==1),
            legend('B&A (E)','B&A (M)','Mapped E','Space B','location','nw');
        end
        if iM==1,
            ylabel(sprintf('D_{LL}^{M}, 1/day, Kp=%g',Kp));
        end
        if iKp==1,
            axis([3 7 1e-4 2e2]);
            set(gca,'ytick',10.^(-4:2),'xtick',3:7);
        else
            axis([3 7 1e-5 2e0]);
            set(gca,'ytick',10.^(-4:0),'xtick',3:7);
        end
    end
end
