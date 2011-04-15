% test Daa_FA_local



%% reproduce some figures from Summers, 2005

hemi = +1; % always in northern hemisphere

%% Figure 1.
% R-mode waves with alpha_star = 0.16 and dB = 0.1 nT

keV = [100 300 1000 3000];
styles = {'r.-','g.-','b.-','k.-'};

DoubleHeightFigure;
alpha = (0.1:0.1:90)';
% dummy field values for realistic frequencies
L = 4.5; % RE
Beq = 30e3/L^3; % equator, nT
% f_ce = 9.54 kHz = Omega_e/2/pi
% Omega_e = qB/m
% B= Omega_e*m/q
m0c2 = 0.510998910; % rest mass, MeV
m0kg = m0c2*1.782661758877380e-030; % rest mass, kg
qC = -1.602176487e-19; % charge, C
B = 9.54e3*2*pi*m0kg/abs(qC)*1e9; % somewhere off the equator
MLT = 6; % hours

clear wave_model
wave_model.mode = 'R';
wave_model.dB = 0.1; % nT
wave_model.alpha_star = 0.16;
wave_model.normalization = 'Omega_e';
wave_model.omega_m = 0.35;
wave_model.domega = 0.15;
wave_model.sigma = 2;
wave_model.directions = 'b'; % backward waves only

% % print diagnostics
% fprintf('3 MeV, alpha=85\n');
% Daa_FA_local('e',3,85,L,MLT,B,Beq,hemi,wave_model);

h = nan(size(keV));
for iE = 1:length(keV),
    MeV = keV(iE)/1e3;
    Daa = nan(size(alpha));
    Dap = nan(size(alpha));
    Dpp = nan(size(alpha));
    for ialpha = 1:length(alpha),
        [Daa(ialpha),Dap(ialpha),Dpp(ialpha)] = Daa_FA_local('e',MeV,alpha(ialpha),L,MLT,B,Beq,hemi,wave_model);
    end
    subplot(3,1,1);
    h(iE) = semilogy(alpha,Daa,styles{iE});
    hold on;
    subplot(3,1,2);
    semilogy(alpha,abs(Dap),styles{iE});
    hold on;
    subplot(3,1,3);
    semilogy(alpha,Dpp,styles{iE});
    hold on;    
end

subplot(3,1,1);
axis([0 90 1e-8 1e-1]);
ylabel('D_{\alpha\alpha} 1/sec');
legend(h,arrayfun(@(x)sprintf('%g',x),keV,'uniform',false),'location','sw');
grid on;

subplot(3,1,2);
axis([0 90 1e-8 1e-1]);
ylabel('|D_{\alphap}|/p 1/sec');
grid on;

subplot(3,1,3);
axis([0 90 1e-8 1e-2]);
ylabel('D_{pp}/p^2 1/sec');
xlabel('Local Pitch Angle, \alpha (deg)');
grid on;


%% Figure 2.
% R-mode waves with alpha_star = 0.16 and dB = 0.1 nT

alpha_star = [0.4444,0.16,0.04,0.0178,0.01];
MeV = 1;
styles = {'r.-','g.-','b.-','k.-','m.-'};

DoubleHeightFigure;
alpha = (0.1:0.1:90)';
% dummy field values for realistic frequencies
L = 4.5; % RE
Beq = 30e3/L^3; % equator, nT
% f_ce = 9.54 kHz = Omega_e/2/pi
% Omega_e = qB/m
% B= Omega_e*m/q
m0c2 = 0.510998910; % rest mass, MeV
m0kg = m0c2*1.782661758877380e-030; % rest mass, kg
qC = -1.602176487e-19; % charge, C
B = 9.54e3*2*pi*m0kg/abs(qC)*1e9; % somewhere off the equator
MLT = 6; % hours

clear wave_model
wave_model.mode = 'R';
wave_model.dB = 0.1; % nT
wave_model.normalization = 'Omega_e';
wave_model.omega_m = 0.35;
wave_model.domega = 0.15;
wave_model.sigma = 2;
wave_model.directions = 'b'; % backward waves only

% % print diagnostics
% fprintf('3 MeV, alpha=85\n');
% Daa_FA_local('e',3,85,L,MLT,B,Beq,hemi,wave_model);

h = nan(size(alpha_star));
for istar = 1:length(alpha_star),
    wave_model.alpha_star = alpha_star(istar);
    Daa = nan(size(alpha));
    Dap = nan(size(alpha));
    Dpp = nan(size(alpha));
    for ialpha = 1:length(alpha),
        [Daa(ialpha),Dap(ialpha),Dpp(ialpha)] = Daa_FA_local('e',MeV,alpha(ialpha),L,MLT,B,Beq,hemi,wave_model);
    end
    subplot(3,1,1);
    h(istar) = semilogy(alpha,Daa,styles{istar});
    hold on;
    subplot(3,1,2);
    semilogy(alpha,abs(Dap),styles{istar});
    hold on;
    subplot(3,1,3);
    semilogy(alpha,Dpp,styles{istar});
    hold on;    
end

subplot(3,1,1);
axis([0 90 1e-8 1e-1]);
ylabel('D_{\alpha\alpha} 1/sec');
legend(h,arrayfun(@(x)sprintf('%g',x),alpha_star,'uniform',false),'location','sw');
grid on;

subplot(3,1,2);
axis([0 90 1e-8 1e-1]);
ylabel('|D_{\alphap}|/p 1/sec');
grid on;

subplot(3,1,3);
axis([0 90 1e-8 1e-2]);
ylabel('D_{pp}/p^2 1/sec');
xlabel('Local Pitch Angle, \alpha (deg)');
grid on;

%% Figure 3.
% R-mode waves with forward and backward chorus

keV = 1000;
alpha = (0.1:0.1:90)';
% dummy field values for realistic frequencies
L = 4; % RE
Beq = 30e3/L^3; % equator, nT
% f_ce = 9.54 kHz = Omega_e/2/pi
% Omega_e = qB/m
% B= Omega_e*m/q
m0c2 = 0.510998910; % rest mass, MeV
m0kg = m0c2*1.782661758877380e-030; % rest mass, kg
qC = -1.602176487e-19; % charge, C
B = 13.65e3*2*pi*m0kg/abs(qC)*1e9; % somewhere off the equator
MLT = 6; % hours
c = 2.99792458e8; % m/s


clear wave_model
wave_model.mode = 'R';
wave_model.dB = 0.1; % nT
wave_model.alpha_star = 1;
wave_model.normalization = 'Omega_e';
wave_model.omega_m = 0.3;
wave_model.domega = 0.2;
wave_model.sigma = 1;

h = nan(2,1);
styles = {'kx','k-'};
directions = {'b','bf'};
figure;
for i = 1:length(directions),
    wave_model.directions = directions{i};
    MeV = keV/1e3;
    Daa = nan(size(alpha));
    Dap = nan(size(alpha));
    Dpp = nan(size(alpha));
    for ialpha = 1:length(alpha),
        [Daa(ialpha),Dap(ialpha),Dpp(ialpha)] = Daa_FA_local('e',MeV,alpha(ialpha),L,MLT,B,Beq,hemi,wave_model);
    end
    h(i) = semilogy(alpha,Dpp,styles{i});
    hold on;
end    
axis([0 90 1e-5 1e-3]);
ylabel('D_{pp}/p^2 1/sec');
legend(h,'backward only, k<0','forward and backward','location','nw');
grid on;

%% Figure 4.
% L-mode waves electrons

keV = [1.25 1.5 2 5 10]*1e3;
styles = {'r-','g--','b:','k.-','r.:'};

figure;
alpha = (0.1:0.1:90)';
% dummy field values for realistic frequencies
L = 4; % RE
Beq = 30e3/L^3; % equator, nT
% f_ce = 9.54 kHz = Omega_e/2/pi
% Omega_e = qB/m
% B= Omega_e*m/q
m0c2 = 0.510998910; % rest mass, MeV
m0kg = m0c2*1.782661758877380e-030; % rest mass, kg
qC = -1.602176487e-19; % charge, C
B = 13.65e3*2*pi*m0kg/abs(qC)*1e9; % somewhere off the equator
MLT = 6; % hours

clear wave_model
wave_model.mode = 'L';
wave_model.dB = 1; % nT
wave_model.alpha_star = 2.3e-3;
wave_model.normalization = 'Omega_p';
wave_model.omega_m = 1/3;
wave_model.domega = 1/6;
wave_model.sigma = 1;
wave_model.directions = 'f'; % forwards only (backwards does not resonate)

% % print diagnostics
% fprintf('3 MeV, alpha=85\n');
% Daa_FA_local('e',3,85,L,MLT,B,Beq,hemi,wave_model);

h = nan(size(keV));
for iE = 1:length(keV),
    MeV = keV(iE)/1e3;
    Daa = nan(size(alpha));
    Dap = nan(size(alpha));
    Dpp = nan(size(alpha));
    for ialpha = 1:length(alpha),
        [Daa(ialpha),Dap(ialpha),Dpp(ialpha)] = Daa_FA_local('e',MeV,alpha(ialpha),L,MLT,B,Beq,hemi,wave_model);
    end
    h(iE) = semilogy(alpha,Daa,styles{iE},'linew',2);
    hold on;
end

axis([0 90 1e-3 1e0]);
ylabel('D_{\alpha\alpha} 1/sec');
legend(h,arrayfun(@(x)sprintf('%g',x/1e3),keV,'uniform',false),'location','sw');
grid on;
xlabel('Local Pitch Angle, \alpha (deg)');


%% Figure 5.
% L-mode waves protons

keV = [25 50 75 110 150];
styles = {'r-','g--','b:','k.-','r.:'};

figure;
alpha = (0.1:0.1:90)';
% dummy field values for realistic frequencies
L = 4; % RE
Beq = 30e3/L^3; % equator, nT
% f_ce = 9.54 kHz = Omega_e/2/pi
% Omega_e = qB/m
% B= Omega_e*m/q
m0c2 = 0.510998910; % electron rest mass, MeV
m0kg = m0c2*1.782661758877380e-030; % rest mass, kg
qC = -1.602176487e-19; % charge, C
B = 13.65e3*2*pi*m0kg/abs(qC)*1e9; % somewhere off the equator
MLT = 6; % hours

clear wave_model
wave_model.mode = 'L';
wave_model.dB = 1; % nT
wave_model.alpha_star = 4.6e-3;
wave_model.normalization = 'Omega_p';
wave_model.omega_m = 0.15;
wave_model.domega = 0.0375;
wave_model.sigma = 4/3;
wave_model.directions = 'bf'; % forwards only (backwards does not resonate)

% % print diagnostics
% fprintf('3 MeV, alpha=85\n');
% Daa_FA_local('e',3,85,L,MLT,B,Beq,hemi,wave_model);

h = nan(size(keV));
for iE = 1:length(keV),
    MeV = keV(iE)/1e3;
    Daa = nan(size(alpha));
    Dap = nan(size(alpha));
    Dpp = nan(size(alpha));
    for ialpha = 1:length(alpha),
        [Daa(ialpha),Dap(ialpha),Dpp(ialpha)] = Daa_FA_local('H+',MeV,alpha(ialpha),L,MLT,B,Beq,hemi,wave_model);
    end
    h(iE) = semilogy(alpha,Daa,styles{iE},'linew',2);
    hold on;
end

axis([0 90 1e-5 1e-3]);
ylabel('D_{\alpha\alpha} 1/sec');
legend(h,arrayfun(@(x)sprintf('%g',x),keV,'uniform',false),'location','sw');
grid on;
xlabel('Local Pitch Angle, \alpha (deg)');


% %% test against UCLA implementation
% 
% stop; % deprecated code hereafter
% 
% % dummy field values for realistic frequencies
% L = 3.5; % RE
% Beq = 31e3/L^3; % equator, nT
% B = Beq;
% MLT = 6; % hours
% 
% clear wave_model
% wave_model.mode = 'R';
% wave_model.dB = 0.1; % nT
% wave_model.f= 2.5;
% wave_model.omega_pe_normalization = 'Omega_e_eq';
% wave_model.omega_pe = 2.5;
% wave_model.normalization = 'Omega_e';
% wave_model.omega_m = 0.35;
% wave_model.domega = 0.15;
% wave_model.sigma = 2;
% wave_model.directions = 'b'; % backward waves only
% 
% constants_Daa; % set UCLA up run
% % f=2.5;
% % d_omega_per= 0.15;
% % omega_m_per=0.35;
% % o_uc=0.65;
% % o_lc=0.05;
% % dB=1.e-6;% G
% % lam_max=15;
% % L=3.5;
% % EMEV=1.;
% %  function y=Daa_local(alpha_eq, E, lambda, L)
% keV = 3000;
% maglat = 0;
% alpha = 1:89;
% for ialpha = 1:length(alpha),
%     Daa_UCLA(ialpha) =Daa_local(alpha(ialpha)*pi/180, keV/511, maglat, L);
%     Daa(ialpha) = Daa_FA_local('e',keV/1e3,alpha(ialpha),L,MLT,B,Beq,hemi,wave_model);
% end
% semilogy(alpha,Daa_UCLA,'k-',alpha,Daa,'r-');
% 
% 
