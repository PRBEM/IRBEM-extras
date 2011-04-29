% test odc_Daa_FA_local
%% reproduce some figures from Summers, 2005, and 2007
% as of 4/25/2011, visual comparison shows apparent success
% computing diffusion coefficients for hydrogen and multispecies plasmas

hemi = +1; % always in northern hemisphere
util = odc_util; % get constants and function handles

%% Figure 1. 2005
% R-mode waves with alpha_star = 0.16 and dB = 0.1 nT

keV = [100 300 1000 3000];
styles = {'r.-','g.-','b.-','k.-'};

DoubleHeightFigure;
alpha = (0.1:0.1:90)';
% dummy field values for realistic frequencies
L = 4.5; % RE
Beq = 30e3/L^3; % equator, nT
B = util.fce2B(9.54e3); % f_ce = 9.54 kHz
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
% odc_Daa_FA_local('e',3,85,L,MLT,B,Beq,hemi,wave_model);

h = nan(size(keV));
for iE = 1:length(keV),
    MeV = keV(iE)/1e3;
    Daa = nan(size(alpha));
    Dap = nan(size(alpha));
    Dpp = nan(size(alpha));
    for ialpha = 1:length(alpha),
        [Daa(ialpha),Dap(ialpha),Dpp(ialpha)] = odc_Daa_FA_local('e',MeV,alpha(ialpha),L,MLT,B,Beq,hemi,wave_model);
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

%% Figure 2. 2005
% R-mode waves with alpha_star = 0.16 and dB = 0.1 nT

alpha_star = [0.4444,0.16,0.04,0.0178,0.01];
MeV = 1;
styles = {'r.-','g.-','b.-','k.-','m.-'};

DoubleHeightFigure;
alpha = (0.1:0.1:90)';
% dummy field values for realistic frequencies
L = 4.5; % RE
Beq = 30e3/L^3; % equator, nT
B = util.fce2B(9.54e3); % f_ce = 9.54 kHz
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
% odc_Daa_FA_local('e',3,85,L,MLT,B,Beq,hemi,wave_model);

h = nan(size(alpha_star));
for istar = 1:length(alpha_star),
    wave_model.alpha_star = alpha_star(istar);
    Daa = nan(size(alpha));
    Dap = nan(size(alpha));
    Dpp = nan(size(alpha));
    for ialpha = 1:length(alpha),
        [Daa(ialpha),Dap(ialpha),Dpp(ialpha)] = odc_Daa_FA_local('e',MeV,alpha(ialpha),L,MLT,B,Beq,hemi,wave_model);
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

%% Figure 3. 2005
% R-mode waves with forward and backward chorus

keV = 1000;
alpha = (0.1:0.1:90)';
% dummy field values for realistic frequencies
L = 4; % RE
Beq = 30e3/L^3; % equator, nT
B = util.fce2B(13.65e3); % f_ce = 13.65 kHz
MLT = 6; % hours

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
        [Daa(ialpha),Dap(ialpha),Dpp(ialpha)] = odc_Daa_FA_local('e',MeV,alpha(ialpha),L,MLT,B,Beq,hemi,wave_model);
    end
    h(i) = semilogy(alpha,Dpp,styles{i});
    hold on;
end    
axis([0 90 1e-5 1e-3]);
ylabel('D_{pp}/p^2 1/sec');
legend(h,'backward only, k<0','forward and backward','location','nw');
grid on;

%% Figure 4. 2005
% L-mode waves electrons

keV = [1.25 1.5 2 5 10]*1e3;
styles = {'r-','g--','b:','k.-','r.:'};

figure;
alpha = (0.1:0.1:90)';
% dummy field values for realistic frequencies
L = 4; % RE
Beq = 30e3/L^3; % equator, nT
B = util.fce2B(13.65e3); % f_ce = 13.65 kHz
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
% odc_Daa_FA_local('e',3,85,L,MLT,B,Beq,hemi,wave_model);

h = nan(size(keV));
for iE = 1:length(keV),
    MeV = keV(iE)/1e3;
    Daa = nan(size(alpha));
    Dap = nan(size(alpha));
    Dpp = nan(size(alpha));
    for ialpha = 1:length(alpha),
        [Daa(ialpha),Dap(ialpha),Dpp(ialpha)] = odc_Daa_FA_local('e',MeV,alpha(ialpha),L,MLT,B,Beq,hemi,wave_model);
    end
    h(iE) = semilogy(alpha,Daa,styles{iE},'linew',2);
    hold on;
end

axis([0 90 1e-3 1e0]);
ylabel('D_{\alpha\alpha} 1/sec');
legend(h,arrayfun(@(x)sprintf('%g',x/1e3),keV,'uniform',false),'location','sw');
grid on;
xlabel('Local Pitch Angle, \alpha (deg)');

%% Figure 5. 2005
% L-mode waves protons

keV = [25 50 75 110 150];
styles = {'r-','g--','b:','k.-','r.:'};

figure;
alpha = (0.1:0.1:90)';
% dummy field values for realistic frequencies
L = 4; % RE
Beq = 30e3/L^3; % equator, nT
B = util.fce2B(13.65e3); % f_ce = 13.65 kHz
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
% odc_Daa_FA_local('e',3,85,L,MLT,B,Beq,hemi,wave_model);

h = nan(size(keV));
for iE = 1:length(keV),
    MeV = keV(iE)/1e3;
    Daa = nan(size(alpha));
    Dap = nan(size(alpha));
    Dpp = nan(size(alpha));
    for ialpha = 1:length(alpha),
        [Daa(ialpha),Dap(ialpha),Dpp(ialpha)] = odc_Daa_FA_local('H+',MeV,alpha(ialpha),L,MLT,B,Beq,hemi,wave_model);
    end
    h(iE) = semilogy(alpha,Daa,styles{iE},'linew',2);
    hold on;
end

axis([0 90 1e-5 1e-3]);
ylabel('D_{\alpha\alpha} 1/sec');
legend(h,arrayfun(@(x)sprintf('%g',x),keV,'uniform',false),'location','sw');
grid on;
xlabel('Local Pitch Angle, \alpha (deg)');

%% Figure 15. 2007
% L-mode waves electrons, H+ band

styles = {'r-','b-','k-','c-','r-'};

DoubleHeightFigure;
alpha = (0.1:0.1:90)';
% dummy field values for realistic frequencies
L = 4; % RE
Beq = 30e3/L^3; % equator, nT
B = Beq;
MLT = 6; % hours
drift_av_weighting = 0.01; % 1% weighting

clear wave_model
wave_model.mode = 'L';
wave_model.dB = 1; % nT
%wave_model.alpha_star = 1e-3; 
wave_model.normalization = 'Omega_p';
wave_model.omega_m = 0.6;
wave_model.domega = 0.1;
wave_model.omega1 = 0.5;
wave_model.omega2 = 0.7;
wave_model.directions = 'f'; % forward
wave_model.composition = [0.85 0.1 0.05]; % H+, He+, O+

alpha_star = [1e-3, 1e-2];
for ia = 1:2,    
    subplot(2,1,ia);
    switch(ia),
        case 1,
            keV = [0.5 1 2 5]*1e3;
        case 2,
            keV = [2 5 10]*1e3;
    end
    wave_model.alpha_star = alpha_star(ia);
    h = nan(size(keV));
    for iE = 1:length(keV),
        MeV = keV(iE)/1e3;
        Daa = nan(size(alpha));
        Dap = nan(size(alpha));
        Dpp = nan(size(alpha));
        for ialpha = 1:length(alpha),
            [Daa(ialpha),Dap(ialpha),Dpp(ialpha)] = odc_Daa_FA_local('e',MeV,alpha(ialpha),L,MLT,B,Beq,hemi,wave_model);
        end
        h(iE) = semilogy(alpha,Daa*drift_av_weighting,styles{iE},'linew',2);
        hold on;
    end    
    axis([0 90 1e-6 1e-2]);
    ylabel('D_{\alpha\alpha} 1/sec');
    legend(h,arrayfun(@(x)sprintf('%g',x/1e3),keV,'uniform',false),'location','sw');
    grid on;
    xlabel('Local Pitch Angle, \alpha (deg)');
    title(sprintf('\\alpha^*=%g',alpha_star(ia)));
end


%% Figure 16. 2007
% L-mode waves electrons, He+ band

styles = {'r-','b-','k-','c-','r-'};

DoubleHeightFigure;
alpha = (0.1:0.1:90)';
% dummy field values for realistic frequencies
L = 4; % RE
Beq = 30e3/L^3; % equator, nT
B = Beq;
MLT = 6; % hours
drift_av_weighting = 0.01; % 1% weighting

clear wave_model
wave_model.mode = 'L';
wave_model.dB = 1; % nT
%wave_model.alpha_star = 1e-3; 
wave_model.normalization = 'Omega_O+';
wave_model.omega_m = 3;
wave_model.domega = 0.5;
wave_model.omega1 = 2.5;
wave_model.omega2 = 3.5;
wave_model.directions = 'f'; % forward
wave_model.composition = [0.7 0.2 0.1]; % H+, He+, O+

alpha_star = [1e-3, 1e-2];
for ia = 1:2,    
    subplot(2,1,ia);
    switch(ia),
        case 1,
            keV = [1 2 5 10]*1e3;
        case 2,
            keV = [5 10]*1e3;
    end
    wave_model.alpha_star = alpha_star(ia);
    h = nan(size(keV));
    for iE = 1:length(keV),
        MeV = keV(iE)/1e3;
        Daa = nan(size(alpha));
        Dap = nan(size(alpha));
        Dpp = nan(size(alpha));
        for ialpha = 1:length(alpha),
            [Daa(ialpha),Dap(ialpha),Dpp(ialpha)] = odc_Daa_FA_local('e',MeV,alpha(ialpha),L,MLT,B,Beq,hemi,wave_model);
        end
        h(iE) = semilogy(alpha,Daa*drift_av_weighting,styles{iE},'linew',2);
        hold on;
    end    
    axis([0 90 1e-6 1e-2]);
    ylabel('D_{\alpha\alpha} 1/sec');
    legend(h,arrayfun(@(x)sprintf('%g',x/1e3),keV,'uniform',false),'location','sw');
    grid on;
    xlabel('Local Pitch Angle, \alpha (deg)');
    title(sprintf('\\alpha^*=%g',alpha_star(ia)));
end

%% Figure 17. 2007
% L-mode waves electrons, O+ band

styles = {'r-','b-','k-','c-','r-'};

DoubleHeightFigure;
alpha = (0.1:0.1:90)';
% dummy field values for realistic frequencies
L = 4; % RE
Beq = 30e3/L^3; % equator, nT
B = Beq;
MLT = 6; % hours
drift_av_weighting = 0.01; % 1% weighting

clear wave_model
wave_model.mode = 'L';
wave_model.dB = 1; % nT
%wave_model.alpha_star = 1e-3; 
wave_model.normalization = 'Omega_O+';
wave_model.omega_m = 0.9;
wave_model.domega = 0.05;
wave_model.omega1 = 0.85;
wave_model.omega2 = 0.95;
wave_model.directions = 'f'; % forward
wave_model.composition = [0.6 0.2 0.2]; % H+, He+, O+

alpha_star = [1e-3, 1e-2];
for ia = 1:2,    
    subplot(2,1,ia);
    switch(ia),
        case 1,
            keV = [2 5 10]*1e3;
        case 2,
            keV = [5 10]*1e3;
    end
    wave_model.alpha_star = alpha_star(ia);
    h = nan(size(keV));
    for iE = 1:length(keV),
        MeV = keV(iE)/1e3;
        Daa = nan(size(alpha));
        Dap = nan(size(alpha));
        Dpp = nan(size(alpha));
        for ialpha = 1:length(alpha),
            [Daa(ialpha),Dap(ialpha),Dpp(ialpha)] = odc_Daa_FA_local('e',MeV,alpha(ialpha),L,MLT,B,Beq,hemi,wave_model);
        end
        h(iE) = semilogy(alpha,Daa*drift_av_weighting,styles{iE},'linew',2);
        hold on;
    end    
    axis([0 90 1e-6 1e-2]);
    ylabel('D_{\alpha\alpha} 1/sec');
    legend(h,arrayfun(@(x)sprintf('%g',x/1e3),keV,'uniform',false),'location','sw');
    grid on;
    xlabel('Local Pitch Angle, \alpha (deg)');
    title(sprintf('\\alpha^*=%g',alpha_star(ia)));
end

