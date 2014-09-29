%
% $Author$
% $LastChangedDate$
% $Revision$
% $Id$
%

clc
clear
clf
set(0,'defaultaxesfontname','times')
set(0,'defaulttextfontname','times')
set(0,'defaulttextfontsize',12)
set(0,'defaultaxesfontsize',12)

%% Help
help ubk.lstar

USE_INTERP_FIELD_WITH_POLY_ORDER = 2;
A = [1.00000,-17.2196,-43.8456,-73.3322,-97.5089,-119.160,44.1685,-28.3264,...
    -10.0919,-6.92003,-9.75857,15.2687,-34.5620,-4.25799,-24.7348,-2.83235,...
    -0.496528,-24.5300,-24.1017,17.9590,18.7953,10.5520,11.1135,-15.2804,...
    -0.598420,4.85943,-10.3045,-6.15506,13.5167,-3.53762,-5.13419,5.62481,...
    11.6524,4.45178,14.5666,-19.2603,-11.6760,5.27181,16.2248,-4.03140,...
    -10.1715,-9.40554,-8.63108,19.2596,8.45173,5.43871,-7.09198,-5.21992,...
    5.06120,8.87807,15.7564,46.7266,-48.0707,-8.33984,1.25046,-23.5229,...
    0.877495,-31.7131,4.87073,-8.23912,1.57636,-33.3436,-21.1656,-3.80171,...
    -14.3754,6.42494,-2.94579,-18.7112,37.4036,10.8904,-1.33068,6.50348,...
    1.51085,-6.63102,4.40149,-2.77839,-3.40469,-4.87805,17.3867,-7.80144,...
    0.872749,9.48005,5.42796,-2.30417,17.9994,0.485395,-8.85228,1.79156,...
    -9.94661,0.396798E-01,7.30606,0.421574,0.575977E-01,-0.246429,...
    0.678498E-01,2.63071,9.15388,32.1001,0.820630,1.71374,0.141481E-01];

%% Ex 1. Ldip v.s. UBKL
if 0
    d = now; % Date
    ioptparmod = 2;
    external = 'none';
    internal = 'dip';
    
    [rsm psm pa0] = meshgrid(2:.1:10, (0:10:359)*pi/180, (10:5:90)*pi/180);
    
    disp('%% L* with DIP+NONE')
    tic
    Ls = ...
        ubk.lstar( rsm(:), psm(:), pa0(:), d, ioptparmod, ...
        external, internal, 'USE_INTERP_FIELD_WITH_POLY_ORDER', USE_INTERP_FIELD_WITH_POLY_ORDER);
    Ls = reshape(Ls, size(rsm));
    toc
    
    Ls = permute(Ls, [2 1 3]);
    Ls = reshape(Ls, size(rsm,2),size(rsm,1)*size(rsm,3));
    Ldip = repmat(rsm(1,:,1)',1,size(Ls,2));
    dL1 = abs(Ls - Ldip)./Ldip;
    
    figure(1)
    subplot(1,1,1)
    h1 = errorbar(Ldip(:,1),nanmean(abs(dL1),2),nanstd(abs(dL1),[],2),'r');
    xlabel('Ldip')
    ylabel('|dL*|')
end

%% Ex 2. Orbit calculation
if 1
    d = now; % Date
    m_threads = 10; % Number of threads to use
    %parmod = [4, -100, 2, -15, 1, 1, 1, 1, 1, 1]';
    parmod = [4, A]';
    external = 'ts07';
    internal = 'igrf';
    ionoR = 1.01; % RE
    
    pa0 = (10:10:90)*pi/180;
    xsm = -6.6 * ones(size(pa0));
    ysm = 0 * zeros(size(pa0));
    
    disp('%% Drift contour')
    tic
    [~, ~, ~, Xc, Yc, ~, ~] = ...
        ubk.lstar( xsm(:), ysm(:), pa0(:), d, parmod, ...
        external, internal, ...
        'ionor',ionoR, 'isCARTESIANGRID',true, ...
        'm_threads',m_threads, 'USE_INTERP_FIELD_WITH_POLY_ORDER', USE_INTERP_FIELD_WITH_POLY_ORDER);
    Xc = reshape(Xc, size(xsm));
    Yc = reshape(Yc, size(xsm));
    toc
    
    figure(2)
    co = jet(length(pa0));
    subplot(1,1,1, 'box','on','colororder',co)
    plot(xsm,ysm,'ok','markerfacecolor','k')
    hold on
    Rc = [Xc(:) Yc(:)]';
    for idx=1:length(pa0)
        plot(Rc{:,idx},'color',co(idx,:))
    end
    hold off
    xlabel('x sm [RE]')
    ylabel('y sm [RE]')
    caxis(minmax(pa0(:)')*180/pi)
    ylabel(colorbar,'pa [deg]')
    axis equal
end

%% Ex 3. Satellite tracking
% Let's assume GEO orbiting satellite.
% If you are patient, try!
if 0
    parmod = [4, -100, 2, -15, 1, 1, 1, 1, 1, 1];
    external = 'none';
    internal = 'dip';
    
    T = datenum([2007 3 23]) + linspace(0,1,1441);
    [psm, pa0] = meshgrid(linspace(0,2*pi,length(T)), linspace(10,90,8)*pi/180);
    rsm = 6.6 * ones(size(psm));
    parmod = 3*ones(size(T)); %repmat(parmod', 1, length(T));
    
    disp('%% L* tracking')
    tic
    Ls = ...
        ubk.lstar( rsm, psm, pa0, ...
        T, parmod, external, internal,...
        'n_threads',4,'m_threads',4, 'USE_INTERP_FIELD_WITH_POLY_ORDER', USE_INTERP_FIELD_WITH_POLY_ORDER);
    toc
    
    figure(3)
    subplot(1,1,1)
    plot(T,Ls)
    datetick
end
