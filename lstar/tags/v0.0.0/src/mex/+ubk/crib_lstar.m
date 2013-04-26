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

%% Ex 1. Ldip v.s. UBKL
if 0
    d = now; % Date
    n_threads = 10; % Number of threads to use
    ioptparmod = 2;
    external = 'none';
    internal = 'dip';
    
    [rsm psm pa0] = meshgrid(2:.1:10, (0:10:359)*pi/180, (10:5:90)*pi/180);
    
    disp('%% L* with DIP+NONE')
    tic
    Ls = ...
        ubk.lstar( rsm(:), psm(:), pa0(:), d, ioptparmod, ...
        external, internal);
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
if 0
    d = now; % Date
    n_threads = 10; % Number of threads to use
    parmod = [4, -100, 2, -15, 1, 1, 1, 1, 1, 1]';
    external = 'ts05';
    internal = 'igrf';
    ionoR = 1.01; % RE
    
    pa0 = (10:10:90)*pi/180;
    xsm = -6.6 * ones(size(pa0));
    ysm = 0 * zeros(size(pa0));
    
    disp('%% Drift contour')
    tic
    [~, ~, ~, Xc, Yc, ~, ~] = ...
        ubk.lstar( xsm(:), ysm(:), pa0(:), d, parmod, ...
        external, internal, ionoR, [], [], [], ...
        [], [], true, true, n_threads);
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
if 1
    parmod = [4, -100, 2, -15, 1, 1, 1, 1, 1, 1];
    external = 't89';
    internal = 'igrf';
    
    T = datenum([2007 3 23]) + linspace(0,1,1441);
    [psm pa0] = meshgrid(linspace(0,2*pi,length(T)), linspace(10,90,8)*pi/180);
    rsm = 6.6 * ones(size(psm));
    parmod = 3*ones(size(T)); %repmat(parmod', 1, length(T));
    
    disp('%% L* tracking')
    tic
    Ls = ...
        ubk.lstar( rsm, psm, pa0, ...
        T, parmod, external, internal);
    toc
    
    figure(3)
    subplot(1,1,1)
    plot(T,Ls)
    datetick
end
