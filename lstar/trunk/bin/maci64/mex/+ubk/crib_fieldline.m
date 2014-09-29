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
help ubk.fieldline

%% Ex 1. Field line info
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

d = datenum([2001 1 1 1 0 0]); % Date
m_threads = 10; % Number of threads to use
% parmod = [2, -10, 3, -10, 0 0 0 0 0 0]';
parmod = [2, A]';
external = 'TS07';
internal = 'IGRF';
ionoR = 1.0; % RE
ds = 0.05; % RE

[r p] = meshgrid(2:6, (0:45:359)*pi/180);
[xsm ysm zsm] = pol2cart(p,r,zeros(size(r)));
% zsm = rand(size(xsm))*2 - 1;

disp('%% Field line info with IGRF+TS07')
tic
[ K Bm Xmeq Ymeq Zmeq Bmeq Xeq Yeq Zeq Beq Xfoot Yfoot Zfoot ...
    Xfl Yfl Zfl] = ...
    ubk.fieldline( xsm(:), ysm(:), zsm(:), d, parmod, ...
    external, internal,...
    'ionor',ionoR, 'ds',ds, 'm_threads',m_threads);
K = reshape(K, size(xsm));
Bm = reshape(Bm, size(xsm));
Xmeq = reshape(Xmeq, size(xsm));
Ymeq = reshape(Ymeq, size(xsm));
Zmeq = reshape(Zmeq, size(xsm));
Bmeq = reshape(Bmeq, size(xsm));
Xeq = reshape(Xeq, size(xsm));
Yeq = reshape(Yeq, size(xsm));
Zeq = reshape(Zeq, size(xsm));
Beq = reshape(Beq, size(xsm));
Xfoot = reshape(Xfoot, size(xsm));
Yfoot = reshape(Yfoot, size(xsm));
Zfoot = reshape(Zfoot, size(xsm));
Xfl = reshape(Xfl, size(xsm));
Yfl = reshape(Yfl, size(xsm));
Yfl = reshape(Yfl, size(xsm));
toc

subplot(2,2,1)
Rfl = [Xfl(:) Yfl(:) Zfl(:)]';
Rfl(4,:) = {'-r'};
h1 = plot3(Rfl{:});
hold on
h2 = plot3(Xmeq, Ymeq, Zmeq,'xb'); % Magnetic equator
h3 = plot3(Xeq, Yeq, Zeq,'.g'); % Coordnate equator
hold off
title('Field lines')

subplot(2,2,2)
[xe ye ze] = sphere(20);
surf(xe,ye,ze,'facecolor','w')
hold on
plot3(Rfl{:});
plot3(Xfoot, Yfoot, Zfoot, '.b')
hold off
title('Foot points')
axis([-2 2 -2 2 -2 2])

subplot(2,2,[3 4])
a = [Bm(:) K(:)]';
semilogx(a{:})
xlabel('Bm [nT]')
ylabel('K')
title('K(Bm)')