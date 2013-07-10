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
help ubk.ts_field

%% Ex 1. T89 field
d = now; % Date
n_threads = 10; % Number of threads to use
iopt= 3; % Kp = 2
external = 't89';
internal = 'igrf';

[x y z] = meshgrid(-10:.05:10, -10:.05:10, 0);

disp('%% T89 field')
tic
[bx by bz] = ubk.ts_field(x(:),y(:),z(:),d,iopt,external,internal,...
    'co_system','gsm','m_threads',n_threads);
bx = reshape(bx, size(x));
by = reshape(by, size(x));
bz = reshape(bz, size(x));
toc

figure(1)
bs = sqrt(bx.^2 + by.^2 + bz.^2);
imagesc(x(1,[1 end]),y([1 end],1),log10(bs))
caxis([0 3])
xlabel('xgsm [RE]')
ylabel('xgsm [RE]')
ylabel(colorbar,'log10 |B| [nT]')
title('T89')

%% Ex 2. TS05 field
% parmod=[PDYN (nPa), DST (nT), BYIMF, BZIMF (nT), W1 W2 W3 W4 W5 W6]
d = [now now-1000];
parmod = repmat([4, -100, 2, -15, 1, 1, 1, 1, 1, 1]', 1,2);
external = 'ts05';
internal = 'igrf';

x = -[6 6];
y = [0 0];
z = [0 0];

disp('%% TS05 field')
tic
[bx by bz] = ubk.ts_field(x,y,z,d,parmod,external,internal,...
    'co_system','gsm','m_threads',n_threads)
toc

%% Ex 2. TS07 field
% parmod=[PDYN (nPa), A (101)]
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
d = now; d = [d d d];
parmod = repmat([4, A]', 1,3);
external = 'ts07';
internal = 'igrf';

x = -[6 6 6];
y = [0 0 0];
z = [0 0 0];

disp('%% TS07 field')
tic
[bx by bz] = ubk.ts_field(x,y,z,d,parmod,external,internal,...
    'co_system','gsm','m_threads',n_threads)
toc
