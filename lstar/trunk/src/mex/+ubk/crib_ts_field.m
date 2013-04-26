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
