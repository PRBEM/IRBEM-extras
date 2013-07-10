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
d = datenum([2001 1 1 1 0 0]); % Date
m_threads = 10; % Number of threads to use
parmod = [2, -10, 3, -10, 0 0 0 0 0 0]';
external = 'TS05';
internal = 'IGRF';
ionoR = 1.0; % RE
ds = 0.05; % RE

[r p] = meshgrid(2:6, (0:45:359)*pi/180);
[xsm ysm zsm] = pol2cart(p,r,zeros(size(r)));
% zsm = rand(size(xsm))*2 - 1;

disp('%% Field line info with IGRF+TS05')
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