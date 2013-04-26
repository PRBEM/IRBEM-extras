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
help ubk.cotrans

%% Ex 1. GSM to SM to GSM
d = [now, now]; % Date
n_threads = 1; % Number of threads to use

xgsm = rand(100000,2)*2 - 1;
ygsm = rand(100000,2)*2 - 1;
zgsm = rand(100000,2)*2 - 1;

xgsm = xgsm * 10;
ygsm = ygsm * 10;
zgsm = zgsm * 10;

tic
[xsm ysm zsm] = ubk.cotrans(xgsm,ygsm,zgsm,d,'gsm2sm',n_threads);
[xgsm1 ygsm1 zgsm1] = ubk.cotrans(xsm,ysm,zsm,d,'sm2gsm',n_threads);
toc

dx = abs(xgsm - xgsm1);
dy = abs(ygsm - ygsm1);
dz = abs(zgsm - zgsm1);

fprintf('|dx| = %e, |dy| = %e, |dz| = %e\n', ...
    nanmean(dx(:)), nanmean(dy(:)), nanmean(dz(:)))