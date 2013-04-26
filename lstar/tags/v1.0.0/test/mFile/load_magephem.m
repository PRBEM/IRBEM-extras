function [ mag_struct ] = load_magephem( magpath )
%LOAD_MAGEPHEM Summary of this function goes here
%   Detailed explanation goes here

%% Date/time
Date = hdf5read(magpath, '/Date'); % Integer of YYYYMMDD
UTC = hdf5read(magpath,'/UTC'); % hours from the start of the day
YYYY = floor(Date/10000); Date = Date - YYYY*10000;
MM = floor(Date/100);
DD = Date - MM*100;

Tday = datenum( double([YYYY MM DD]) ) + UTC/24;

%% Data
Alpha = hdf5read(magpath, '/Alpha');
Alen = length(Alpha);
PArad = repmat(Alpha*pi/180, 1, length(Tday));

Rsm = hdf5read(magpath, '/Rsm');
Xsm = repmat(Rsm(1,:), Alen,1);
Ysm = repmat(Rsm(2,:), Alen,1);
Zsm = repmat(Rsm(3,:), Alen,1);

GtonT = 1/10^(-9/2+1);
K = hdf5read(magpath, '/K') * GtonT;
Bm = hdf5read(magpath, '/Bm');
Ls = hdf5read(magpath, '/Lstar');
Kp = hdf5read(magpath, '/Kp');

%% Wrap
mag_struct = struct(...
    'Tday', Tday, ...
    'PArad', PArad, ...
    'Xsm', Xsm, ...
    'Ysm', Ysm, ...
    'Zsm', Zsm, ...
    'K', K, ...
    'Ls', Ls, ...
    'Kp', Kp, ...
    'IOPT', round(Kp+1), ...
    'Bm', Bm);

end