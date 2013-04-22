
clc
clear
clf
colormap([[1 1 1]*0; jet(63)])
set(gcf,'papersize',[6 6],'paperposition',[.1 .1 5.8 5.8])
set(0,'defaultaxesfontname','times')
set(0,'defaulttextfontname','times')
set(0,'defaulttextfontsize',12)
set(0,'defaultaxesfontsize',12)

magroot = '../tmp/MagEphem';

%% Loop over files
list = dir( fullfile(magroot, '*.h5') );
for entry = list'
    name = entry.name;
    magpath = fullfile(magroot, name);
    
    %% Load h5
    [ mag_struct ] = load_magephem( magpath );
    
    %% Calculate field line parameters
    [ K, Bm, ~, ~, ~, Bmeq, ~, ~, ~, Beq ] = ...
        ubk.fieldline(mag_struct.Xsm(1,:),mag_struct.Ysm(1,:),mag_struct.Zsm(1,:),...
        mag_struct.Tday,mag_struct.IOPT,'t89','igrf',[],[],4);
    
    break
end