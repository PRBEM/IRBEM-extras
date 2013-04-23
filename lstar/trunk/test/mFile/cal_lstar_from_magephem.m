
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
resultroot = '../tmp/result';

%% Loop over files
list = dir( fullfile(magroot, '*.h5') );
for entry = list'
    name = entry.name;
    magpath = fullfile(magroot, name);
    disp(['%% Processing ' name])
    
    %% Load h5
    [ mag ] = load_magephem( magpath );
    palen = size(mag.PArad,1);
    
    %% Calculate field line parameters
    [ K0, Bm0, ~, ~, ~, Bmeq, Xeq, Yeq, Zeq, Beq ] = ...
        ubk.fieldline(mag.Xsm(1,:),mag.Ysm(1,:),mag.Zsm(1,:),...
        mag.Tday,mag.IOPT,'t89','igrf',[],[],1,8);
    [ Bxsc Bysc Bzsc ] = ...
        ubk.ts_field(mag.Xsm(1,:),mag.Ysm(1,:),mag.Zsm(1,:),...
        mag.Tday,mag.IOPT,'t89','igrf','sm');
    Bsc = sqrt(Bxsc.^2 + Bysc.^2 + Bzsc.^2);
    
    %% Calculate pitch angle at the magnetic equator
    sin2ameq = repmat(Bmeq./Bsc,palen,1) .* ...
        sin(mag.PArad).^2;
    sin2ameq(sin2ameq>1) = NaN;
    PAmeq = asin(sqrt(sin2ameq));
    
    %% Calculate L*
    [p0 r0] = cart2pol(Xeq, Yeq);
    p0 = repmat(p0, palen,1);
    r0 = repmat(r0, palen,1);
    disp('%% Calculating L*')
    tic
    [ Ls K ] = ...
        ubk.lstar(r0,p0,PAmeq,mag.Tday,mag.IOPT,...
        't89','igrf',[], [], [], [], [], [], ...
        [], 4, 4);
    elapsedT = toc;
    disp(['%% elapsedT is ' num2str(elapsedT) ' seconds.'])
    lstar = struct(...
        'p0',p0,...
        'r0',r0,...
        'pa0',PAmeq,...
        'Ls',Ls,...
        'K',K,...
        'Tday',mag.Tday,...
        'IOPT',mag.IOPT,...
        'elapsedT',elapsedT,...
        'magname',name);
    
    clf
    dLs = abs(Ls - mag.Ls)./Ls;
    imagesc(mag.Tday, mag.PArad(:,1)*180/pi, log10(dLs))
    set(gca,'ydir','normal')
    datetick
    caxis([-3 -1])
    colorbar
    title(strrep(name,'_',' '))
    drawnow
    
    %% Save
    newname = [name(1:5) '_' datestr(mag.Tday(1),'yyyymmdd') '.mat'];
    newpath = fullfile(resultroot, newname);
    save(newpath,'mag','lstar')
    
    disp(['%% Processed ' name])
    
    %break
end