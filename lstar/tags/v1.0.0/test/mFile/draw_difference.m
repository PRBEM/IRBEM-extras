
clc
clear
clf
colormap([[1 1 1]*0; graphics.rainbow(63)])
set(gcf,'papersize',[6 5],'paperposition',[.1 .1 5.8 4.8])
set(0,'defaultaxesfontname','times')
set(0,'defaulttextfontname','times')
set(0,'defaulttextfontsize',11)
set(0,'defaultaxesfontsize',12)

resultroot = '../tmp/result';

%% Loop over files
list = dir( fullfile(resultroot, '*.mat') );
eT = NaN(length(list), 2);

idx = 0;
for entry = list'
    idx = idx + 1;
    name = entry.name;
    Lspath = fullfile(resultroot, name);
    disp(['%% Processing ' name])
    
    %% Load
    load( Lspath );
    eT(idx, 1) = lstar.Tday(1);
    eT(idx, 2) = lstar.elapsedT;
    
    %% Analyze
    dLs = abs(lstar.Ls - mag.Ls)./lstar.Ls;
    
    %% Draw
    clf
    gaps = [.02 .03];
    position = [.09 .09 .8 .84];
    clim1 = [0.8 6];
    clim2 = [-4.01 -1];
    
    hca = [];
    hca(1) = graphics.subpanel(2,1,1, gaps,position);
    imagesc(mag.Tday, mag.PArad(:,1)*180/pi, lstar.Ls)
    datetick('keepticks')
    caxis(clim1)
    ylabel('\alpha [\circ]')
    
    hca(2) = graphics.subpanel(2,1,2, gaps,position);
    imagesc(mag.Tday, mag.PArad(:,1)*180/pi, log10(dLs))
    set(gca, 'xtick', linspace(mag.Tday(1),mag.Tday(end),9))
    datetick('keepticks')
    caxis(clim2)
    ylabel('\alpha [\circ]')
    
    title(hca(1),...
        [datestr(mag.Tday(1), 'yyyy-mm-dd') ...
        ', Elapsed T=' sprintf('%.1f sec',lstar.elapsedT)])
    yticks = [5 15:15:90];
    set(hca,'ydir','normal', 'ytick',yticks,...
        'xtick', linspace(mag.Tday(1),mag.Tday(end),9),...
        'tickdir','out')
    set(hca(1:end-1),'xticklabel','')
    text(round(mag.Tday(1)),-.16*diff(ylim)+min(ylim),...
        datestr(mag.Tday(1), 'yyyy-mm-dd'),'horizontalalignment','center')
    
    % colorbars
    cbarpos = [sum(position([1 3]))+.02 position(2) gaps(1) position(4)];
    graphics.subpanel(8,1,1:3, gaps,cbarpos)
    caxis(clim1)
    axis off
    xlabel(colorbar('position',get(gca,'position'),...
        'ytick',[1 2 4 6],'fontsize',10,'tickdir','out'),...
        '{\itL}*',...
        'horizontalalignment','left','fontsize',10)
    
    graphics.subpanel(8,1,5:7, gaps,cbarpos)
    caxis(10.^clim2)
    axis off
    xlabel(colorbar('position',get(gca,'position'),...
        'yscale','log','ytick',logspace(-4,-1,4),...
        'fontsize',10,'tickdir','out'),...
        {'|\Delta{\itL}*|' '----' '  {\itL}*'},...
        'horizontalalignment','left','fontsize',10)
    return
    %% Save
    epsname = name(1:end-4);
    print('-depsc2',['../tmp/eps/' epsname])
    
    disp(['%% Processed ' name])
    
    %break
end

%% Plot elapsed T
clf
set(gcf,'papersize',[6 4],'paperposition',[.1 .1 5.8 3.8])
position = [.1 .09 .87 .86];
graphics.subpanel(1,1,1, gaps,position)
plot(eT(:,1),eT(:,2),'o-b','markersize',4)
xlim(eT([1 end],1)+[-10; 10])
datetick('keeplimits')
ylabel('Elapsed T [s]')

print('-depsc2','../tmp/eps/elapsedT')