% test odc_util functions

util = odc_util(true);

%% Generate T/Y/D/Q Figures from Schulz and Lanzerotti

y = linspace(0,1,1000);
T = util.T(y);
Y = util.Y(y);
D = util.D(y);
Q = util.Q(y);

% Fig 7, page 21
figure;
plot(y,Y,'k-',y,2*T,'r-',y,2*D,'k--',y,3*D./T,'r--');
xlabel('y=sin\alpha_{eq}');
axis([0 1 0 3]);
legend('Y','2*T','2*D','3*D/T');
title('S&L fig. 7, pg 21');

% Fig 25, page 90
figure;
plot(y,abs(Q./D)/180,'.-');
xlabel('y=sin\alpha_{eq}');
axis([0 1 0 1]);
ylabel('|Q(y)/D(y)/180| = [D_{LL}(y)./D_{LL}(1)]^{1/2}');
title('S&L fig. 25, pg 90');

%% Generate frequencies figure from S&L (Fig 6, pg 13)

L = 2.^(0:0.1:3); % L= 1 - 8
MeV = 10.^(-3:0.1:3);
[L,MeV] = meshgrid(L,MeV);
PitchAngle = 90;
MagLat = 0;

figure;
subplot(3,2,1);
Tg = util.GyroPeriod('p',MeV,MagLat,L);
[cs,h] = contour(L,MeV,1./Tg,unique([10.^(0:2),3*10.^(0:2)]),'k-');
set(gca,'xscale','log','yscale','log','xtick',[1 2 4 8],'ytick',10.^(-2:2:2),'xticklabel',[]);
hlab = clabel(cs);
set(hlab(end-2),'string',[get(hlab(end-2),'string'),' Hz']);
title('PROTONS');

subplot(3,2,2);
Tg = util.GyroPeriod('e',MeV,MagLat,L);
[cs,h] = contour(L,MeV,1e-3./Tg,unique([10.^(-2:2),3*10.^(-2:2)]),'k-');
set(gca,'xscale','log','yscale','log','xtick',[1 2 4 8],'ytick',10.^(-2:2:2),'xticklabel',[],'yticklabel',[]);
set(gca,'YAxisLocation','right');
ylabel('GYRATION');
hlab = clabel(cs);
set(hlab(end),'string',[get(hlab(end),'string'),' kHz']);
title('ELECTRONS');

subplot(3,2,3);
Tb = util.BouncePeriod('p',MeV,PitchAngle,L);
[cs,h] = contour(L,MeV,1./Tb,unique([10.^(-2:1),3*10.^(-2:0)]),'k-');
set(gca,'xscale','log','yscale','log','xtick',[1 2 4 8],'ytick',10.^(-2:2:2),'xticklabel',[]);
hlab = clabel(cs);
set(hlab(end-6),'string',[get(hlab(end-6),'string'),' Hz']);

subplot(3,2,4);
Tb = util.BouncePeriod('e',MeV,PitchAngle,L);
[cs,h] = contour(L,MeV,1./Tb,unique([10.^(0:1),3*10.^(-1:0)]),'k-');
set(gca,'xscale','log','yscale','log','xtick',[1 2 4 8],'ytick',10.^(-2:2:2),'xticklabel',[],'yticklabel',[]);
set(gca,'YAxisLocation','right');
ylabel('BOUNCE');
hlab = clabel(cs);
set(hlab(end-2),'string',[get(hlab(end-2),'string'),' Hz']);

subplot(3,2,5);
Td = util.DriftPeriod('p',MeV,PitchAngle,L);
[cs,h] = contour(L,MeV,1e3./Td,unique([10.^(-2:3),3*10.^(-3:2)]),'k-');
set(gca,'xscale','log','yscale','log','xtick',[1 2 4 8],'ytick',10.^(-2:2:2));
hlab = clabel(cs);
set(hlab(12),'string',[get(hlab(12),'string'),' mHz']);

subplot(3,2,6);
Td = util.DriftPeriod('e',MeV,PitchAngle,L);
[cs,h] = contour(L,MeV,1e3./Td,unique([10.^(-2:3),3*10.^(-3:2)]),'k-');
set(gca,'xscale','log','yscale','log','xtick',[1 2 4 8],'ytick',10.^(-2:2:2),'yticklabel',[]);
set(gca,'YAxisLocation','right');
ylabel('DRIFT');
hlab = clabel(cs);
set(hlab(12),'string',[get(hlab(12),'string'),' mHz']);
