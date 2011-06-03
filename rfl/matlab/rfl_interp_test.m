% test RFL interpolation routines

xgrid = linspace(0,1,10)';
ygrid = linspace(0,2,10)';
zgrid = linspace(0,3,10)';

Nrand = 1000;
xhat = rand(Nrand,1);
yhat = rand(Nrand,1)*2;
zhat = rand(Nrand,1)*3;

% test 1-D interp

v0 = xgrid*0.5+4;
v1 = interp1(xgrid,v0,xhat,'linear');
H = rfl_interp_weights_1d(xgrid,xhat);
v2 = H*v0;
figure;
plot(v1-v2,'.');
ylabel('RFL - interp1');
title('1-D');


% test 2-D interp
[XI,YI] = ndgrid(xgrid,ygrid);

v0 = XI*0.5-0.7*YI+4;
v1 = interpn(xgrid,ygrid,v0,xhat,yhat,'linear');
H = rfl_interp_weights_2d(xgrid,ygrid,xhat,yhat);
v2 = H*v0(:);
figure;
plot(v1-v2,'.');
ylabel('RFL - interpn');
title('2-D');

% test 3-D interp
[XI,YI,ZI] = ndgrid(xgrid,ygrid,zgrid);

v0 = XI*0.5-0.7*YI+0.01*ZI+4;
v1 = interpn(xgrid,ygrid,zgrid,v0,xhat,yhat,zhat,'linear');
H = rfl_interp_weights_3d(xgrid,ygrid,zgrid,xhat,yhat,zhat);
v2 = H*v0(:);
figure;
plot(v1-v2,'.');
ylabel('RFL - interpn');
title('3-D');

