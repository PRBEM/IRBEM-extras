% example use of rfl and invlib
% energy spectral inversion for ICO
% angular inversion for SAMPEX/PET

%%%%% Spectral inversion %%%%%%%%

%% load ICO response
cdfname = which('ico.cdf');
ico = rfl_load_inst_info(cdfname);

%% set up problem
Egrid = 10.^(-1:0.1:1)'; % MeV
tgrid = 1; % use 1 to leave time integral out of H
dt = 120; % typical ICO integration time, seconds
Nchans = length(ico.CHANNEL_NAMES);
NE = length(Egrid);
H = sparse(Nchans,NE);
b = zeros(1,Nchans);
dy = zeros(1,Nchans);
options = [];
for i = 1:Nchans,
    resp = ico.(ico.CHANNEL_NAMES{i}).ELE;
    hE = resp.make_hE(resp,Egrid,options);
    H(i,:) = hE;
    dy(i) = resp.XCAL_RMSE;
end

%% make simulated data
flux = 100*(Egrid/1).^-3; % power law spectrum
lambda = H*flux*dt;
try % with stats toolbox
    y = round(poissrnd(lambda'.*exp(randn(1,Nchans).*dy))); % counts with Poisson and lognormal errors
catch % without, no poissrnd
    y = round(lambda'.*exp(randn(1,Nchans).*dy)); % counts with lognormal errors
end

%% do the energy inversion
fit = invlib('ana_spec_inv',y,dy,Egrid,H,dt,b,Egrid,'G=GdE','PL');

%% plot
figure;
pos = get(gcf,'pos');
set(gcf,'pos',pos+[0 -pos(4)/2 0 pos(4)/2]);
subplot(2,1,1);
loglog(Egrid,flux,'k-',Egrid,fit.flux,'r.-',Egrid,fit.flux.*exp(fit.dlogflux),'r-',Egrid,fit.flux.*exp(-fit.dlogflux),'r-');
xlabel('MeV');
ylabel('Flux, #/cm^2/s/sr/MeV');
title('ICO spectral inversion');
legend('Truth','Fit','\pm 1 stdev');
grid on;

subplot(2,1,2);
loglog(Egrid,H);
xlabel('MeV');
ylabel('Response, cm^2 sr MeV');
grid on

%%%%% Angular inversion %%%%%%%%

%% load PET response - pretend PROT is for electrons
cdfname = which('sampex_pet.cdf');
pet = rfl_load_inst_info(cdfname);

%% set up problem
Egrid = 10.^(1:0.1:3)'; % energy grid, MeV
PAgrid = (0:2:180)';
tgrid = 1; % use 1 to leave time integral out of H
options = [];
resp = pet.(pet.CHANNEL_NAMES{5}).PROT; % pick a random channel
% make up orientation angles
alpha0 = 75;
beta0 = 125;
phib = 211;
warning('Obtaining halpha from make_hEalpha because halpha doesn''t exist yet');
H = resp.make_hEalpha(resp,Egrid,PAgrid,tgrid,alpha0,beta0,phib,options);
H = sum(H,1); % integral over energy
H = H/sum(H); % angular response should sum to 1 if wideflux is in same units as uniflux

%% make simulated data
flux = sind(PAgrid).^2;
dlogwideflux = resp.XCAL_RMSE; % usually would get error from some other source
wideflux = H*flux*exp(randn(1)*dlogwideflux);
keV = 500; % fake energy value
method = 'Vampola';
L = 5;
BB0 = 2;
alpha = 90; % target pitch angle for inverted uniflux

%% do the angular inversion
[uniflux,dloguniflux,result] = invlib('wide2uni',wideflux,dlogwideflux,PAgrid,H,alpha,method,'keV',keV,'Lm',L,'B/B0',BB0);

%% plot
figure;
pos = get(gcf,'pos');
set(gcf,'pos',pos+[0 -pos(4)/2 0 pos(4)/2]);
plot(PAgrid,flux,'k-',PAgrid,H/max(H),'b-',alpha0,wideflux,'bx',alpha,uniflux,'ro',alpha,uniflux*exp([-1 1]*dloguniflux),'r.-');
xlabel('Pitch Angle, degrees');
set(gca,'xtick',0:30:180);
legend('Truth: sin^2\alpha','Normalized angular response','wideflux@\alpha_0','uniflux (fit)','\pm 1 stdev','location','so');
grid on;
title(sprintf('Angular inversion using %s and PET angular response',method));

