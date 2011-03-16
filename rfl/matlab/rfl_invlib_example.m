% example use of rfl and invlib
% energy spectral inversion for ICO
% wide2uni inversion for SAMPEX/PET
% angular inversion for simulated spinning PET using PC's

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
H = resp.make_halpha(resp,PAgrid,tgrid,alpha0,beta0,phib,options);
H = H'/sum(H); % angular response should sum to 1 if wideflux is in same units as uniflux

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
plot(PAgrid,flux,'k-',PAgrid,H/max(H),'b-',alpha0,wideflux,'bx',alpha0*[1 1],wideflux*exp([-1 1]*dlogwideflux),'b.-',...
    alpha,uniflux,'ro',alpha*[1 1],uniflux*exp([-1 1]*dloguniflux),'r.-');
xlabel('Pitch Angle, degrees');
set(gca,'xtick',0:30:180);
legend('Truth: sin^2\alpha','Normalized angular response','wideflux@\alpha_0','\pm 1 stdev','uniflux (fit)','\pm 1 stdev','location','so');
grid on;
title(sprintf('Angular inversion using %s and PET angular response',method));

%% local pitch angle inversion using PC's in equatorial pitch angle
% fake spinner
BB0 = 2; % B/Bequator
PAgrid = (0:2:180)'; % local pitch angle, LPA
%% make up orientation angles
dt = 1;
times = (0:dt:30)'; % time tags
Nchans = length(times);
alpha0 = times*360/times(end); % one spin of fake alpha0 values
beta0 = 125;
phib = 211;
resp = pet.(pet.CHANNEL_NAMES{5}).PROT; % pick a random channel
local_flux = 100*sind(PAgrid).^2.5; % slightly more peaked than sin^2

%% generate H for one spin
options = struct('Nint',1e3); % better angular integrals in beta
H = nan(length(times),length(PAgrid));
for i = 1:length(times),
    H(i,:) = resp.make_halpha(resp,PAgrid,dt,alpha0(i),beta0,phib,options);
end

dy = repmat(resp.XCAL_RMSE,[1,Nchans]);
b = zeros(size(dy));

%% generate fake counts
lambda = H*local_flux;
try % with stats toolbox
    y = round(poissrnd(lambda'.*exp(randn(1,Nchans).*dy))); % counts with Poisson and lognormal errors
catch % without, no poissrnd
    y = round(lambda'.*exp(randn(1,Nchans).*dy)); % counts with lognormal errors
end


%% define fake principal components in equatorial pitch angle
pc.EPA = (0:90)'; % equatorial pitch angle
pc.mean_log_flux = sind(pc.EPA).^2;
pc.Nq = 10;
pc.basis = nan(length(pc.EPA),pc.Nq);
for i = 1:(pc.Nq/2),
    j = i*2-1;
    pc.basis(:,j) = sind(2*pc.EPA*i);
    pc.basis(:,j+1) = cosd(2*pc.EPA*i);
end
pc.variance = (1:pc.Nq)'.^-2;

%% get interpolation wieghts from EPA to local PA
local_epa = abs(asind(sind(PAgrid)/sqrt(BB0)));
Hinterp = rfl_interp_weights_1d(pc.EPA,local_epa); % projects flux vs EPA to flux vs LPA
Hepa = H*Hinterp;
    
%% cajole pc_spec_inv to do angular inversion
fit = invlib('pc_spec_inv',y,dy,ones(size(PAgrid)),Hepa,1,b,pc.mean_log_flux,pc.basis,pc.variance,'G=GdE');
% note, to get pc_spec_inv to do generic PC inversion, 
% replace Egrid with ones, pass dt = 1, and set G=GdE, and 
fit.local_flux = Hinterp*fit.flux';
fit.local_flux_1stdev = Hinterp*[fit.flux'.*exp(-fit.dlogflux') fit.flux'.*exp(fit.dlogflux')];

%% plot
figure;
subplot(2,1,1);
plot(PAgrid,H);
ylabel('Sector Response Weights');
subplot(2,1,2);
plot(PAgrid,local_flux,'k-',PAgrid,fit.local_flux,'r-',PAgrid*[1 1],fit.local_flux_1stdev,'r--');
xlabel('Local Pitch Angle, deg');
ylabel('Flux');
legend('True Flux','PC inversion','\pm 1 stdev');


