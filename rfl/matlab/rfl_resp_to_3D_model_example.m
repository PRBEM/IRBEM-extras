% example PET response to a 3D model in 
% Energy, equatorial pitch angle, and L

%% load PET response
outpath = fileparts(which('rfl_resp_to_3D_model_example.m')); % directory containing this file
cdfname = [outpath,filesep,'sampex_pet.cdf'];
pet = rfl_load_inst_info(cdfname);

%% define local response grid
local.L = 1.7; % random choice, Lshell, interesting place in inner zone
local.BB0 = 2; % random choice, B/B0, somewhere off equator
local.alpha0 = 35; % random choice, pitch angle of boresight
local.beta0 = 235; % random choice, gyrophase of boresight
local.phib = 75; % random choice, longitude of B in instrument coordinates
local.E = 10.^(0:0.05:3)'; % MeV grid
local.alpha = (2:2:178)'; % local pitch angle grid, degrees
local.tgrid = 6; % single entry means dt, seconds

%% define model grid
model.E = 10.^(0:0.1:3.5); % MeV
model.alpha_eq = (5:5:90)'; % equatorial pitch angle, degrees
model.L = (1:0.1:3)'; % inner zone

%% get pet response versus E, alpha, beta
Nchans = numel(pet.CHANNEL_NAMES);
H_pet = nan(Nchans,length(local.E),length(local.alpha));
options = [];

for i = 1:Nchans,
    resp = pet.(pet.CHANNEL_NAMES{i}).PROT;
    [hEalpha,result_code] = resp.make_hEalpha(resp,local.E,local.alpha,local.tgrid,local.alpha0,local.beta0,local.phib,options);
    H_pet(i,:,:) = hEalpha;
end

%% plot instrument response
figure;
subplot(2,1,1);
loglog(local.E,squeeze(sum(H_pet,3)),'.-');
xlabel(pet.E_UNIT);
ylabel('sum H over alpha');

subplot(2,1,2);
semilogy(local.alpha,squeeze(sum(H_pet,2)),'.-');
xlabel('Local Pitch Angle, deg');
ylabel('sum H over E');

%% compute weights to interpolate to model grid

local.alpha_eq = asind(sind(local.alpha)/sqrt(local.BB0)); % compute model coordinates from local coordinates
[Ehat,alpha_eqhat,Lhat] = ndgrid(local.E,local.alpha_eq,local.L);
H_interp = rfl_interp_weights_3d(model.E,model.alpha_eq,model.L,Ehat(:),alpha_eqhat(:),Lhat);

%% plot interpolation weights -- can't do this, it overflows memory
% figure;
% subplot(2,1,1);
% loglog(Ehat(:),H_interp,'.');
% set(gca,'ylim',[0.01 1]);
% xlabel(pet.E_UNIT);
% subplot(2,1,2);
% semilogy(alpha_eqhat(:),H_interp,'.');
% xlabel('Equatorial Pitch Angle, deg');
% set(gca,'ylim',[0.01 1]);

%% combine weights from response function and model interpolation

H = sparse(H_pet(:,:))*H_interp; % make H_pet sparse to save memory
% y = H*flux(:)
% y(i) is expected counts in i'th channel

%% plot interpolation weights
Hmax = max(H(:));
Hax = max(H(:))*[1e-3 1];
[EI,AI,LI] = ndgrid(model.E,model.alpha_eq,model.L);
figure;
subplot(3,1,1);
loglog(EI(:),H,'.');
set(gca,'ylim',Hax);
xlabel(pet.E_UNIT);
ylabel('H');
subplot(3,1,2);
semilogy(AI(:),H,'.');
xlabel('Equatorial Pitch Angle, deg');
set(gca,'ylim',Hax);
ylabel('H');
subplot(3,1,3);
semilogy(LI(:),H,'.');
xlabel('L');
set(gca,'ylim',Hax);
ylabel('H');

