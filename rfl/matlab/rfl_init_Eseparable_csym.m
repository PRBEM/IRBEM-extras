function inst_info = rfl_init_Eseparable_csym(inst_info)
% inst_info = rfl_init_Eseparable_csym(inst_info)
% populate methods and extra properties for
% an [E],[TH] response type
% inst_info is assumed to be a fully-initialized Eseparable sensor
% this function just overloads certain methods

inst_info.internal.hAthetaphi = @make_hAthetaphi_csym;
inst_info.internal.hAtheta = @make_hAtheta_csym;

switch(inst_info.RESP_TYPE),
    case '[E],[TH]',
        switch(inst_info.TH_TYPE),
            case 'TBL',
                inst_info.internal.A = @(inst_info,theta,phi)interpn(inst_info.TH_GRID,inst_info.A,theta,'linear',0); % zero outside grid, ignore phi
            case 'PINHOLE',
                inst_info = rfl_init_pinhole(inst_info);
            case 'CYL_TELE', %%
                inst_info = init_telescope_csym(inst_info);
            case 'SLAB',
                inst_info = init_single_element(inst_info,inst_info.H1*inst_info.W1);
            case 'DISK',
                inst_info = init_single_element(inst_info,pi*inst_info.R1^2);
            otherwise
                error('TH_TYPE %s not defined yet',inst_info.TH_TYPE);
        end
    case '[E]', 
        inst_info = rfl_init_omni(inst_info);
    otherwise
        error('RESP_TYPE %s not defined yet',inst_info.RESP_TYPE);
end

%%% accidentally wrote this code that duplicates rfl_init_pinhole.m. Might
%%% actually be better code ;)
% function inst_info = init_pinhole(inst_info)
% % initialize pinhole angular response
% % G, BIDIRECTIONAL
% 
% if inst_info.internal.bidirectional,
%     inst_info.internal.A = @(inst_info,theta,phi)(theta==0)*inst_info.G; % delta function
% else
%     inst_info.internal.A = @(inst_info,theta,phi)((theta==0)|(theta==180))*inst_info.G/2; % delta function
% end
% inst_info.internal.hAalphabeta = @make_hAlphabeta_pinhole;
% inst_info.internal.hAthetaphi = @make_hAthetaphi_pinhole;
% inst_info.internal.hA0 = inst_info.G; % already includes any bidirectional effects
% 
% function [hAthetaphi,result_code] = make_hAthetaphi_pinhole(inst_info,thetagrid,phigrid,options)
% result_code = 1; % success
% h = zeros(length(thetagrid),length(phigrid));
% h(thetagrid==0) = inst_info.G; % delta function at theta=0
% if inst_info.internal.bidirectional,
%     h(thetagrid==180) = inst_info.G; % add backward response
%     h = h/2; % G includes bidirectional effect
% end
% hAthetaphi = h; % success, returns h
% 
% 
% function [hAalphabeta,result_code] = make_hAalphabeta_pinhole(inst_info,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options)
% if isfield(options,'acute'),
%     acute = options.acute;
% else
%     acute = alphagrid(end)<=90;
% end
% dt = rfl_make_deltas(tgrid,options);
% % prepare to loop
% result_code = 1; % success
% hAalphabeta = nan; % on error, returns a nan
% h = 0;
% for it = 1:length(tgrid),
%     alpha_pinhole = alpha0(it);
%     beta_pinhole = beta0(it);
%     if acute,
%         alpha_pinhole = min(alpha_pinhole,180-alpha_pinhole);
%     end
%     if result_code ~= 1,
%         return
%     end
%     va = rfl_lbf_eval(alphagrid,alpha_pinhole);
%     vb = rfl_lbf_eval(betagrid,beta_pinhole,{'period',360});
%     vab = va*vb';
% 
%     h = h+inst_info.G*dt(it)*vab;
%     
%     if inst_info.bidirectional,
%         alpha_pinhole = 180-alpha0(it); % opposite hemisphere on [0,180]
%         beta_pinhole = 180+beta0(it); % opposite hemisphere on [0,360]
%         if acute,
%             alpha_pinhole = min(alpha_pinhole,180-alpha_pinhole);
%         end
%         if result_code ~= 1,
%             return
%         end
%         va = rfl_lbf_eval(alphagrid,alpha_pinhole);
%         vb = rfl_lbf_eval(betagrid,beta_pinhole,{'period',360});
%         vab = va*vb';
%         
%         h = h+inst_info.G*dt(it)*vab;
%     end
% end
% if inst_info.bidirectional,
%     h = h/2; % G already includes bidirectional effects
% end
% hAalphabeta = h; % success, returns h
% 
% 
function inst_info = init_telescope_csym(inst_info)
% initialize cylindrically symmetric 2-element telescope
% overloads inst_info.internal.A, sets some inst_info.internal fields, too
% R1, R2, D, BIDIRECTIONAL

% Using Sullivan's 1971 paper as updated in ~2010
inst_info.internal.Rs = min(inst_info.R1,inst_info.R2);
inst_info.internal.thetac = atand(abs(inst_info.R1-inst_info.R2)/inst_info.D);
inst_info.internal.thetam = atand((inst_info.R1+inst_info.R2)/inst_info.D);
inst_info.internal.A = @telescope_csym_A;
tmp = inst_info.R1^2+inst_info.R2^2+inst_info.D^2;
inst_info.internal.hA0 = pi^2/2*(tmp - sqrt(tmp^2-4*inst_info.R1^2*inst_info.R2^2));  % (eqn 8)
if inst_info.internal.bidirectional,
    inst_info.internal.hA0 = inst_info.internal.hA0*2; % double-headed telescope has twice the area/geometric factor
end


function A = telescope_csym_A(inst_info,theta,phi)
% effective area for 2-element cylindrical telescope
% ignores phi
% right out of Sullivan's updated paper, equation (10)
A = zeros(size(theta));

if inst_info.internal.bidirectional,
    theta = min(theta,180-theta);
end

f = (theta <= inst_info.internal.thetac);
if any(f(:)),
    A(f) = pi*inst_info.internal.Rs^2*cosd(theta(f));
end

f = (theta > inst_info.internal.thetac) & (theta < inst_info.internal.thetam);
if any(f(:)),
    tantheta = tand(theta(f));
    Psi1 = acos((inst_info.R1^2 + inst_info.D^2*tantheta.^2 - inst_info.R2^2) ./ (2*inst_info.D*inst_info.R1*tantheta)); % rads
    Psi2 = acos((inst_info.R2^2 + inst_info.D^2*tantheta.^2 - inst_info.R1^2) ./ (2*inst_info.D*inst_info.R2*tantheta)); % rads
    A(f) = cosd(theta(f))./2.*( inst_info.R1^2*(2*Psi1 - sin(2*Psi1)) + inst_info.R2^2*(2*Psi2-sin(2*Psi2)) );
end

function inst_info = init_single_element(inst_info,area)
% initialize angle methods for single element telescopes
% (DISK or SLAB)
% for a single element with specified area

% Overload .A method (effective area)

inst_info.internal.area = area;

if inst_info.internal.bidirectional,
    inst_info.internal.A = @(inst_info,theta,phi)inst_info.internal.area*cosd(theta);
    inst_info.internal.hA0 = area*2*pi;
else % this is the case Sullivan treats
    inst_info.internal.A = @(inst_info,theta,phi)inst_info.internal.area*cosd(theta).*(theta<=90);
    inst_info.internal.hA0 = area*pi;
end


function inst_info = rfl_init_omni(inst_info)
% populate methods and variables for
% make_hA* part of response function with
% idealized omnidirectional channel

if inst_info.internal.bidirectional,
    inst_info.internal.A = @(inst_info,theta,phi)inst_info.G/4/pi;
    inst_info.internal.hAthetaphi = @fullomni_hAthetaphi;
    inst_info.internal.hAtheta = @fullomni_hAtheta;
    inst_info.internal.hAalphabeta = @(inst_info,alphagrid,betagrid,alpha0,beta0,phib,options)inst_info.internal.hAthetaphi(inst_info,alphagrid,betagrid,options);
    inst_info.internal.hAalpha = @(inst_info,alphagrid,alpha0,beta0,phib,options)inst_info.internal.hAtheta(inst_info,alphagrid,options);
else % this is the case Sullivan treats
    inst_info.internal.A = @(inst_info,theta,phi)inst_info.G*(theta<=90)/2/pi;
    inst_info.internal.hAthetaphi = @halfomni_hAthetaphi;
    inst_info.internal.hAtheta = @halfomni_hAtheta;
    % keep standard csym separable hAalphabeta and hAalpha
end
inst_info.internal.hA0 = inst_info.G;  % already includes any bidirectional effects

function [hAthetaphi,result_code] = fullomni_hAthetaphi(inst_info,thetagrid,phigrid,options)
result_code = 1;
hAthetaphi = inst_info.G/4/pi*rfl_make_deltas(-cosd(thetagrid(:)),options)*rfl_make_deltas(phigrid(:),options)';

function [hAtheta,result_code] = fullomni_hAtheta(inst_info,thetagrid,options)
result_code = 1;
hAtheta = inst_info.G/2*rfl_make_deltas(-cosd(thetagrid(:)),options);

function [hAthetaphi,result_code] = halfomni_hAthetaphi(inst_info,thetagrid,phigrid,options)
result_code = 1;
hAthetaphi = inst_info.G/2/pi*((thetagrid(:)<=90).*rfl_make_deltas(-cosd(thetagrid(:)),options))*rfl_make_deltas(phigrid(:),options)';

function [hAtheta,result_code] = halfomni_hAtheta(inst_info,thetagrid,options)
result_code = 1;
hAtheta = inst_info.G*(thetagrid(:)<=90).*rfl_make_deltas(-cosd(thetagrid(:)),options);

function [hAthetaphi,result_code] = make_hAthetaphi_csym(inst_info,thetagrid,phigrid,options)
dcost = rfl_make_deltas(-cosd(thetagrid),options);
dp = rfl_make_deltas(phigrid,options)*pi/180;

result_code = 1; % success
A = inst_info.internal.A(inst_info,thetagrid,0);
h = (A(:).*dcost(:))*dp(:)'; % outer product, column x row - NTHETA x NPHI matrix
hAthetaphi = h; % success, returns h

function [hAtheta,result_code] = make_hAtheta_csym(inst_info,thetagrid,options)
dcost = rfl_make_deltas(-cosd(thetagrid),options);
result_code = 1; % success
A = inst_info.internal.A(inst_info,thetagrid,0);
h = (A(:).*dcost(:))*2*pi;
hAtheta = h; % success, returns h
