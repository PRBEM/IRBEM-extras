function inst_info = rfl_init_inseparable(inst_info)
% inst_info = rfl_init_inseparable(inst_info)
% initialize sensor with inseparable E,theta,phi response
% RESP_TYPE = '[E,TH,PH]';
% generic response function in E, theta, phi coordinates

% min/max forces plateau extrapolation in E
inst_info.R = @(inst_info,E,theta,phi)interpn(inst_info.E_GRID,inst_info.TH_GRID,inst_info.PH_GRID,inst_info.R,min(max(E,inst_info.E_GRID(1)),inst_info.E_GRID(end)),theta,phi,'linear',0); % set to zero outside grid

% make_h for grids in instrument coordinates
inst_info.make_hEthetaphi = @inseparable_hEthetaphi; % on grid in E,theta,phi
inst_info.make_hEtheta = @inseparable_hEtheta; % on grid in E,theta
inst_info.make_hE = @inseparable_hE; % on grid in E
% make_h for grids in magnetic coordinates over time
inst_info.make_hEalphabeta = @inseparable_hEalphabeta; % on grid in E,alpha,beta
inst_info.make_hEalpha = @inseparable_hEalpha; % on grid in E,alpha
inst_info.make_hEiso = @inseparable_hEiso; % on grid in E

% one private variable
inst_info.internal.bidirectional = strcmpi(inst_info.BIDIRECTIONAL,'TRUE'); % convert to logical for faster execution

function [hEthetaphi,result_code] = inseparable_hEthetaphi(inst_info,Egrid,thetagrid,phigrid,options)
dE = rfl_make_deltas(Egrid,options);
dcost = rfl_make_deltas(-cosd(thetagrid),options);
dp = rfl_make_deltas(phigrid,options)*pi/180;
[dE,dcost,dp] = ndgrid(dE,dcost,dp);
[E,theta,phi] = ndgrid(Egrid,thetagrid,phigrid);

result_code = 1; % success
R = inst_info.R(inst_info,E,theta,phi);
h = R.*dE.*dcost.*dp/inst_info.XCAL;
hEthetaphi = h; % success, returns h

function [hEtheta,result_code] = inseparable_hEtheta(inst_info,Egrid,thetagrid,options)
phigrid = rfl_make_grid(0,360,'phi',options);
[hEthetaphi,result_code] = inseparable_hEthetaphi(inst_info,Egrid,thetagrid,phigrid,options);
if result_code == 1,
    hEtheta = sum(hEthetaphi,3);
else
    hEtheta = nan;
end

function [hE,result_code] = inseparable_hE(inst_info,Egrid,options)
thetagrid = rfl_make_grid(0,180,'theta',options);
[hEtheta,result_code] = inseparable_hEtheta(inst_info,Egrid,thetagrid,options);
if result_code == 1,
    hE = sum(hEtheta,2);
else
    hE = nan;
end

function [hEalphabeta,result_code] = inseparable_hEalphabeta(inst_info,Egrid,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options)
dE = rfl_make_deltas(Egrid,options);
if isfield(options,'acute'),
    acute = options.acute;
else
    acute = alphagrid(end)<=90;
end
dcosa = rfl_make_deltas(-cosd(alphagrid),options);
db = rfl_make_deltas(betagrid,options)*pi/180;
dt = rfl_make_deltas(tgrid,options);
[dE,dcosa,db] = ndgrid(dE,dcosa,db);
[E,alpha,beta] = ndgrid(Egrid,alphagrid,betagrid);
% set up smaller alpha/beta grid for rotation to theta/phi
NE = length(Egrid);
[aI,bI] = ndgrid(alphagrid,betagrid);
theta = shiftdim(nan(size(aI)),-1);
phi = theta;
% prepare to loop
result_code = 1; % success
hEalphabeta = nan; % on error, returns a nan
h = 0;
for it = 1:length(tgrid),
    [theta(:),phi(:),result_code] = rfl_alphabeta2thetaphi(aI(:),bI(:),alpha0(it),beta0(it),phib(it));
    if result_code ~= 1,
        return
    end
    h = h+inst_info.R(inst_info,E,repmat(theta,[NE,1,1]),repmat(phi,[NE,1,1]))*dt(it);
    if acute, % add pi-alpha point as well
        [theta(:),phi(:),result_code] = rfl_alphabeta2thetaphi(180-aI(:),bI(:),alpha0(it),beta0(it),phib(it));
        if result_code ~= 1,
            return
        end
        h = h+inst_info.R(inst_info,E,repmat(theta,[NE,1,1]),repmat(phi,[NE,1,1]))*dt(it);
    end
end
h = h.*dE.*dcosa.*db/inst_info.XCAL;
hEalphabeta = h; % success, returns h

function [hEalpha,result_code] = inseparable_hEalpha(inst_info,Egrid,alphagrid,tgrid,alpha0,beta0,phib,options)
betagrid = rfl_make_grid(0,360,'beta',options);
[hEalphabeta,result_code] = inseparable_hEalphabeta(inst_info,Egrid,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options);
if result_code == 1,
    hEalpha = sum(hEalphabeta,3);
else
    hEalpha = nan;
end

function [hEiso,result_code] = inseparable_hEiso(inst_info,Egrid,tgrid,options)
[hE,result_code] = inseparable_hE(inst_info,Egrid,options);
if result_code == 1,
    dt = rfl_make_deltas(tgrid,options);
    hEiso = hE*sum(dt);
else
    hEiso = nan;
end
