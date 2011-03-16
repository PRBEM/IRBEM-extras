function inst_info = rfl_init_pinhole(inst_info)
% initialize a pinhole angular response
% assumes inst_info is initialized to a generic cylindrically symmetric response
% will overload lots of inst_info.internal functions

% uses G, BIDIRECTIONAL

inst_info.internal.A = @(inst_info,theta,phi)(theta==0)*inst_info.G; % bad, but unaviodable

inst_info.internal.hAthetaphi = @make_hAthetaphi_pinhole;
inst_info.internal.hAtheta = @make_hAtheta_pinhole;
inst_info.internal.hAalphabeta = @make_hAalphabeta_pinhole;
inst_info.internal.hAalpha = @make_hAalpha_pinhole;
inst_info.internal.hA0 = inst_info.G;

function [hAthetaphi,result_code] = make_hAthetaphi_pinhole(inst_info,thetagrid,phigrid,options)
dp = rfl_make_deltas(phigrid,options); % no need to convert to radians, see next line
result_code = 1; % success
h = zeros(length(thetagrid),length(phigrid));
[dtheta,i0] = min(abs(thetagrid)); % find theta nearest 0
h(i0,:) = dp;
if inst_info.internal.bidirectional,
    [dtheta,i0] = min(abs(thetagrid-180)); % find theta nearest 180
    if dtheta < 90,  % just in case user is implicitly assuming not bidirectional
        h(i0,:) = dp;
    end
end
h = inst_info.G/sum(h)*h; % force sum(h) = G
hAthetaphi = h; % success, returns h

function [hAtheta,result_code] = make_hAtheta_pinhole(inst_info,thetagrid,options)
result_code = 1; % success
h = zeros(length(thetagrid),1);
[dtheta,i0] = min(abs(thetagrid)); % find theta nearest 0
h(i0) = 1;
if inst_info.internal.bidirectional,
    [dtheta,i0] = min(abs(thetagrid-180)); % find theta nearest 180
    if dtheta < 90, % just in case user is implicitly assuming not bidirectional
        h(i0) = 1;
    end
end
h = inst_info.G/sum(h)*h; % force sum(h) = G
hAtheta = h; % success, returns h

function [hAalphabeta,result_code] = make_hAalphabeta_pinhole(inst_info,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options)

if inst_info.internal.bidirectional,
    tmpinst_info = inst_info;
    tmpinst_info.internal.bidirectional=false;
    [hAalphabeta1,result_code] = make_hAalphabeta_pinhole(tmpinst_info,alphagrid,betagrid,tgrid,180-alpha0,180+beta0,phib,options);
    if result_code ~= 1,
        return
    end
    [hAalphabeta,result_code] = make_hAalphabeta_pinhole(inst_info,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options);
    if result_code ~= 1,
        return
    end
    hAalphabeta = (hAalphabeta1+hAalphabeta)/2;
    return
end

if isfield(options,'acute'),
    acute = options.acute;
else
    acute = alphagrid(end)<=90;
end

if acute,
    alpha0 = min(alpha0,180-alpha0);
end

dt = rfl_make_deltas(tgrid,options);

NA = length(alphagrid);
NB = length(betagrid);
h = zeros(NA,NB); % full matrix
for it = 1:length(tgrid),
    IA = rfl_get_list_neighbors(alphagrid,alpha0(it));
    IB = rfl_get_list_neighbors(betagrid,beta0(it));
    dh = sparse(NA,NB);
    j = IA(1); % alpha(j) <= alpha0 < alpha(j+1)
    if (j>=1) && (j < NA),
        da = (alphagrid(j+1)-alpha0(it))/(alphagrid(j+1)-alphagrid(j));
        k = IB(1); % beta(k) <= beta0 < beta(k+1)
        if (k>=1) && (k < NB),
            db = (betagrid(k+1)-beta0(it))/(betagrid(k+1)-betaagrid(k));
            dh(j,k) = da*db;
        end
        k = IB(2); % beta(k-1) < beta0 < beta(k)
        if (k>=2) && (k <= NB),
            db = (beta0(it)-betagrid(k-1))/(betagrid(k)-betaagrid(k-1));
            dh(j,k) = da*db;
        end
    end
    j = IA(2); % alpha(j-1) < alpha0 < alpha(j)
    if (j>=2) && (j <= NA),
        da = (alpha0(it)-alphagrid(j))/(alphagrid(j)-alphagrid(j-1));
        k = IB(1); % beta(k) <= beta0 < beta(k+1)
        if (k>=1) && (k < NB),
            db = (betagrid(k+1)-beta0(it))/(betagrid(k+1)-betaagrid(k));
            dh(j,k) = da*db;
        end
        k = IB(2); % beta(k-1) < beta0 < beta(k)
        if (k>=2) && (k <= NB),
            db = (beta0(it)-betagrid(k-1))/(betagrid(k)-betaagrid(k-1));
            dh(j,k) = da*db;
        end
    end    
    h = h+dt(it)*dh;
end
h = h*inst_info.G*sum(dt)/sum(h); % force sum(h) = G*sum(dt)
result_code = 1; % success
hAalphabeta = h; % success, returns h

function [hAalpha,result_code] = make_hAalpha_pinhole(inst_info,alphagrid,tgrid,alpha0,beta0,phib,options)

if inst_info.internal.bidirectional,
    tmpinst_info = inst_info;
    tmpinst_info.internal.bidirectional=false;
    [hAalpha1,result_code] = make_hAalpha_pinhole(tmpinst_info,alphagrid,tgrid,180-alpha0,180+beta0,phib,options);
    if result_code ~= 1,
        return
    end
    [hAalpha,result_code] = make_hAalpha_pinhole(tmpinst_info,alphagrid,tgrid,180-alpha0,180+beta0,phib,options);
    if result_code ~= 1,
        return
    end
    hAalpha = (hAalpha+hAlpha1)/2;
    return
end

if isfield(options,'acute'),
    acute = options.acute;
else
    acute = alphagrid(end)<=90;
end

if acute,
    alpha0 = min(alpha0,180-alpha0);
end

dt = rfl_make_deltas(tgrid,options);

NA = length(alphagrid);
h = zeros(NA,1); % full matrix
for it = 1:length(tgrid),
    IA = rfl_get_list_neighbors(alphagrid,alpha0(it));
    dh = sparse(NA,1);
    j = IA(1); % alpha(j) <= alpha0 < alpha(j+1)
    if (j>=1) && (j < NA),
        dh(j) = (alphagrid(j+1)-alpha0(it))/(alphagrid(j+1)-alphagrid(j));
    end
    j = IA(2); % alpha(j-1) < alpha0 < alpha(j)
    if (j>=2) && (j <= NA),
        dh(j) = (alpha0(it)-alphagrid(j))/(alphagrid(j)-alphagrid(j-1));
    end    
    h = h+dt(it)*dh;
end
h = h*inst_info.G*sum(dt)/sum(h); % force sum(h) = G*sum(dt)
result_code = 1; % success
hAalpha = h; % success, returns h

