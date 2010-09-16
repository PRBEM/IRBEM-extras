function resp = rfl_init_omni(resp)
% populate methods and variables for
% make_h<angles> part of response function with
% idealized omnidirectional channel

% halphabeta = resp.make_alphabeta(resp,alphagrid,betagrid,alpha0,beta0,phib,options);
% halpha = resp.make_alpha(resp,alphagrid,alpha0);
% hthetaphi = resp.make_thetaphi(resp,thetagrid,phigrid,options);
% htheta = resp.make_thetaphi(resp,thetagrid,options);
resp.make_htheta = @rfl_make_htheta_omni;
% homni = resp.make_homni(resp,options);
resp.make_homni = @rfl_make_homni_omni;
error('halpha, halphabeta, and hthetaphi methods not populated');

function [homni,result_code] = rfl_make_homni_omni(resp,options)
% [homni,result_code] = rfl_make_homni(inst_info,options) Compute geometric
% factor omnidirectional instrument (replaces integral over theta/phi or
% alpha/beta)

if isfield(resp,'G'),
    result_code = 1;
    homni = resp.G;
else
    if nargout >= 2,
        result_code = -2;
        homni = [];
    else
        error('Requested omnidirectional approximation for channel with no "G" property');
    end
end

function [htheta,result_code] = rfl_make_htheta_omni(resp,thetagrid,options)
% [htheta,result_code] = rfl_make_htheta_omni(resp,thetagrid,options)

if isfield(resp,'G'),
    result_code = 1;
    dT = rfl_make_deltas(thetagrid,options);
    htheta = dT*sind(thetagrid)*G/2; % G = 2*pi*A*int sin(theta)dtheta = 4*pi*A. so, htheta = 2*pi*A*sin(theta)dtheta = G/2*dtheta*sin(theta)
else
    if nargout >= 2,
        result_code = -2;
        htheta= [];
    else
        error('Requested omnidirectional approximation for channel with no "G" property');
    end
end

function [hthetaphi,result_code] = rfl_make_hthetaphi_omni(resp,thetagrid,phigrid,options)
% [hthetaphi,result_code] = rfl_make_hthetaphi_omni(resp,thetagrid,phigrid,options)

if isfield(resp,'G'),
    result_code = 1;
    dT = rfl_make_deltas(thetagrid,options);
    htheta = dT*sind(thetagrid)*G/2; % G = 2*pi*A*int sin(theta)dtheta = 4*pi*A. so, htheta = 2*pi*A*sin(theta)dtheta = G/2*dtheta*sin(theta)
else
    if nargout >= 2,
        result_code = -2;
        htheta= [];
    else
        error('Requested omnidirectional approximation for channel with no "G" property');
    end
end
