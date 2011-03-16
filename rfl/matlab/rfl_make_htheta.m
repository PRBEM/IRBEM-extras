function [htheta,result_code] = rfl_make_htheta(inst_info,thetagrid,options)
% [htheta,result_code] = rfl_make_hthetaphi(inst_info,thetagrid,options)
% Compute weights for numerical integral over theta
% arbitrary geometries   
% (response is integrated over phi)
% (response is integrated over energy, except for E_TYPE=INT energy
% channels, for which the energy bandwith is assumed to be 1 to avoid inf)

[htheta,result_code] = inst_info.make_htheta(inst_info,thetagrid,options);

