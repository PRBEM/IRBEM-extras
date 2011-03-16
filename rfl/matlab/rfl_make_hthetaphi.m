function [hthetaphi,result_code] = rfl_make_hthetaphi(inst_info,thetagrid,phigrid,options)
% [hthetaphi,result_code] = rfl_make_hthetaphi(inst_info,thetagrid,phigrid,options)
% Compute weights for double numerical integral over theta, phi,
% arbitrary geometries   
% (response is integrated over energy, except for E_TYPE=INT energy
% channels, for which the energy bandwith is assumed to be 1 to avoid inf)

[hthetaphi,result_code] = inst_info.make_hthetaphi(inst_info,thetagrid,phigrid,options);

