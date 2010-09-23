function [hEthetaphi,result_code] = rfl_make_hEthetaphi(inst_info,Egrid,thetagrid,phigrid,options)
% [hEthetaphi,result_code] = rfl_make_hEthetaphi(inst_info,Egrid,thetagrid,phigrid,options)
% Compute weights for triple numerical integral over E, theta, phi,
% arbitrary geometries   
[hEthetaphi,result_code] = inst_info.make_hEthetaphi(inst_info,Egrid,thetagrid,phigrid,options);

