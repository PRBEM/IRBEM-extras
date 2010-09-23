function [hEtheta,result_code] = rfl_make_hEtheta(inst_info,Egrid,thetagrid,options)
% [hEtheta,result_code] = rfl_make_hEtheta(inst_info,Egrid,thetagrid,options)
% Compute weights for double numerical integral over E, theta, 
% arbitrary geometries   
[hEtheta,result_code] = inst_info.make_hEtheta(inst_info,Egrid,thetagrid,options);


