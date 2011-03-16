function [halphabeta,result_code] = rfl_make_halphabeta(inst_info,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options)
% function [halphabeta,result_code] = rfl_make_halphabeta(inst_info,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options)
% Compute weights for double numerical integral over alpha, beta,
% arbitrary geometries 
% (response is integrated over energy, except for E_TYPE=INT energy
% channels, for which the energy bandwith is assumed to be 1 to avoid inf)

[halphabeta,result_code] = inst_info.make_halphabeta(inst_info,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options);

