function [halpha,result_code] = rfl_make_halpha(inst_info,alphagrid,tgrid,alpha0,beta0,phib,options)
% function [halpha,result_code] = rfl_make_halpha(inst_info,alphagrid,tgrid,alpha0,beta0,phib,options)
% Compute weights for numerical integral over alpha 
% arbitrary geometries 
% (response is integrated over beta)
% (response is also integrated over energy, except for E_TYPE=INT energy
% channels, for which the energy bandwith is assumed to be 1 to avoid inf)

[halpha,result_code] = inst_info.make_halpha(inst_info,alphagrid,tgrid,alpha0,beta0,phib,options);

