function [hEalpha,result_code] = rfl_make_hEalpha(inst_info,Egrid,alphagrid,tgrid,alpha0,beta0,phib,options)
% function [hEalpha,result_code] = rfl_make_hEalpha(inst_info,Egrid,alphagrid,tgrid,alpha0,beta0,phib,options)
% Compute weights for double numerical integral over E, alpha, 
% arbitrary geometries 

[hEalpha,result_code] = inst_info.make_hEalpha(inst_info,Egrid,alphagrid,tgrid,alpha0,beta0,phib,options);

