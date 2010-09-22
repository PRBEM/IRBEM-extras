function [hEalphabeta,result_code] = rfl_make_hEalphabeta(inst_info,Egrid,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options)
% function [hEalphabeta,result_code] = rfl_make_hEalphabeta(inst_info,Egrid,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options)
% Compute weights for triple numerical integral over E, alpha, beta,
% arbitrary geometries 

[hEalphabeta,result_code] = inst_info.make_hEalphabeta(inst_info,Egrid,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options);

