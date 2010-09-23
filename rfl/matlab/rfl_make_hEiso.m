function [hEiso,result_code] = rfl_make_hEiso(inst_info,Egrid,tgrid,options)
% [hE,result_code] = rfl_make_hEiso(inst_info,Egrid,tgrid,options)
% Compute weights for numerical integral over energy

[hEiso,result_code] = inst_info.make_hEiso(inst_info,Egrid,tgrid,options);


