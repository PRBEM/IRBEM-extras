function [hE,result_code] = rfl_make_hE(inst_info,Egrid,options)
% [hE,result_code] = rfl_make_hE(inst_info,Egrid,options)
% Compute weights for numerical integral over energy

[hE,result_code] = inst_info.make_hE(inst_info,Egrid,options);

