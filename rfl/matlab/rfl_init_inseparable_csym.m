function inst_info = rfl_init_inseparable_csym(inst_info)
% inst_info = rfl_init_inseparable_csym(inst_info)
% initialize cylindrically symmetric sensor with inseparable E,theta response
% RESP_TYPE = '[E,TH]';
% inst_info is assumed to be a fully-initialized inseparable R(E,theta,phi) sensor
% this function just overloads the definition of R

% overload
% generic response function in E, theta coordinates, ignores phi input
inst_info.R = @(inst_info,E,theta,phi)interpn(inst_info.E_GRID,inst_info.TH_GRID,inst_info.R,E,theta,'linear',0)/inst_info.XCAL; % set to zero outside grid

% All other methods are the same as in the R(E,theta,phi) inseparable case.
