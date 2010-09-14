function [homni,result_code] = rfl_make_homni(inst_info,options)
% [homni,result_code] = rfl_make_homni(inst_info,options) Compute geometric
% factor omnidirectional instrument (replaces integral over theta/phi or
% alpha/beta)

if isfield(inst_info,'G'),
    result_code = 1;
    homni = inst_info.G;
else
    if nargout >= 2,
        result_code = -2;
        homni = [];
    else
        error('Requested omnidirectional approximation for channel with no "G" property');
    end
end
