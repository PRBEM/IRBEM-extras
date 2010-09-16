function resp = rfl_init_DIFF(resp)
% populate methods and variables for
% make_hE part of response function with
% idealized differential energy channel

resp.make_hE = @rfl_make_hE_diff;

function [hE,result_code] = rfl_make_hE_diff(resp,Egrid,options)
result_code = 1;
hE = zeros(size(Egrid));
i = rfl_get_list_neighbors(Egrid,resp.E0);
if all(i>0),
    hE(i(1)) = (Egrid(i(2))-resp.E0);
    hE(i(2)) = (resp.E0-Egrid(i(1)));
    hE(i) = hE(i)/diff(Egrid(i))*resp.DE*resp.EPS;
end
