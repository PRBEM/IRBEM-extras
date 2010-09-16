function resp = rfl_init_WIDE(resp)
% populate methods and variables for
% make_hE part of response function with
% idealized wide differential energy channel

resp.make_hE = @rfl_make_hE_wide;

function [hE,result_code] = rfl_make_hE_wide(resp,Egrid,options)
% treat as difference between two integral channels
% at E1 and E0
result_code = 1;
resptmp = resp;
resptmp.E0 = resptmp.E1;
hE = rfl_make_hE_int(resp,Egrid,options)-rfl_make_hE_int(resptmp,Egrid,options);
hE = hE/(resp.E1-resp.E0);
