function resp = rfl_init_INT(resp)
% populate methods and variables for
% make_hE part of response function with
% idealized integral energy channel

resp.make_hE = @rfl_make_hE_int;

function [hE,result_code] = rfl_make_hE_int(resp,Egrid,options)
result_code = 1;
hE = zeros(size(Egrid));

i = find(resp.E0 <= Egrid);
if ~isempty(i),
    dE = rfl_make_deltas(Egrid,options);
    if (i(1)==1),
        % all channels fully involved
        hE = dE;
    else
        hE(i) = dE(i); % all i channels fully involved
        Em = Egrid(i(1)-1); % Eminus
        Ep = Egrid(i(1)); % Eplus
        hE(i(1)-1) = (Ep-resp.E0)^2/(2*(Ep - Em));
        hE(i(1)) = ((Ep-resp.E0)*(resp.E0  + Ep - 2*Em))/(2*(Ep - Em));
    end
    hE = hE*resp.EPS;
end

