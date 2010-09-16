function i = rfl_get_list_neighbors(list,q)
% returns i = Nx2, where list(i(:,1))<=q, and list(i(:,2)>=q
% i==0 implies q out of range

N = length(list);
i = nan(length(q),2);
q = q(:);
iq = find(isfinite(q));
list = list(:);
i(iq,1) = interp1(list,(1:N)',q(iq),'nearest');
f = isfinite(i(:,1)); % nearest can return nan for out of range
if any(f),
    i(f,1) = i(f,1) - (q(f)<list(i(f,1))); % force to greatest list <= q
    i(f,2) = i(f,1)+1; % next neighbor to right
end
% limit checks
i(i > N) = 0;
i(i < 1) = 0;
i(~isfinite(i)) = 0;

