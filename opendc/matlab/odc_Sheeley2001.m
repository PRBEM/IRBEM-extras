function N0 = odc_Sheeley2001(which,L,MLT,apply_limits)
% N0 = odc_Sheeley2001(which,L,MLT)
% Return density model from Sheeley 2001
% which is 'plasmasphere' or 'trough'
% L is L shell in Re (Probably Lm in IGRF or OPQ)
% MLT is Magnetic local time in hours (LT in original paper)
% N0 is in #/cm^3
% outside 3 <= L <= 7, returns NaN unless...
% N0 = Sheeley2001(which,L,MLT,false) # ignor limits

if nargin < 4,
    apply_limits = true;
end

switch(lower(which)),
    case {'p','plasmasphere'},
        N0 = 1390*(3./L).^4.83; % n_p in paper
    case {'t','trough'},
        N0 = 124*(3./L).^4.0 + 36*(3./L).^3.5.*cos((MLT-(7.7*(3./L).^2.0+12))*pi/12); % n_t in paper
    otherwise
        error('Unknown "which" = "%s"',which);
end

if apply_limits,
    N0((L<3) | (L>7)) = NaN; % apply L limits from paper;
end

