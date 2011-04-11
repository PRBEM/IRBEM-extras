function maglat = BB0toMagLat(BB0,units)
% maglat = BB0toMagLat(BB0) % returns maglat in deg
% maglat = BB0toMagLat(BB0,'rad') % returns maglat in radians
%
% BB0 is B / Bmin
% where B is the local magnetic field strength
% and Bmin is the minimum (i.e., equatorial) magnetic field strength on the
% same field line.
%
% maglat is the unsigned dipole latitude

% we invert the B/Bmin vs maglat equation via a look-up table
% we use a persistent variable so we don't have to regenerate the table
% for each call to this function
persistent maglat_table
if isempty(maglat_table),
    maglat_table.deg = (0:0.1:90)';
    maglat_table.rads = maglat_table.deg*pi/180;
    % this next bit is the expression for B/Bmin vs magnetic latitude for a dipole
    maglat_table.bb0 = (1+3*sin(maglat_table.rads).^2).^(1/2)./cos(maglat_table.rads).^6;
end

if nargin < 2,
    units = 'deg';
end
switch(lower(units(1))),
    case 'd',
        maglat = interp1(maglat_table.bb0,maglat_table.deg,BB0,'linear');
    case 'r',
        maglat = interp1(maglat_table.bb0,maglat_table.rads,BB0,'linear');
    otherwise
        error('Unknown units "%s"',units);
end
