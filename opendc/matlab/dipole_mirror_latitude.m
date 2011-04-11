function mirror_lat = dipole_mirror_latitude(alpha0,units)
% mirror_lat = dipole_mirror_latitude(alpha0)
% compute dipolar mirror latitude of particle with equatorial pitch angle
% alpha0, angles in degrees
% mirror_lat = dipole_mirror_latitude(alpha0,'rad')
% angles in radians

if nargin < 2,
    units = [];
end

if isempty(units),
    units = 'deg';
end

if lower(units(1))=='d',
    torad = pi/180;
else    
    torad = 1;
end

sina0 = sin(alpha0*torad);
% note, below fixes error in Shprits thesis, eqn  F13
% which has sina0^2. Instead use sina0^4 from Shprits 2006 eqn 10.
Px = [1, 0, 0, 0, 0, 3*sina0^4, -4*sina0^4]; 
xroots = roots(Px);
xroot = xroots((imag(xroots)==0) & (xroots>0));
mirror_lat = acos(sqrt(xroot))/torad; % mirror latitude

% % check against Schulz & Lanzerotti, 1.25
% thetam = pi/2-mirror_lat*torad;
% y = (1+3*cos(thetam).^2).^(-1/4)*sin(thetam)^3;
% [y sina0] % should be equal

