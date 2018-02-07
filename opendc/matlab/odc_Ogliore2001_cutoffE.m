function E = odc_Ogliore2001_cutoffE(L,species,Dst)
% E = odc_Ogliore2001_cutoffE(L,species)
% return cutoff energy in MeV from Ogliore et al., 2001
% at specified L
% Ogliore 2001: Rc = 15.062/Lc^2 - 0.363
% E = odc_Ogliore2001_Cutoff(L,species,Dst)
% adds Dst correction
% from Leske et al. 2001 
% which is Dst/19.11 degrees of invariant latitude


% Ogliore 2001: Rc = 15.062/Lc^2 - 0.363 A/L^2 - B
% for Rc in GV
% (Rc+B)*L^2 = A
% L = sqrt(A/(Rc+B))

util = odc_util;

if nargin >= 3,
    % Leske et al 2001
    % adds Dst/19.11 to Ogliore result, reverse that here
    invlat = acosd(L.^-0.5)-Dst/19.11;
    L = cosd(invlat).^-2;
end

R = 15.062./L.^2 - 0.363; % GV
R(R<0) = nan;
R = R*1e9; % V
% R = p/q*c % kg m^2 /s^2 / C = V 
sp = util.SelectSpecies(species);
mks = util.mks;
p = R/mks.c*abs(mks.(sp).q);
% E = sqrt(m0^2*c^4 + p^2*c^2)-m0*c^2
E = sqrt(mks.(sp).m0^2*mks.c^4 + p.^2*mks.c^2)-mks.(sp).m0*mks.c^2; % J
E = E./(mks.eV*1e6); % MeV


