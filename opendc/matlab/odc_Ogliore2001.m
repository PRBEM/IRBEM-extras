function Lc = odc_Ogliore2001(MeV,species,Dst)
% L = odc_Ogliore2001(MeV,species)
% return L cutoff from Ogliore et al., 2001
% Ogliore 2001: Rc = 15.062/Lc^2 - 0.363
% L = odc_Ogliore2001(MeV,species,Dst)
% adds Dst correction
% from Leske et al. 2001 
% which is Dst/19.11 degrees of invariant latitude


% Ogliore 2001: Rc = 15.062/Lc^2 - 0.363 A/L^2 - B
% for Rc in GV
% (Rc+B)*L^2 = A
% L = sqrt(A/(Rc+B))

util = odc_util;
% R = util.rigidity(MeV,species) % returns rigidity (p/q) in GV
Lc_func = @(Rc)sqrt(15.062./(Rc+0.363));
Lc = Lc_func(util.rigidity(MeV,species));

if nargin >=3,
% Leske et al 2001
% add Dst/19.11 to invlat
    invlat = acosd(Lc.^-0.5)+Dst/19.11;
    Lc = cosd(invlat).^-2;
end