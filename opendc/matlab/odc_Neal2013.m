function Lc = odc_Neal2013(channel,varName,varValue)
% Lc = odc_Neal2013(channel,varName,varValue)
% return L cutoff from Neal et al., 2013
% channel is 1, 2, or 3
% correpsonding to 24.3 MeV, 51.5 MeV, and 101.0 MeV at POES
%  (or 37.2, 76.7, and 151 MeV at 100 km alt)
% varName is 'Dst' or 'Kp'
% varValue is the value of Dst or Kp
% Note Lc is IGRF L value
% Also, Kp should be shifted to be 3 hours old, since apparently it
%   preceds the change in the cutoff

% Energy 
% at 100km
% (MeV) Indices A           B           C       r       R^2     Delta   Data Points
% 37.2  Kp      0.057912    0.38237     63.1626 0.70679 0.50154 1.72    7683
% 76.7  Kp      0.08087     0.14163     61.712  0.74373 0.6216  1.3243  4620
% 151   Kp      0.083756    0.06691     59.8825 0.80126 0.71039 1.1081  4547
% 37.2  Dst                 0.031679    62.5344 0.62563 0.46114 1.7938  7791
% 76.7  Dst                 0.029931    61.3043 0.68625 0.54862 1.4467  4653
% 151   Dst                 0.028514    59.5979 0.73982 0.64016 1.2349  4581

DATA = [%  var=1=Kp, var=2=Dst
% channel, var, A,          B,          C
1          1    0.057912    0.38237     63.1626
2          1    0.08087     0.14163     61.712
3          1    0.083756    0.06691     59.8825
1          2    0           0.031679    62.5344
2          2    0           0.029931    61.3043
3          2    0           0.028514    59.5979
];


if lower(varName(1)) == 'k',
    ivar = 1;
else
    ivar = 2;
end

irow = find((DATA(:,1)==channel) & (DATA(:,2)==ivar));
A = DATA(irow,3);
B = DATA(irow,4);
C = DATA(irow,5);

% IGRF invlat = A*var^2+B*var+C

invlat = A*varValue.^2 + B*varValue+C;
Lc = cosd(invlat).^-2;
