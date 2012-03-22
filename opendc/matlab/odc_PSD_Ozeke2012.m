function [PSDM,PSDE] = odc_PSD_Ozeke2012(param,value,L,f)
% [PSDM,PSDE] = odc_PSD_Ozeke2012(param,value,L,f)
% returns power spectral density
% for magnetic and electric fields
% from Ozeke et al., 2012
% PSDM units (nT)^2/mHz
% PSDE units (mV/m)^2/mHz
% param - 'Kp' or 'Vsw'
% value - Kp or Vsw (in km/s)
% L - L shell
% f - frequency, mHz
% only f can be non-scalar

persistent PSDdata

if ~ismember(param,{'Kp','Vsw'}),
    error('param="%s" unrecognized. Expected "Kp" or "Vsw"',param);
end

if ~isfield(PSDdata,param),
    PSDdata.(param) = load_coefs(param);
end

coefs = PSDdata.(param);

% q = log10(PSD/f^2) = a3*exp(-z^2/2)+a2+a1*y, y = log10(f), z = (y-a3)/a5
y = log10(f);

xbin = find((value>=coefs.xmin) & (value <= coefs.xmax));
if isempty(xbin),
    error('%s=%g out of range',param,value);
end
xbin = xbin(1); % chose lower bin to achieve a <= x < b
iL1 = max(1,sum(coefs.L>=L));
iL2 = min(iL1 + 1,length(coefs.L));

error('working here');


function coefs = load_coefs(param)

% 'C:\work\CrossTerms\OzekeMann';
% Ozeke2012_PSDvsKp.csv
% Ozeke2012_PSDvsVsw.csv

filename = which(['Ozeke2012_PSDvs',param,'.csv']);
if ~exist(filename,'file'),
    error('Unable to locate %s',filename);
end
data = csvread(filename);

coefs = struct('param',param);

if isequal(param,'Kp'),
    x = data(1,3:end);
    coefs.xmin = max(0,x-0.5);
    coefs.xmax = min(9,x+0.5);
    coefs.labels = arrayfun(@(x)sprintf('Kp~%g',x),x,'uniform',false);
    data = data(2:end,:);
else
    coefs.xmin = data(1,3:end);
    coefs.xmax = data(2,3:end);
    coefs.labels = cell(1,length(coefs.xmin));
    for i = 1:(length(coefs.xmin)-1),
        coefs.labels{i} = sprintf('%g <= Vsw < %g km/s',coefs.xmin(i),coefs.xmax(i));
    end
    if coefs.xmin(1)==0,
        coefs.labels{1} = sprintf('Vsw < %g km/s',coefs.xmax(1));
    end
    coefs.labels{end} = sprintf('Vsw > %g km/s',coefs.xmin(end));
    data = data(3:end,:);
end

[NL,Nx] = size(data);
Nx = Nx-1;

coefs.L = data(1:5:end,1);

for ia = 1:5,
    coefs.(sprintf('a%d',ia)) = data(ia:5:end,3:end);
end

