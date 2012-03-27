function [PSDM,PSDE] = odc_PSD_Ozeke2012(param,value,L,f)
% [PSDM,PSDE] = odc_PSD_Ozeke2012(param,value,L,f)
% returns power spectral density
% for magnetic and electric fields
% from:
% Ozeke, L. G., I. R. Mann, K. R. Murphy, I. J. Rae, 
% D. K. Milling, S. R. Elkington, A. A. A. Chan, and H. J. Singer (2012),
% ULF Wave Derived Radiation Belt Radial Diffusion Coefficients,
% J. Geophys. Res., doi:10.1029/2011JA017463, in press. 
% (accepted 23 February 2012) 
%
% PSDM units (nT)^2/mHz
% PSDE units (mV/m)^2/mHz
% param - 'Kp' or 'Vsw'
% value - Kp or Vsw (in km/s)
% L - L shell
% f - frequency, mHz
% only f can be non-scalar
% Interpolation/Extrapolation notes:
%  For PSDM L>6.5, GEO values are used with no L scaling
%  For PSDM 2.5<=L<6.5, is constant within a bin, with no L interpolation
%  For both PSDM and PSDE, Kp and Vsw are constant within a bin
%  For PSDM, the highest Kp bin of 6 is applied for all Kp>=6
%  For PSDE, interpolation in L is provided for 2.55 <= L <= 7.94.
%    Extrapolation is flat.
%
% run with no input arguments to generate test figures from paper
%
% requires two data tables in search path:
%  Ozeke2012_PSDvsKp.csv, and Ozeke2012_PSDvsVsw.csv

persistent PSDdata

if nargin == 0, % run test case
    PSDM = [];
    PSDE = [];
    test;
    return
end

if ~ismember(param,{'Kp','Vsw'}),
    error('param="%s" unrecognized. Expected "Kp" or "Vsw"',param);
end

if ~isfield(PSDdata,param),
    PSDdata.(param) = load_coefs(param);
end

% do PSDE
coefs = PSDdata.(param);
xbin = find((value>=coefs.xmin) & (value <= coefs.xmax));
if isempty(xbin),
    error('%s=%g out of range',param,value);
end
xbin = xbin(1); % chose lower bin to achieve a <= x < b
% find lower and upper L bins
iL1 = max(1,sum(coefs.E.L<=L));
iL2 = min(iL1 + 1,length(coefs.E.L));

% q = log10(PSD/f^2) = a3*exp(-z^2/2)+a2+a1*y, y = log10(f), z = (y-a4)/a5
y = log10(f);
z = (y-coefs.E.a4(iL1,xbin))./coefs.E.a5(iL1,xbin);
% q at L1
q1 = coefs.E.a3(iL1,xbin)*exp(-z.^2/2)+coefs.E.a2(iL1,xbin)+coefs.E.a1(iL1,xbin).*y;
if iL2 == iL1, % no need to interpolate
    q = q1;
else
    % compute q at L2
    z = (y-coefs.E.a4(iL2,xbin))./coefs.E.a5(iL2,xbin);
    q2 = coefs.E.a3(iL2,xbin)*exp(-z.^2/2)+coefs.E.a2(iL2,xbin)+coefs.E.a1(iL2,xbin).*y;
    % interpolate
    C = (L-coefs.E.L(iL1))/(coefs.E.L(iL2)-coefs.E.L(iL1));
    q = q2*C+(1-C)*q1;
end
PSDE = (10.^q).*(f.^2); % undo log(PSD/f^2)

% do PSDM
xbin = min(xbin,length(coefs.M.A)); % Kp only goes to 6, not 9. Otherwise, same bins

% PSD = A*f^-B @ GEO
iL = find((L>=coefs.M.Lmin) & (L <= coefs.M.Lmax));
if L>coefs.M.Lmax(end), % use GEO values
    PSDM = coefs.M.A(xbin).*f.^-coefs.M.B(xbin); % nT^2/mHz
elseif L<coefs.M.Lmin(1),
    warning('L<%g, PSDM is NaN',coefs.M.Lmin(1));
    PSDM = nan(size(PSDE));
else
    iL = iL(1); % chose lower bin to achieve a <= L < b
    % compute F, X from ratios in Tables 3 and 4
    F = coefs.M.A(xbin)/coefs.M.Aref;
    X = coefs.M.B(xbin)/coefs.M.Bref;
    % evaulate expression from Tables 3 and 4
    PSDM = F.*coefs.M.prefactor(iL).*f.^(-coefs.M.exponent(iL).*X); % nT^2/Hz
    % convert from /Hz to /mHz
    PSDM = PSDM/1e3; % nT^2/mHz
end


function coefs = load_coefs(param)
% read coefficients for Kp or Vsw
% Ozeke2012_PSDvsKp.csv
% Ozeke2012_PSDvsVsw.csv

filename = which(['Ozeke2012_PSDvs',param,'.csv']); % search path
if ~exist(filename,'file'),
    error('Unable to locate %s',filename);
end
data = csvread(filename);

coefs = struct('param',param);

if isequal(param,'Kp'),
    x = data(1,3:end); % remove left "headers"
    % construct Kp bin boundaries
    coefs.xmin = max(0,x-0.5);
    coefs.xmax = min(9+1/3,x+0.5);
    % construct labels
    coefs.labels = arrayfun(@(x)sprintf('Kp~%g',x),x,'uniform',false);
    % remove header row
    data = data(2:end,:);
else
    % construct Vsw bins
    coefs.xmin = data(1,3:end); % remove left "headers"
    coefs.xmax = data(2,3:end); % remove left "headers"
    coefs.xmax(end) = inf; % no upper bound on last bin
    % construct labels
    coefs.labels = cell(1,length(coefs.xmin));
    for i = 1:(length(coefs.xmin)-1),
        coefs.labels{i} = sprintf('%g <= Vsw < %g km/s',coefs.xmin(i),coefs.xmax(i));
    end
    if coefs.xmin(1)==0,
        coefs.labels{1} = sprintf('Vsw < %g km/s',coefs.xmax(1));
    end
    coefs.labels{end} = sprintf('Vsw > %g km/s',coefs.xmin(end));
    % remove header rows
    data = data(3:end,:);
end

[coefs.E.L,iLsort] = sort(data(1:5:end,1));

for ia = 1:5,
    tmp = data(ia:5:end,3:end);
    coefs.E.(sprintf('a%d',ia)) = tmp(iLsort,:);
end

coefs.M.Lmin = (2.5:5.5)';
coefs.M.Lmax = (3.5:6.5)';
coefs.M.L = (4:6)';

switch(param),
    case 'Kp',
        % Table 1
        coefs.M.A = [0.02	0.054	0.163	0.445	1.021	2.018	4.293]';
        coefs.M.B = [2.37	2.54	2.4	2.31	2.24	2.18	2.06]';
        coefs.M.Aref = coefs.M.A(3); % Kp==2
        coefs.M.Bref = coefs.M.B(3); % Kp==2
        % Table 3
        coefs.M.prefactor = 10.^[1;1.2;1.3;1.5];
        coefs.M.exponent = [1.3;1.9;1.9;1.9];
    case 'Vsw',
        % Table 2
        coefs.M.A = [0.031	0.074	0.168	0.38	0.671	1.292]';
        coefs.M.B = [2.47	2.45	2.43	2.32	2.22	2.1]';
        coefs.M.Aref = coefs.M.A(3); % Vsw ~ 450
        coefs.M.Bref = coefs.M.B(3); % Vsw ~ 450
        % Table 4
        coefs.M.prefactor = 10.^[1;1.2;1.3;1.5];
        coefs.M.exponent = [1.3;1.9;1.9;1.9];
end

function test

% make PSD figures from paper

figure(4);
sp = [1,3,5,7,4,6,8];
Ls = [2.55,2.98,4.21,4.26,5.4,6.51,7.94];
Kps = [1,3,6];
styles = {'go-','b^-','rs-'};
f = (0.5:0.2:20)';
for i = 1:length(sp),
    subplot(4,2,sp(i));
    L = Ls(i);
    leg = {};
    for j = 1:length(Kps),
        Kp = Kps(j);
        [PSDM,PSDE] = odc_PSD_Ozeke2012('Kp',Kp,L,f);
        loglog(f,PSDE,styles{j});
        leg{j} = sprintf('Kp=%g',Kp);
        hold on;
    end
    title(sprintf('L=%g',L));
    axis([0.1 20 1e-4 1]);
    if i==1,
        legend(leg{:},'location','sw');
    end
    if rem(sp(i),2) == 1,
        ylabel('(mV/m)^2/mHz');
    end
    grid on
end

figure(5); % panel c
f = (1:0.1:9)';
Kps = 0:6;
L = 6.6;
leg = {};
styles = {'--','+-','^-','d-','*-','o-','s-'};
for j = 1:length(Kps),
    Kp = Kps(j);
    [PSDM,PSDE] = odc_PSD_Ozeke2012('Kp',Kp,L,f);
    loglog(f,PSDM,styles{j});
    leg{j} = sprintf('Kp=%g',Kp);
    hold on;
end
title(sprintf('L=%g',L));
axis([1 10 1e-4 1e1]);
legend(leg{:},'location','sw');
ylabel('(nT)^2/mHz');
grid on


figure(6); % panel c
f = (1:0.1:9)';
Vsws = [150,350:100:750];
L = 6.6;
leg = {};
styles = {'--','+-','^-','d-','*-','o-','s-'};
for j = 1:length(Vsws),
    Vsw = Vsws(j);
    [PSDM,PSDE] = odc_PSD_Ozeke2012('Vsw',Vsw,L,f);
    loglog(f,PSDM,styles{j});
    leg{j} = sprintf('Vsw=%g',Vsw);
    hold on;
end
title(sprintf('L=%g',L));
axis([1 10 1e-4 1e1]);
legend(leg{:},'location','sw');
ylabel('(nT)^2/mHz');
grid on

figure (8); % panels a and b
Kps = [2,4,6];
Ls = [2.55,2.98,4.21,5.4,6.51,7.94];
f = [2.5;6.11;8.06];
styles = {'bd-','g^-','rs-'};
leg = cell(1,3);
for i = 1:length(Kps),
    Kp = Kps(i);
    subplot(3,2,i*2-1);
    y = zeros(length(Ls),length(f));
    for j = 1:length(Ls),
        [PSDM,PSDE] = odc_PSD_Ozeke2012('Kp',Kp,Ls(j),f);
        y(j,:) = PSDE;
    end
    for j = 1:length(f),
        semilogy(Ls,y(:,j),styles{j});
        hold on;
        leg{j} = sprintf('%.2f mHz',f(j));
    end
    xlabel('L Shell');
    ylabel(sprintf('(mV/m)^2/mHz, Kp=%g',Kp));
    axis([2 8 0.5e-4 2e1]);
    grid on;
    if i == 1,
        legend(leg{:},'location','nw');
    end
end

Kps = 0:8;
Vsws = [350,379,415,463,520,563,564,558,662];
Ls = [4.22,5.4,6.51];
styles = {'bd-','g^-','rs-'};
for i = 1:length(Ls),
    L = Ls(i);
    subplot(3,2,i*2);
    y = zeros(length(Vsws),length(f));
    for j = 1:length(Vsws),
        [PSDM,PSDE] = odc_PSD_Ozeke2012('Kp',Kps(j),L,f);
        y(j,:) = PSDE;
    end
    for j = 1:length(f),
        semilogy(Kps,y(:,j),styles{j});
        hold on;
    end
    xlabel('Kp');
    ylabel(sprintf('(mV/m)^2/mHz, L=%g',L));
    axis([0 9 0.5e-4 2e1]);
    grid on;
end



