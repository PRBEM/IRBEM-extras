function mat2cdf(matfile,cdffile)
% mat2cdf(matfile,cdffile);
% mat2cdf(matfile) - guesses cdffile name
% convert .mat file to .cdf
% the .cdf file can be read by cdf2var

if nargin < 2,
    cdffile = strrep(matfile,'.mat','.cdf'); % try to change .mat to .cdf
    if strcmp(matfile,cdffile),
        cdffile = [matfile,'.cdf']; % just append .cdf to end of matfile
    end
end
if ~exist(matfile,'file'),
    error('file not found: "%s"',matfile);
end

s = whos('-file',matfile);
if length(s) == 1,
    var = load(matfile);
    var = var.(s(1).name);
else
    var = load(matfile);
end

var2cdf(cdffile,var);
