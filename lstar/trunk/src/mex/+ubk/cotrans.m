function [ xout yout zout psi ] = cotrans(xin, yin, zin, datenum_, from_to, ...
    M_threads)
%UBK.COTRANS GSM <--> SM
%   [ xout yout zout psi ] = cotrans(xin, yin, zin, datenum_, from_to, ...
%   M_threads)
%   Coordinate-transform between GSM and SM coordinate at a given date
%   using GEOPACK library (part of Tsyganenko model,
%   http://geo.phys.spbu.ru/~tsyganenko/modeling.html).
%
%   REQUIRED INPUTS
%   * [xin yin zin]: Matrices of cartesian coordinates in RE. M by N size.
%   * datenum_: N-element time vector returned by datenum.
%   * from_to: 'GSM2SM' for GSM-->SM or 'SM2GSM' for SM-->GSM.
%   Case-insensitive.
%
%   OPTIONAL INPUTS (pass [] to use default)
%   * M_threads: The number of threads to loop over the first dimension
%   (positive integer). Default is 8.
%
%   OUTPUTS
%   * [xout yout zout]: Transformed coordinates with size M by N.
%   * psi: Dipole tilt angle in radian obtained from GEOPACK. Size 1 by N.

%
% $Author$
% $LastChangedDate$
% $Revision$
% $Id$
%

%% Argument check
error(nargchk(5,6,nargin))

MxN = size(xin);

% Size check
if ~isequal(size(xin), size(yin), size(zin), [MxN(1) length(datenum_)])
    error('cotrans:InvalidArgument',...
        'size of xin, yin and zin is different.')
end

% Date check
if ~isnumeric(datenum_)
    error('cotrans:InvalidArgument',...
        'Non-numeric datenum_ input.')
end
if length(datenum_)~=numel(datenum_)
    error('cotrans:InvalidArgument',...
        'Invalid datenum_ dimension.')
end
[Y, ~, ~, H, MN, S] = datevec(datenum_);
tmp = datenum([Y(:), zeros(numel(Y),1), ones(numel(Y),1)]);
DOY = floor( datenum_ - reshape(tmp, size(datenum_)) ) + 1;

% from_to check
if ~ischar(from_to)
    error('cotrans:InvalidArgument',...
        'from_to is not character array.')
end
switch lower(from_to)
    case 'gsm2sm'
        to_co_system = 1;
    case 'sm2gsm'
        to_co_system = 0;
    otherwise
        error('cotrans:InvalidArgument',...
        '%s is unknown.', from_to)
end

% M_threads check
if nargin<6 || isempty(M_threads)
    M_threads = 8;
end
if M_threads < 1
    warning('cotrans:InvalidArgument',...
        'M_threads should be greater than or equal to 1. Set to default value.')
    M_threads = 8;
end

%% Call MEX entry
[xout yout zout psi] = ...
    UBKGeopackCotranscxx(xin, yin, zin, ...
    [Y(:), DOY(:), H(:), MN(:), S(:)]', ...
    to_co_system, M_threads, 1);

%% Post-processing

end
