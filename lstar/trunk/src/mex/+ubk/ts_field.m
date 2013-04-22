function [ bx by bz ] = ts_field( x, y, z, datenum_, ioptparmod, ...
    external, internal, co_system, n_threads)
%UBK.TS_FIELD Tsynenko Magnetic Field
%   [ bx by bz ] = ts_field( x, y, z, datenum_, ioptparmod, ...
%   external, internal, co_system, n_threads)
%   Calculates magnetic field using Tsyganenko magnetic field models.
%   Please refer to http://geo.phys.spbu.ru/~tsyganenko/modeling.html and
%   references therein.
%
%   References: Tsyganenko 1987, 1989, 1995, 2002a, 2002b;
%   Tsyganenko and Sitnov 2005;
%
%   REQUIRED INPUTS
%   * [x y z]: Matrices of cartesian coordinates in RE. M by N size.
%   * datenum_: N-element time vector returned by datenum.
%   * ioptparmod: External field parameters.
%   IF external=='T89', THEN: ioptparmod = IOPT (length N)
%       T89 input equivalent to Kp.
%       IOPT= 1       2        3        4        5        6        7
%       KP  = 0,0+  1-,1,1+  2-,2,2+  3-,3,3+  4-,4,4+  5-,5,5+  >=6-
%   ELSE IF external!='T89' && external!='NONE', THEN: ioptparmod = PARMOD
%       PARMOD=[PDYN (nPa), DST (nT), BYIMF, BZIMF (nT), W1 W2 W3 W4 W5 W6]
%       The dimension should be [10, N].
%   * external: External field model. 'NONE' for no external model, 'T89' for
%   T89 model, 'T96' for T96 model, 'T02' for T02 model, or 'TS05' for TS05
%   model. Default is 'NONE'.
%   * internal: Internal field model. 'DIP' for dipole field or 'IGRF' for
%   IGRF model (case-insensitive). Default is 'IGRF'.
%
%   OPTIONAL INPUTS (pass [] to use default)
%   * co_system: 'GSM' for GSM or 'SM' for SM coordinate system.
%   Case-insentive. Default is 'GSM'.
%   * n_threads: The number of concurrent executions (positive integer).
%   Default is 8.
%
%   OUTPUTS
%   * [bx by bz]: Calculated magnetic field vectors of size M by N.

%
% $Author$
% $LastChangedDate$
% $Revision$
% $Id$
%

%% Argument check
error(nargchk(7,9,nargin))

MxN = size(x);

% Size check
if ~isequal(size(x), size(y), size(z), [MxN(1) length(datenum_)])
    error('cotrans:InvalidArgument',...
        'size of x, y and z is different.')
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

% External field
if ~ischar(external)
    error('cotrans:InvalidArgument',...
        'external is not character array.')
end
switch lower(external)
    case 'none'
        external = 0;
    case 't89'
        external = 1;
    case 't96'
        external = 2;
    case 't02'
        external = 3;
    case 'ts05'
        external = 4;
    otherwise
        error('cotrans:InvalidArgument',...
        'Coordinate system %s is unknown.', external)
end

if external==1
    if length(ioptparmod)~=MxN(2)
        error('cotrans:InvalidArgument',...
            'Incompatible IOPT length.')
    end
    ioptparmod = ioptparmod(:)';
elseif (external~=1&&external~=0)
    if ~isequal(size(ioptparmod),[10 MxN(2)])
        error('cotrans:InvalidArgument',...
            'Incompatible PARMOD size.')
    end
end

% Internal field
if ~ischar(internal)
    error('cotrans:InvalidArgument',...
        'internal is not character array.')
end
switch lower(internal)
    case 'dip'
        internal = 0;
    case 'igrf'
        internal = 1;
    otherwise
        error('cotrans:InvalidArgument',...
        'Coordinate system %s is unknown.', internal)
end

% Coordinate system check
if nargin<8 || isempty(co_system)
    co_system = 'gsm';
end
if ~ischar(co_system)
    error('cotrans:InvalidArgument',...
        'co_system is not character array.')
end
switch lower(co_system)
    case 'sm'
        co_system = 1;
    case 'gsm'
        co_system = 0;
    otherwise
        error('cotrans:InvalidArgument',...
        'Coordinate system %s is unknown.', co_system)
end

% n_thread check
if nargin<9 || isempty(n_threads)
    n_threads = 8;
end
if n_threads < 1
    warning('cotrans:InvalidArgument',...
        'n_threads should be greater than or equal to 1. Set to default value.')
    n_threads = 8;
end

%% Call MEX entry
[bx by bz] = ...
    UBKTSFieldcxx(x, y, z, ...
    [Y(:), DOY(:), H(:), MN(:), S(:)]', ...
    ioptparmod, external, internal, co_system, ...
    n_threads);

%% Post-processing

end
