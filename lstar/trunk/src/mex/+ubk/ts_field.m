function [ bx by bz ] = ts_field( x, y, z, datenum_, ioptparmod, ...
    external, internal, varargin)
%UBK.TS_FIELD Tsynenko Magnetic Field
%   [ bx by bz ] = ts_field( x, y, z, datenum_, ioptparmod, ...
%       external, internal, param1, value1, ...)
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
%   IF external=='T89', THEN: ioptparmod = IOPT (length N), where
%       T89 input equivalent to Kp.
%       IOPT= 1       2        3        4        5        6        7
%       KP  = 0,0+  1-,1,1+  2-,2,2+  3-,3,3+  4-,4,4+  5-,5,5+  >=6-
%   ELSE IF external!='T89' && external!='NONE', THEN: ioptparmod = PARMOD,
%   where
%       PARMOD=[PDYN (nPa), DST (nT), BYIMF, BZIMF (nT), W1 W2 W3 W4 W5 W6]
%       The dimension should be [10, N].
%   * external: External field model. 'NONE' for no external model, 'T89' for
%   T89 model, 'T96' for T96 model, 'T02' for T02 model, or 'TS05' for TS05
%   model. Default is 'NONE'.
%   * internal: Internal field model. 'DIP' for dipole field or 'IGRF' for
%   IGRF model (case-insensitive). Default is 'IGRF'.
%
%   OPTIONAL INPUTS
%   * CO_SYSTEM: 'GSM' for GSM or 'SM' for SM coordinate system.
%   Case-insentive. Default is 'GSM'.
%   * M_THREADS: The number of threads to loop over the first dimension
%   (positive integer). Default is 8.
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
error(nargchk(7,nargin,nargin))

%% Required
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
        'External model %s is unknown.', external)
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
        'Internal model %s is unknown.', internal)
end

%% Options
opts = ubk.optset(varargin{:});

% Coordinate system check
if isempty(opts.CO_SYSTEM)
    opts.CO_SYSTEM = 'gsm';
end
if ~ischar(opts.CO_SYSTEM)
    error('cotrans:InvalidArgument',...
        'CO_SYSTEM is not character array.')
end
switch lower(opts.CO_SYSTEM)
    case 'sm'
        opts.CO_SYSTEM = 1;
    case 'gsm'
        opts.CO_SYSTEM = 0;
    otherwise
        error('cotrans:InvalidArgument',...
        'Coordinate system %s is unknown.', opts.CO_SYSTEM)
end

% M_threads check
if isempty(opts.M_THREADS)
    opts.M_THREADS = 8;
end
if opts.M_THREADS < 1
    warning('cotrans:InvalidArgument',...
        'M_THREADS should be greater than or equal to 1. Set to default value.')
    opts.M_THREADS = 8;
end

%% Call MEX entry
[bx by bz] = ...
    UBKTSFieldcxx(double(x), double(y), double(z), ...
    [Y(:), DOY(:), H(:), MN(:), S(:)]', ...
    double(ioptparmod), external, internal, opts.CO_SYSTEM, ...
    opts.M_THREADS, 1);

%% Post-processing

end
