function [ K Bm Xmeq Ymeq Zmeq Bmeq Xeq Yeq Zeq Beq Xfoot Yfoot Zfoot ...
    Xfl Yfl Zfl] = ...
    fieldline( x0, y0, z0, datenum_, ioptparmod, ...
    external, internal, ionoR, ds, n_threads)
%UBK.FIELDLINE Field Line Tracing and Bm(K) Evaluator
%   [ K Bm Xmeq Ymeq Zmeq Bmeq Xeq Yeq Zeq Beq Xfoot Yfoot Zfoot ...
%   Xfl Yfl Zfl] = ...
%   fieldline( x0, y0, z0, date_, ioptparmod, ...
%   external, internal, ionoR, ds, n_threads)
%   Trace magnetic field lines (bidirectional) and evaluate the modified
%   2nd invariant, K (Roederer 1970).
%   The tracing uses fixed step size RK4 (Numerical Recipe, Press).
%
%   References: Sheldon and Gaffey 1993;
%
%   REQUIRED INPUTS (All in SM coordinate system)
%   * [x0 y0 z0]: Matrices of SM cartesian coordinates in RE. M by N size.
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
%   * ionoR: Ionospheric boundary (radial distance to the ionosphere) in RE.
%   Default is 1.015 RE (~ 100 km above Earth).
%   * ds: Field line integration step size. Default is 0.05 RE.
%   * n_threads: The number of concurrent executions (positive integer).
%   Default is 8.
%
%   OUTPUTS (M by N)
%   * [K Bm]: Cell matrices of modified 2nd invariant and magnetic mirror
%   magnitude. Each element holds a vector of values evaluated at the
%   discrete points long the field line. K in nT^.5 RE, Bm in nT.
%   * [Xmeq Ymeq Zmeq]: SM magnetic equator coordinate matrices. The magnetic
%   equator is the location where the field strength is minimum along the
%   field line in RE.
%   * Bmeq: Matrix of magnetic field magnitude at the magnetic equator in
%   nT.
%   * [Xeq Yeq Zeq]: SM coordinates of intersection that the field line
%   crosses on the x-y plane in RE.
%   * Beq: Matrix of magnetic field magnitude at the coordinate equator in
%   nT.
%   * [Xfoot Yfoot Zfoot]: Matrices of the magnetic foot points of field
%   lines. SM cartesian in RE.
%   * [Xfl Yfl Zfl]: Cell matrices of each field line coordinate. RE.
%
%   NOTE
%   K and Bm are only evaluated from Bmin to either northern ionosphere if
%   Bmin is on the northern hemisphere or southern ionosphere if Bmin is on
%   the southern hemisphere.
%   The vectors in each cell element of K, Bm, Xfl, Yfl and Zfl have the
%   same length.

%
% $Author$
% $LastChangedDate$
% $Revision$
% $Id$
%

%% Argument check
error(nargchk(7,10,nargin))

MxN = size(x0);

% Size check
if ~isequal(size(x0), size(y0), size(z0), [MxN(1) length(datenum_)])
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


% Ionospheric boundary check
if nargin<8 || isempty(ionoR)
    ionoR = 1.015;
end
if ionoR < 1
    error('cotrans:InvalidArgument',...
        'ionoR should be greater than or equal to 1.')
end

% Integration step size check
if nargin<9 || isempty(ds)
    ds= .05;
end
if ds <= 0
    error('cotrans:InvalidArgument',...
        'step_size should be greater than 0.')
end

% n_thread check
if nargin<10 || isempty(n_threads)
    n_threads = 8;
end
if n_threads < 1
    warning('cotrans:InvalidArgument',...
        'n_threads should be greater than or equal to 1. Set to default value.')
    n_threads = 8;
end

%% Call MEX entry
[ K Bm ...
    Xmeq Ymeq Zmeq Bmeq ...
    Xeq Yeq Zeq Beq ...
    Xfoot Yfoot Zfoot ...
    Xfl Yfl Zfl] =...
    UBKFieldLinecxx(x0, y0, z0, ...
    [Y(:), DOY(:), H(:), MN(:), S(:)]', ...
    ioptparmod, external, internal, ...
    ionoR, ds, n_threads);

%% Post-processing

end