function [ Ls K Phi0 XorRc YorPhic Phif Thetaf ] = ...
    lstar( x0, y0, pa0, datenum_, ioptparmod, external, internal, ...
    ionoR, ds, dXorR, dYorPhi, n_phi, n_theta, ...
    isCatesianGrid, M_threads, N_threads)
%UBK.LSTAR_CART Lstar in Cartesian Grid
%   [ Ls K Phi0 XorRc YorPhic Phif Thetaf ] = ...
%   lstar( x0, y0, pa0, datenum_, ioptparmod, external, internal, ...
%   ionoR, ds, dXorR, dYorPhi, n_phi, n_theta, ...
%   shouldKeepContour, isCatesianGrid, M_threads, N_threads)
%   Traces particle's drift contour and evaluates generalized L value, L*
%   (Roederer 1970), using built-in Tsyganenko field models.
%   SM coordinate system is default (and the only choice).
%
%   Please refer to http://geo.phys.spbu.ru/~tsyganenko/modeling.html and
%   references therein for Tsyganenko field models.
%
%   References: Tsyganenko 1987, 1989, 1995, 2002a, 2002b;
%   Tsyganenko and Sitnov 2005; Sheldon and Gaffey 1993; Min et al. 2012,
%   2013;
%
%   REQUIRED INPUTS
%   * [x0 y0]: Matrices of SM cartesian/cylindrical coordinates MAPPED ON
%   THE X-Y PLANE ALONG THE FIELD LINE. RE. M by N where N ==
%   length(datenum_).
%   * pa0: Initial pitch angle in radian AT THE MAGNETIC EQUATOR. M by N.
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
%   * ds: Field line integration step size. Default is 0.1 RE.
%   * dXorR: X/R grid resolution in RE. Default is 0.2 RE.
%   * dYorPhi: Y/Phi grid resolution in RE/radian. Y default is 0.2 RE and
%   Phi default is 2pi / n_phi.
%   * n_phi: The number of azimuthal grids for magnetic flux evaluation.
%   Default is 350/5.
%   * n_theta: The number of polar grids for magnetic flux evaluation.
%   Default is 180*2.
%   * isCartesianGrid: Boolean to indicate which grid scheme to use.
%   If true, cartesian grid scheme is used, else cylindrical grid shceme is
%   used. Input/output coordinate system depends on this setting except for
%   the foot point coordinates which is always in SM spherical coordinate
%   system. Default is false.
%   * M_threads: The number of threads to loop over the first dimension
%   (positive integer). Default is 8.
%   * N_threads: The number of threads to loop over the second dimension
%   (positive integer). Default is 2.
%
%   OUTPUTS (M by N otherwise specified)
%   * Ls: Matrix of L*.
%   * K: Matrix of modified 2nd invariant in nT^.5 RE.
%   * Phi0: N-element vector of Earth's magnetic moment times 2*pi.
%   Magnetic flux defined by the drift contour is Phi0/Ls (you have to
%   work out the sign). nT RE^2.
%   * [XorRc YorPhic]: Cell matrices of the drift contour coordinates.
%   If shouldKeepContour==false, cell elements are empty but the cell is
%   still M by N. RE and radian in SM coordinate system.
%   * [Phif Thetaf]: Cell matrices of the foot points of the drift
%   contour. If shouldKeepContour==false, cell elements are empty but
%   the cell is still M by N. RE and radian in SM spherical coordinate
%   system.

%
% $Author$
% $LastChangedDate$
% $Revision$
% $Id$
%

%% Argument check
error(nargchk(7,16,nargin))

MxN = size(x0);

% Size check
if ~isequal(size(x0), size(y0), size(pa0), [MxN(1) length(datenum_)])
    error('cotrans:InvalidArgument',...
        'size of x, y and pa is different.')
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
    ds= .1;
end
if ds <= 0
    error('cotrans:InvalidArgument',...
        'step_size should be greater than 0.')
end

% dx/dr check
if nargin<10 || isempty(dXorR)
    dXorR = .2;
end
if dXorR <= 0
    error('cotrans:InvalidArgument',...
        'Invalid dx.')
end

% n_phi check
if nargin<12 || isempty(n_phi)
    n_phi = 360/5;
end
if n_phi <= 10
    error('cotrans:InvalidArgument',...
        'Invalid n_phi.')
end

% isCatesianGrid check
if nargin<14 || isempty(isCatesianGrid)
    isCatesianGrid = false;
end

% dy/dphi check
if nargin<11 || isempty(dYorPhi)
    if isCatesianGrid
        dYorPhi = .2;
    else
        dYorPhi = 2*pi / n_phi;
    end
end
if dYorPhi <= 0
    error('cotrans:InvalidArgument',...
        'Invalid dy.')
end

% n_theta check
if nargin<13 || isempty(n_theta)
    n_theta = 180*2;
end
if n_theta <= 10
    error('cotrans:InvalidArgument',...
        'Invalid n_theta.')
end

% M_threads check
if nargin<15 || isempty(M_threads)
    M_threads = 8;
end
if M_threads < 1
    warning('cotrans:InvalidArgument',...
        'M_threads should be greater than or equal to 1. Set to default value.')
    M_threads = 8;
end

% N_threads check
if nargin<16 || isempty(N_threads)
    N_threads = 2;
end
if N_threads < 1
    warning('cotrans:InvalidArgument',...
        'N_threads should be greater than or equal to 1. Set to default value.')
    N_threads = 2;
end

% shouldKeepContour check
shouldKeepContour = nargout>3;

%% Call MEX entry
[ Ls, K, Phi0, XorRc, YorPhic, Phif, Thetaf ] =...
    UBKLstarcxx(x0, y0, pa0, ...
    [Y(:), DOY(:), H(:), MN(:), S(:)]', ...
    ioptparmod, external, internal, ...
    ionoR, ds, dXorR, dYorPhi, n_phi, n_theta, ...
    shouldKeepContour, isCatesianGrid, M_threads, N_threads);

%% Post-processing

end
