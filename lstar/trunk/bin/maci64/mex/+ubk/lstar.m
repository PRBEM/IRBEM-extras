function [ Ls K Phi0 XorRc YorPhic Phif Thetaf ] = ...
    lstar( x0, y0, pa0, datenum_, ioptparmod, external, internal, ...
    varargin)
%UBK.LSTAR_CART Lstar in Cartesian Grid
%   [ Ls K Phi0 XorRc YorPhic Phif Thetaf ] = ...
%      lstar( x0, y0, pa0, datenum_, ioptparmod, external, internal, ...
%      param1, value1, ...)
%   Traces particle's drift contour and evaluates generalized L value, L*
%   (Roederer 1970), using built-in Tsyganenko field models.
%   SM coordinate system is default (and the only choice).
%
%   References: Tsyganenko 1987, 1989, 1995, 2002a, 2002b;
%   Tsyganenko and Sitnov 2005; Tsyganenko and Sitnov 2007;
%   Sitnov et al. 2008; Sheldon and Gaffey 1993; Min et al. 2013a,
%   2013b;
%   Please refer to http://geo.phys.spbu.ru/~tsyganenko/modeling.html,
%   http://geomag_field.jhuapl.edu/model/ and references therein.
%
%   REQUIRED INPUTS
%   * [x0 y0]: Matrices of SM cartesian/cylindrical coordinates MAPPED ON
%   THE X-Y PLANE ALONG THE FIELD LINE. RE. M by N where N ==
%   length(datenum_).
%   * pa0: Initial pitch angle in radian AT THE MAGNETIC EQUATOR. M by N.
%   * datenum_: N-element time vector returned by datenum.
%   * ioptparmod: External field parameters.
%   IF external=='T89', THEN: ioptparmod = IOPT (length N), where
%       T89 input equivalent to Kp.
%       IOPT= 1       2        3        4        5        6        7
%       KP  = 0,0+  1-,1,1+  2-,2,2+  3-,3,3+  4-,4,4+  5-,5,5+  >=6-
%   ELSE IF external=='TS07', THEN: ioptparmod = PARMOD, where
%       PARMOD is a matrix of size [102, N] whose columns are
%       [Pdyn (nPa); TS07_COEF (101, fitting coefficients)].
%   ELSE IF external~='T89' && external~='NONE' && external~='TS07', THEN:
%   ioptparmod = PARMOD, where
%       PARMOD is a matrix of size [10, N] whose columns are
%       [PDYN (nPa), DST (nT), BYIMF, BZIMF (nT), W1 W2 W3 W4 W5 W6].
%   * external: External field model. 'NONE' for no external model, 'T89' for
%   T89 model, 'T96' for T96 model, 'T02' for T02 model, 'TS05' for TS05
%   model, or 'TS07' for TS07D model (case-insensitive). Default is 'NONE'.
%   * internal: Internal field model. 'DIP' for dipole field or 'IGRF' for
%   IGRF model (case-insensitive). Default is 'IGRF'.
%
%   OPTIONAL INPUTS
%   * IONOR: Ionospheric boundary (radial distance to the ionosphere) in RE.
%   Default is 1.015 RE (~ 100 km above Earth).
%   * DS: Field line integration step size. Default is 0.1 RE.
%   * DXorDR: X/R grid resolution in RE. Default is 0.2 RE.
%   * DYorDPHI: Y/Phi grid resolution in RE/radian. Y default is 0.2 RE and
%   Phi default is 2pi / n_phi.
%   * N_PHI: The number of azimuthal grids for magnetic flux evaluation.
%   Default is 350/5.
%   * N_THETA: The number of polar grids for magnetic flux evaluation.
%   Default is 180*2.
%   * isCARTESIANGRID: Boolean to indicate grid scheme to use.
%   If true, cartesian grid is used, else cylindrical grid is
%   used. Input/output coordinate system depends on this setting except for
%   the foot point coordinates which is always in SM spherical coordinate
%   system. Default is false.
%   * M_THREADS: The number of threads to loop over the first dimension
%   (positive integer). Default is 8.
%   * N_THREADS: The number of threads to loop over the second dimension
%   (positive integer). Default is 2.
%   * USE_INTERP_FIELD_WITH_POLY_ORDER: Either 1 (linear) or 2 (2nd order).
%   If set, caches the field vectors at discrete points and interpolates the
%   field vector at the given location from the cached vectors.
%   The interpolation may increase the calculation speed if the underlying
%   field model is expensive, but the result includes the numerical error
%   from the interpolation and the interpolated vector violates divergence
%   free contidion. Default is unset (equivalently pass 0).
%
%   OUTPUTS (M by N otherwise specified)
%   * Ls: Matrix of L*.
%   * K: Matrix of modified 2nd invariant in nT^.5 RE.
%   * Phi0: N-element vector of Earth's magnetic moment times 2*pi.
%   Magnetic flux defined by the drift contour is Phi0/Ls (positive).
%   nT RE^2.
%   * [XorRc YorPhic]: Cell matrices of the drift contour coordinates.
%   RE and radian in SM coordinate system.
%   * [Phif Thetaf]: Cell matrices of the foot points of the drift
%   contour. RE and radian in SM spherical coordinate system.
%
%   See also ubk.fieldline

%
% $Author$
% $LastChangedDate$
% $Revision$
% $Id$
%

%% Argument check
error(nargchk(7,nargin,nargin))

%% Required
MxN = size(x0);

% Size check
if ~isequal(size(x0), size(y0), size(pa0), [MxN(1) length(datenum_)])
    error('ubkL8:InvalidArgument',...
        'size of x, y and pa is different.')
end

% Date check
if ~isnumeric(datenum_)
    error('ubkL8:InvalidArgument',...
        'Non-numeric datenum_ input.')
end
if length(datenum_)~=numel(datenum_)
    error('ubkL8:InvalidArgument',...
        'Invalid datenum_ dimension.')
end
[Y, ~, ~, H, MN, S] = datevec(datenum_);
tmp = datenum([Y(:), zeros(numel(Y),1), ones(numel(Y),1)]);
DOY = floor( datenum_ - reshape(tmp, size(datenum_)) ) + 1;

% External field
if ~ischar(external)
    error('ubkL8:InvalidArgument',...
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
    case 'ts07'
        external = 5;
    otherwise
        error('ubkL8:InvalidArgument',...
        'External model %s is unknown.', external)
end

if external==1
    if length(ioptparmod)~=MxN(2)
        error('ubkL8:InvalidArgument',...
            'Incompatible IOPT length.')
    end
    ioptparmod = ioptparmod(:)';
elseif external==5
    if ~isequal(size(ioptparmod),[102 MxN(2)])
        error('ubkL8:InvalidArgument',...
            'Incompatible PARMOD length.')
    end
elseif (external~=1&&external~=0)
    if ~isequal(size(ioptparmod),[10 MxN(2)])
        error('ubkL8:InvalidArgument',...
            'Incompatible PARMOD size.')
    end
end

% Internal field
if ~ischar(internal)
    error('ubkL8:InvalidArgument',...
        'internal is not character array.')
end
switch lower(internal)
    case 'dip'
        internal = 0;
    case 'igrf'
        internal = 1;
    otherwise
        error('ubkL8:InvalidArgument',...
        'Internal model %s is unknown.', internal)
end

%% Options
opts = ubk.optset(varargin{:});

% Ionospheric boundary check
if isempty(opts.IONOR)
    opts.IONOR= 1.015;
end
if opts.IONOR < 1
    error('ubkL8:InvalidArgument',...
        'IONOR should be greater than or equal to 1.')
end

% Integration step size check
if isempty(opts.DS)
    opts.DS = .1;
end
if opts.DS <= 0
    error('ubkL8:InvalidArgument',...
        'DS should be greater than 0.')
end

% dx/dr check
if isempty(opts.DXORDR)
    opts.DXORDR = .2;
end
if opts.DXORDR <= 0
    error('ubkL8:InvalidArgument',...
        'Invalid DXorDR.')
end

% n_phi check
if isempty(opts.N_PHI)
    opts.N_PHI = 360/5;
end
if opts.N_PHI <= 10
    error('ubkL8:InvalidArgument',...
        'Invalid N_PHI.')
end

% isCaRtesianGrid check
if isempty(opts.ISCARTESIANGRID)
    opts.ISCARTESIANGRID = false;
end

% dy/dphi check
if isempty(opts.DYORDPHI)
    if opts.ISCARTESIANGRID
        opts.DYORDPHI = .2;
    else
        opts.DYORDPHI = 2*pi / opts.N_PHI;
    end
end
if opts.DYORDPHI <= 0
    error('ubkL8:InvalidArgument',...
        'Invalid DYorDPHI.')
end

% n_theta check
if isempty(opts.N_THETA)
    opts.N_THETA = 180*2;
end
if opts.N_THETA <= 10
    error('ubkL8:InvalidArgument',...
        'Invalid N_THETA.')
end

% M_threads check
if isempty(opts.M_THREADS)
    opts.M_THREADS = 8;
end
if opts.M_THREADS < 1
    warning('ubkL8:InvalidArgument',...
        'M_THREADS should be greater than or equal to 1. Set to default value.')
    opts.M_THREADS = 8;
end

% N_threads check
if isempty(opts.N_THREADS)
    opts.N_THREADS = 2;
end
if opts.N_THREADS < 1
    warning('ubkL8:InvalidArgument',...
        'N_THREADS should be greater than or equal to 1. Set to default value.')
    opts.N_THREADS = 2;
end

% USE_INTERP_FIELD_WITH_POLY_ORDER check
if isempty(opts.USE_INTERP_FIELD_WITH_POLY_ORDER)
    opts.USE_INTERP_FIELD_WITH_POLY_ORDER = 0;
end
if opts.USE_INTERP_FIELD_WITH_POLY_ORDER<0 || ...
        opts.USE_INTERP_FIELD_WITH_POLY_ORDER>2
    error('ubkL8:InvalidArgument',...
        'USE_INTERP_FIELD_WITH_POLY_ORDER option should be between 0 and 2.')
end

% shouldKeepContour check
shouldKeepContour = nargout>3;

%% Call MEX entry
[ Ls, K, Phi0, XorRc, YorPhic, Phif, Thetaf ] =...
    UBKLstarcxx(double(x0), double(y0), double(pa0), ...
    [Y(:), DOY(:), H(:), MN(:), S(:)]', ...
    double(ioptparmod), external, internal, ...
    opts.IONOR, opts.DS, opts.DXORDR, opts.DYORDPHI, ...
    opts.N_PHI, opts.N_THETA, ...
    shouldKeepContour, opts.ISCARTESIANGRID, ...
    opts.M_THREADS, opts.N_THREADS, opts.USE_INTERP_FIELD_WITH_POLY_ORDER);

%% Post-processing

end
