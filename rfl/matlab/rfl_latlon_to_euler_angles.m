function varargout = rfl_latlon_to_euler_angles(Blat,Blon,Clat,Clon,S0lat,S0lon,S1lat,S1lon)
% [alpha0,beta0,phib,result_code] = rfl_latlon_to_euler_angles(Blat,Blon,Clat,Clon,S0lat,S0lon,S1lat,S1lon)
% Same as vectors_to_euler_angles, except B,C, S0, S1 specified by angle pairs.
% All three angle pairs are degrees latitude and longitude in some common spherical
% coordinate system (e.g., spherical GEO).
% all inputs must be the same size; alpha0,beta0,phib will have same size
% as inputs, and will be in degrees

varargout = cell(1,nargout);

N = numel(Blat);
Nd = ndims(Blat);
sz = size(Blat);
if (numel(Blon) ~= N) || (numel(Clat) ~= N) || (numel(Clon) ~= N) || ...
        (numel(S0lat) ~= N) || (numel(S0lon) ~= N) || (numel(S1lat) ~= N) || (numel(S1lon) ~= N) || ...
        (ndims(Blon) ~= Nd) || (ndims(Clat) ~= Nd) || (ndims(Clon) ~= Nd) || ...
        (ndims(S0lat) ~= Nd) || (ndims(S0lon) ~= Nd) || (ndims(S1lat) ~= Nd) || (ndims(S1lon) ~= Nd) || ...
        any(size(Blon) ~= sz) || any(size(Clat) ~= sz) || any(size(Clon) ~= sz) || ...
        any(size(S0lat) ~= sz) || any(size(S0lon) ~= sz) || any(size(S1lat) ~= sz) || any(size(S1lon) ~= sz),
    if nargout >= 2,
        alpha0 = [];
        beta0 = [];
        phib = [];
        result_code = -1; % incorrect argument sz
        varargout = {alpha0,beta0,phib,result_code};
        return;
    else
        error('Incorrect size of inputs');
    end
end

% convert to column vectors
Blat = Blat(:);
Blon = Blon(:);
Clat = Clat(:);
Clon = Clon(:);
S0lon = S0lon(:);
S1lon = S1lon(:);

% convert to unit vectors
B = [cosd(Blat).*cos(Blon),cosd(Blat).*sin(Blon),sin(Blat)];
C = [cosd(Clat).*cos(Clon),cosd(Clat).*sin(Clon),sin(Clat)];
S0 = [cosd(S0lat).*cos(S0lon),cosd(S0lat).*sin(S0lon),sin(S0lat)];
S1 = [cosd(S1lat).*cos(S1lon),cosd(S1lat).*sin(S1lon),sin(S1lat)];

% use unit vector routine
[varargout{:}] = rfl_vectors_to_euler_angles(B,C,S0,S1);

% reshape outputs to common shape (sz)
for i = 1:length(varargout),
    if numel(varargout{i})==N,
        varargout{i} = reshape(varargout{i},sz);
    end
end
