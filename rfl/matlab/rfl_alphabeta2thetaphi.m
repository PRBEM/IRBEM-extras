function [theta,phi,result_code] = rfl_alphabeta2thetaphi(alpha,beta,alpha0,beta0,phib)
% [theta,phi,result_code] = rfl_alphabeta2thetaphi(alpha,beta,alpha0,beta0,phib)
% Convert pitch angle and gyrophase to instrument angles. All angles
% degrees. When theta=0,180 phi=0.

% columns are coefficients of c,d,b terms of inst basis vectors
% rows are coefficients of s1,s2,s0 terms of mag basis vectors
result_code = 1;

if (numel(alpha0) ~= 1) || (numel(beta0) ~= 1) || (numel(phib) ~= 1) || ...
        (numel(alpha) ~= numel(beta)) || (ndims(alpha) ~= ndims(beta)) || any(size(alpha)~=size(beta)),
    if nargout >= 2,
        theta = [];
        phi = [];
        result_code = -1; % incorrect argument size
        return;
    else
        error('Incorrect size of inputs');
    end
end

sz = size(alpha);

alpha = alpha(:);
beta = beta(:);

R = [
    - sind(beta0)*sind(phib) - cosd(alpha0)*cosd(beta0)*cosd(phib),   cosd(phib)*sind(beta0) - cosd(alpha0)*cosd(beta0)*sind(phib), cosd(beta0)*sind(alpha0)
    cosd(beta0)*sind(phib) - cosd(alpha0)*cosd(phib)*sind(beta0), - cosd(beta0)*cosd(phib) - cosd(alpha0)*sind(beta0)*sind(phib), sind(alpha0)*sind(beta0)
    cosd(phib)*sind(alpha0),                                     sind(alpha0)*sind(phib),            cosd(alpha0)
    ];

sa = sind(alpha);
ca = cosd(alpha);
sb = sind(beta);
cb = cosd(beta);

theta = acosd(R(3,3)*ca+R(1,3)*sa.*cb + R(2,3)*sa.*sb);
phi = 180/pi*atan2(R(3,2)*ca+R(1,2)*sa.*cb+R(2,2)*sa.*sb,R(3,1)*ca+R(1,1)*sa.*cb+R(2,1)*sa.*sb);

% when theta = 0,180, phi=0
phi((theta==0) | (theta==180)) = 0;

theta = reshape(theta,sz);
phi = reshape(phi,sz);
