function [alpha0,beta0,phib,result_code] = rfl_vectors_to_euler_angles(B,C,S0,S1)
% [alpha0,beta0,phib,result_code] = rfl_vectors_to_euler_angles(B,C,S0,S1)
% Compute pitch angle (alpha0) and gyrophase angle (beta0) of boresight
% and longitude of B (phib), the Euler angles of the Euler
% rotations for converting between instrument coordinates 
% and magnetic coordinates. All four inputs are 3-vectors in some common 
% Cartesian coordinate system (e.g., Cartesian GEI). B is parallel to the 
% magnetic field, C defines the B-C plane, in which beta=0. 
% S0 points into the instrument, parallel to normally incident particles. 
% S1 defines the S0-S1 plane in
% which phi=0. Vectors need not have unit length, and neither B-C nor S0-S1
% need to be orthogonal
% expect input vectors to be Nx3
% output angles will be Nx1, in degrees

result_code = 1; % SUCCESS

if (numel(B)==3) && (numel(C)==3) && (numel(S0)==3) && (numel(S1)==3),
    B = B(:)';
    C = C(:)';
    S0 = S0(:)';
    S1 = S1(:)';
end

[NB,B3] = size(B);
[NC,C3] = size(C);
[N0,S03] = size(S0);
[N1,S13] = size(S1);

if (NC ~= NB) || (N0 ~= NB) || (N1 ~= NB) || (B3 ~= 3) || (C3 ~= 3) || (S03 ~= 3) || (S13 ~= 3),
    if nargout >= 2,
        alpha0 = [];
        beta0 = [];
        phib = [];        
        result_code = -1; % incorrect argument size
        return;
    else
        error('Incorrect size of inputs');
    end
end

hat = @(x)x./repmat(sqrt(sum(x.^2,2)),1,3);

bhat = hat(B);
dhat = hat(cross(B,C));
chat = cross(dhat,bhat);
S0hat = hat(S0);
S2hat = hat(cross(S0,S1));
S1hat = cross(S2hat,S0hat);
cosa0 = max(-1,min(1,dot(bhat,S0hat,2))); % bound cos in [-1 1]
alpha0 = acosd(cosa0);

phib = nan(NB,1); % default phib is NaN
i0 = alpha0==0;
if any(i0(:)),
    phib(i0) = 180/pi*atan2(dot(S1hat(i0,:),dhat(i0,:),2),-dot(S2hat(i0,:),dhat(i0,:),2));
end
beta0 = zeros(NB,1); % beta0=0 when alpha0=0
beta0(~isfinite(alpha0)) = nan; % beta0=NaN when alhpa0=NaN

i0 = alpha0>0;
if any(i0(:)),
    beta0(i0) = 180/pi*atan2(dot(S0hat(i0,:),dhat(i0,:),2),-dot(S0hat(i0,:),chat(i0,:),2));
    phib(i0) = 180/pi*atan2(dot(S2hat(i0,:),bhat(i0,:),2),-dot(S1hat(i0,:),bhat(i0,:),2));
end
