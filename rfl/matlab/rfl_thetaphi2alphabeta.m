function varargout = rfl_thetaphi2alphabeta(theta,phi,alpha0,beta0,phib)
% [alpha,beta,result_code] = thetaphi2alphabeta(theta,phi,alpha0,beta0,phib);
% Convert instrument angles to pitch angle and gyrophase. All angles
% degrees. When alpha=0,180 beta=0

varargout = cell(1,nargout);
[varargout{:}] = rfl_alphabeta2thetaphi(theta,phi,alpha0,phib,beta0); % note phib <-> beta0

% the alpha,beta -> theta,phi transform is identical under
% alpha <-> theta
% beta <-> phi
% beta0 <-> phib (makes R <-> R')
