function response = impulse(By,smat,nstep)
% function response=impulse(By,smat,nstep), C. Sims' code.
% smat is a square matrix of initial shock vectors.  To produce "orthogonalized
% impulse responses" it should have the property that smat'*smat=sigma, where sigma
% is the Var(u(t)) matrix and u(t) is the residual vector.  One way to get such a smat
% is to set smat=chol(sigma).  To get the smat corresponding to a different ordering,
% use smat=chol(P*Sigma*P')*P, where P is a permutation matrix.
% By is a neq x nvar x nlags matrix.  neq=nvar, of course, but the first index runs over 
% equations. In response, the first index runs over variables, the second over 
% shocks (in effect, equations).

% Reduced form AR matrix is By:
%   Equations are in rows
%   Variables are in columns
%   Lags are in dimension 3
[neq,nvar,nlag]=size(By);

% Response matrix
%   Variables are in rows
%   Equations (shocks) are in columns
%   Lags are in dimension 3
response=zeros(nvar,neq,nstep);

% note that the 'chol' command produces an UT matrix s.t. if C = chol(A)
% then A = C'*C. Need a LT matrix, "last innovation untransformed(?)", so
% we take the transpose here to get that. As we usually think of impulse
% vectors as columns of C' (i.e. the lower triangular choleski factor),
% make sure that if a single impulse vector is supplied, it will be a
% column vector following this transpose operation!
response(:,:,1)=smat'; 
for it=2:nstep
   for ilag=1:min(nlag,it-1)
      response(:,:,it)=response(:,:,it)+By(:,:,ilag)*response(:,:,it-ilag);
   end
end
