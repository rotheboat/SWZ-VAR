function [Xlag] = mlagrm(X, nlag)
%MLAG2 Summary of this function goes here
%   The dimension of the matrix being passed is T*-h
[T,dimvar]=size(X);

Xlag=zeros(T-nlag+1,dimvar*nlag);

for lag_i = 1:nlag,
    Xlag(:, dimvar*(lag_i-1)+1:dimvar*lag_i) = ...
						X(nlag+1-lag_i:T+1-lag_i, 1:dimvar);
end