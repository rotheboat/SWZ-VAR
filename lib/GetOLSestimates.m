% Roland Meeks (2010-11)
%
% Given matrices Y and X, computes coefficients A in Y = AX + e, and
% SIGMA_OLS the covariance matrix of e.
function s = GetOLSestimates( Y, X )

s = struct( 'A_OLS', [], 'SSE', [], 'SIGMA_OLS', [], 'XpX', [], ...
    'XpY', [], 'Y', [], 'X', [], 'Y_fitted', [] );
[T K] = size(X);

s.XpX = (X'*X);

s.XpY = (X'*Y);

s.A_OLS = s.XpX\s.XpY; % This is the matrix of regression coefficients

s.SSE = (Y - X*s.A_OLS)'*(Y - X*s.A_OLS);

s.SIGMA_OLS = (s.SSE)./(T-K);

s.V_OLS = kron( s.SIGMA_OLS, inv( s.XpX ) );

% polite to also store the data in case others need it
s.Y = Y;
s.X = X;
s.Y_fitted = X*(s.A_OLS);
