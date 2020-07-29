function y = mvtpdfrm(X, C, df)
%MVTPDF Multivariate t probability density function (pdf) without normalization.
%   Y = MVTPDF(X,C,DF) returns the probability density of the multivariate t
%   distribution with correlation parameters C and degrees of freedom DF,
%   evaluated at each row of X.  Rows of the N-by-D matrix X correspond to
%   observations or points, and columns correspond to variables or
%   coordinates.  Y is an N-by-1 vector.
%
%   C is a symmetric, positive definite, D-by-D matrix, typically a
%   correlation matrix.  If its diagonal elements are not 1, MVTPDF scales
%   C to correlation form.  DF is a scalar, or a vector with N elements.
%
%   MVTPDFRM omits the normalization step, but is otherwise identical
% 
%   Example:
%
%      C = [1 .4; .4 1]; df = 2;
%      [X1,X2] = meshgrid(linspace(-2,2,25)', linspace(-2,2,25)');
%      X = [X1(:) X2(:)];
%      p = mvtpdf(X, C, df);
%      surf(X1,X2,reshape(p,25,25));
%
%   See also MVNPDF, MVTCDF, MVTRND, TPDF.

%   Copyright 2005-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2006/11/11 22:55:36 $

if nargin<3
    error('stats:mvtpdf:TooFewInputs','Requires three inputs.');
elseif ndims(X)~=2
    error('stats:mvtpdf:InvalidData','X must be a matrix.');
end

% Get size of data.  Column vectors provisionally interpreted as multiple scalar data.
[n,d] = size(X);
if d<1
    error('stats:mvtpdf:TooFewDimensions','X must have at least one column.');
end

% Special case: try to interpret X as a row vector if it was a column.
if isvector(X) && (size(C,1) == n)
    X = X';
    [n,d] = size(X);
end

sz = size(C);
if sz(1) ~= sz(2)
    error('stats:mvtpdf:BadCorrelation',...
          'C must be a square matrix.');
elseif ~isequal(sz, [d d])
    error('stats:mvtpdf:InputSizeMismatch',...
          'C must be a square matrix with size equal to the number of columns in X.');
end

% *************************************************************************
% Standardize C to correlation if necessary.  This does NOT standardize X.
s = sqrt(diag(C));
if (any(s~=1))
    C = C ./ (s * s');
end
% *************************************************************************

% Make sure C is a valid covariance matrix
[R,err] = cholcov(C,0); % R is the Choleski factor of C
if err ~= 0
    error('stats:mvtpdf:BadCorrelation',...
        'C must be symmetric and positive definite.');
end

if ~(isscalar(df) || (isvector(df) && length(df) == n))
    error('stats:mvtpdf:InputSizeMismatch', ...
          'DF must be a scalar or a vector with one element for each row in X.');
elseif any(df <= 0)
    error('stats:mvtpdf:InvalidDF','DF must be positive.');
end
df = df(:);

% Create array of standardized data, and compute log(sqrt(det(Sigma)))
Z = X / R; % mrdivide X*R^{-1}
logSqrtDetC = sum(log(diag(R)));

logNumer = -((df+d)/2) .* log(1+sum(Z.^2, 2)./df);
logDenom = logSqrtDetC + (d/2)*log(df*pi);
y = exp(gammaln((df+d)/2) - gammaln(df/2) + logNumer - logDenom);
