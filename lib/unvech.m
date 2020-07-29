function Matrix = unvech(Vector, varargin)
% This function implements the unvech operator.
%
% INPUTS
%   Vector             [double]   a m*1 vector.
%   fill_zeros         [logical]  an optional true/false flag
%
% OUTPUTS
%   Matrix             [double]   a n*n symetric matrix, where n solves n*(n+1)/2=m.
%                                 optionally a n*n matrix with zeros in the
%                                 LT part. [RM]

% Copyright (C) 2010 Dynare Team
% Modified 2014 by Roland Meeks
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

% RM check for additional arguments
nargin = length( varargin );
if 0 < nargin
    fill_zeros = varargin{1};
else
    fill_zeros = 0;
end;

% RM error catching
if ~any( 1 == size( Vector ) )
    warning( 'unvech: Expecting a vector as input' );
end

m = length(Vector);
n = (sqrt(1+8*m)-1)/2;

if ~isint(n)
    warning( 'unvech: Vector must correspond to some square matrix' );
end;

b = 0;
Matrix = zeros(n,n); % RM changed from NaN
for col = 1:n
    idx = 1:col;
    Matrix(1:col,col) = Vector(b+idx);
    if ~fill_zeros
        Matrix(col,1:col) = transpose(Matrix(1:col,col));
    end;
    b = b+length(idx);
end

% if fill_zeros
%     Matrix( logical( tril( ones( n, n ), -1 ) ) ) = 0;
% end;