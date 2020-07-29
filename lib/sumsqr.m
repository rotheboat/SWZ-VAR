function [s, n] = sumsqr( x )
% For the finite elements of an array x, returns the sum of squared 
% elements, and the number of elements that are finite. Ignores NaN values.
% Author: Roland Meeks
% Date: 5 June 2020

% cat all elements of x
y = x(:);
% Find inf values
b_isinf = isinf( y );
% Count the number of finite elements
n = numel(x) - sum( b_isinf );

% square all elements that are finite
z = y(~b_isinf).^2;
% sum of squares
s = nansum( z );

end