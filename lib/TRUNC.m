function [y] = trunc(x)
% PURPOSE: 
% Converts numbers to integers by truncating the fractional portion. 
% 
% FORMAT:
% y = trunc( x); 
% 
% INPUT:
% x   NxK matrix or N-dimensional array. 
% 
% OUTPUT:
% y   NxK matrix or N-dimensional array containing the truncated elements of x. 
% 
% EXAMPLE:
%  
%        x = 100*rndn(2,2);
%        y = trunc(x);
% 
%                                    x =    77.68     ?14.10 
%                                            4.73    ?158.88 
%   
% 
%                                    y =    77.00     ?14.00 
%                                            4.00    ?158.00 
%   
% See Also 
% ceil, floor, round 
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

if nargin ~= 1
    error('Wrong number of input arguments')
end

y = fix(x);