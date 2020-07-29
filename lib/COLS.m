function [y] = cols(x)

% PURPOSE:
% Returns the number of columns in a matrix. 
% 
% FORMAT:
% y = rows(x); 
% 
% INPUT:
% x     NxK matrix or sparse matrix. 
%
% OUTPUT: 
% y     scalar, number of rows in the specified matrix. 
%
% REMARKS: 
% If x is an empty matrix, rows(x) and cols(x) return 0. 
% 
% EXAMPLE: 
%      x = ones(3,5);
%      y = cols(x);
% 
% x =    1    1    1    1    1 
%        1    1    1    1    1 
%        1    1    1    1    1 
%   
% y = 5
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

y=size(x,2);