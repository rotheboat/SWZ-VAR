function [y] = stof(x)

% PURPOSE:
% Converts a string to floating point. 
% 
% FORMAT: 
% y = stof( x); 
% 
% INPUT:
% x     string or NxK matrix containing character elements to be converted. 
% 
% OUTPUT: 
% y     matrix, the floating point equivalents of the ASCII numbers in x. 
% 
% REMARKS: 
% If x is a string containing “1 2 3”, then stof will return a 3x1 matrix
% containing the numbers 1, 2 and 3. If x is a null string, stof will
% return a 0. This uses the same input conversion routine as loadm and let.
% It will convert character elements and missing values. stof also converts
% complex numbers in the same manner as let.
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

y = str2num(x); %#ok<ST2NM>