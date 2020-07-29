function [y] = ln(x)
% PURPOSE:
% Computes the natural log of all elements of x. 
% 
% FORMAT:
% y = ln( x); 
% 
% INPUT:
% x    NxK matrix or N-dimensional array. 
%
% OUTPUT: 
% y    NxK matrix or N-dimensional array containing the natural log values of the elements of x. 
% 
% REMARKS:
% ln is defined for x ? 0. 
% If x is negative, complex results are returned. 
% You can turn the generation of complex numbers for negative inputs on or
% off in the GAUSS configuration file, and with the sysstate function, case
% 8. If you turn it off, ln will generate an error for negative inputs.If x is already 
% complex, the complex number state doesn’t matter; ln will compute a complex result. 
% x can be any expression that returns a matrix. 
% 
% EXAMPLE: 
% 
%      y = ln(16);
% 
%      y = 2.7725887
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

y = log(x);