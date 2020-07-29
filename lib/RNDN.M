function [y] = rndn(r,c)
% PURPOSE:
% Creates a matrix of standard Normal (pseudo) random numbers. 
% 
% FORMAT:
% y = rndn( r, c); 
% 
% INPUT:
% r    scalar, row dimension. 
% 
% c    scalar, column dimension. 
% 
% OUTPUT:
% y    rx c matrix of Normal random numbers having a mean of 0 and standard deviation of 1
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

y=randn(r,c);