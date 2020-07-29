function [y] = rndu(r,c)
% PURPOSE:
% Creates a matrix of uniform (pseudo) random variables. 
% 
% FORMAT:
% y = rndu(r,c); 
% 
% INPUT:
% r    scalar, row dimension. 
% c    scalar, column dimension. 
% 
% OUTPUT:
% y    rx c matrix of uniform random variables between 0 and 1
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

y=rand(r,c);