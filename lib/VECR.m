function [y]=vecr(x)
% PURPOSE:
% Creates a column vector by appending the rows of a matrix to each other. 
% 
% FORMAT:
% yr = vecr( x); 
% 
% INPUT:
% x     NxK matrix. 
% 
% OUTPUT 
% yr   (N*K)x1 vector, the rows of x appended to each other and the result transposed. 
% 
% REMARKS:
% 
% EXAMPLE: 
% 
%      x = [ 1 2,
%            3 4 ];
%      yr = vecr(x);
% 
% x =    1.000000    2.000000 
%        3.000000    4.000000 
%
% yr =    1.000000 
%         2.000000 
%         3.000000 
%         4.000000 
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

y=reshape(x',size(x,1)*size(x,2),1);