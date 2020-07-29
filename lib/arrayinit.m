function [y] = arrayinit(o,v)

% PURPOSE: 
% Creates an N-dimensional array with a specified fill value. 
% 
% FORMAT: 
% y = arrayinit( o, v); 
% 
% INPUT: 
% o        Nx1 vector of orders, the sizes of the dimensions of the array. 
% v        scalar, value to initialize. If v is complex the result will be complex. 
%
% OUTPUT: 
% y        N-dimensional array with each element equal to the value of v. 
%
% EXAMPLE:
%
%      orders = { 2,3,4 };
%      y = arrayinit(orders, 0);
% 
%      y will be a 2x3x4 array of zeros. 
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

y = v*ones(o);