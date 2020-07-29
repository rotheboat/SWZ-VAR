function [y] = rev(x)

% PURPOSE: 
% Reverses the order of the rows in a matrix. 
% 
% FORMAT:
% y = rev( x); 
% 
% INPUT:
% x   NxK matrix. 
% 
% OUTPUT:
% y   NxK matrix containing the reversed rows of x. 
% 
% REMARKS: 
% The first row of y will be where the last row of x was and the last row will be where the
% first was and so on. This can be used to put a sorted matrix in descending order. 
% 
% EXAMPLE:
% x = round(rndn(5,3)*10);
% y = rev(x);
% 
% x =    10     7     8 
%         7     4    ?9 
%       ?11     0    ?3 
%         3    18     0 
%         9    ?1    20 
%   
% y =    9    ?1    20 
%        3    18     0 
%      ?11     0    ?3 
%        7     4    ?9 
%       10     7     8 
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

y = zeros(size(x,1),size(x,2));
for i=0:size(x,1)-1
    y(i+1,:)=x(size(x,1)-i,:);
end