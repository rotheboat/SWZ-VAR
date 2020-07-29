function [y] = stdc(x)

% PURPOSE:
% Computes the standard deviation of the elements in each column of a matrix. 
% 
% FORMAT:
% y = stdc(x); 
% 
% INPUT: 
% x NxK matrix. 
% 
% OUTPUT: 
% y Kx1 vector, the standard deviation of each column of x. 
% 
% REMARKS: 
% This function essentially computes:
%
%             sqrt(1/(N-1)*sumc(( x-meanc( x)')2 )) 
%
% Thus, the divisor is N-1 rather than N, where N is the number of
% elements being summed. To convert to the alternate definition, multiply 
% by sqrt((N-1)/N).
%
% EXAMPLE:
% 
%      y = rndn(8100,1);
%      std = stdc(y);
% 
% std = 1.008377
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

y = std(x,0,1);