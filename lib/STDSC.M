function [y] = stdsc(x)

% PURPOSE: 
% Computes the standard deviation of the elements in each column of a matrix. 
% 
% FORMAT:
% y = stdsc(x); 
% 
% INPUT:
% x    NxK matrix. 
% 
% OUTPUT:
% y    Kx1 vector, the standard deviation of each column of x. 
% 
% REMARKS: 
% This function essentially computes: 
%
%            sqrt(1/(N)*sumc(( x-meanc( x)’)2 )) 
%
% Thus, the divisor is N rather than N-1, where N is the number of elements
% being summed. See stdc for the alternate definition. 
% 
% EXAMPLE: 
% 
%      y = rndn(8100,1);
%      std = stdsc(y);
% 
%                                        std = 1.0151475 
% 
% In this example, 8100 standard Normal random variables are generated, and
% their standard deviation is computed. 
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

y = std(x,1,1);