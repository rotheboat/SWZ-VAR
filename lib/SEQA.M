function y = seqa(start,inc,n)

% PURPOSE: 
% seqa creates an additive sequence. 
% 
% FORMAT: 
% y = seqa( start, inc, n); 
% 
% INPUT: 
% start     scalar specifying the first element. 
% inc       scalar specifying increment. 
% n         scalar specifying the number of elements in the sequence. 
%
% OUTPUT: 
% y         nx1 vector containing the specified sequence. 
%
% REMARKS:
% For seqa, y will contain a first element equal to start, the second
% equal to start+ inc, and the last equal to start+ inc*( n-1). 
% 
% For instance, 
%          seqa(1,1,10); 
% will create a column vector containing the numbers 1, 2, ...10. 
%
% EXAMPLE: 
% 
%      a = seqa(2,2,10)';
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

finish = start + (n-1)*inc;
y = [start:inc:finish]';