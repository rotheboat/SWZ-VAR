function [y]=counts(x,v)
% PURPOSE:
% Counts the numbers of elements of a vector that fall into specified ranges. 
% 
% FORMAT:
% c = counts( x, v); 
% 
% INPUT:
% x   Nx1 vector containing the numbers to be counted. 
% 
% v   Px1 vector containing breakpoints specifying the ranges within which counts are to be
%     made. The vector v MUST be sorted in ascending order. 
% 
% OUTPUT: 
% 
% c 
% Px1 vector, the counts of the elements of x that fall into the regions: 
% 
%   x ?v[1],   
% v[1] < x ?v[2],   
%        
% v[p - 1] < x ?v[p]   
% 
% 
% 
% 
% REMARKS:
% If the maximum value of x is greater than the last element (the maximum value) of v, the 
% sum of the elements of the result, c, will be less than N, the total number of elements in x. 
% 
%            1 
%            2 
%            3 
%            4                4
% If    x =  5     and   v =  5 
%            6                8
%            7 
%            8 
%            9 
%   
% 
%            4 
% then  c =  1 
%            3 
% 
% The first category can be a missing value if you need to count missings directly. 
% Also + ? or - ? are allowed as breakpoints. The missing value must be the first breakpoint
% if it is included as a breakpoint and infinities must be in the proper location depending
% on their sign. - ? must be in the [2,1] element of the breakpoint vector if there is 
% a missing value as a category as well, otherwise it has to be in the [1,1] element. 
% If + ? is included, it must be the last element of the breakpoint vector. 
% 
% EXAMPLE: 
% 
%      x = [ 1, 3, 2,
%            4, 1, 3 ];
%      v = [ 0, 1, 2, 3, 4 ];
%      c = counts(x,v);
% 
% c =    0.0000000 
%        2.0000000 
%        1.0000000 
%        2.0000000 
%        1.0000000 


if nargin ~= 2
    error('Wrong number of inputs in counts')
end

% Check dimensions of input matrix
r=size(v,1);
c=size(v,2);
if r>c && c==1
    n=r;
    v=v';
elseif c>r && r==1
    n=c;
else
    error('Input object in counts determining thresholds is not a vector')
end

if isreal(v) == 0 || isreal(x)==0
    error('Procedure counts cannot be implemented for complex inputs')
end

y=[];
if isnan(v(1))==1
    y(1)=size(find(isnan(x)),2);
elseif v(1) == -Inf
    y(1)=size(find(isinf(x)),2);
else
    y(1)=size(find(x<v(1)),2);
end

for i=2:n
    y(i,1)=size(find(x<=v(i) & x>v(i-1)),2);
end