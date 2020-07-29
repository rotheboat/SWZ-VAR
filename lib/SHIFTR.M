function [y] = shiftr(x,s,f)

% PURPOSE: 
% Shifts the rows of a matrix. 
% 
% FORMAT:% 
% y = shiftr( x, s, f); 
% 
% INPUT: 
% x     NxK matrix to be shifted. 
% s     scalar or Nx1 vector specifying the amount of shift. 
% f     scalar or Nx1 vector specifying the value to fill in. 
% 
% OUTPUT: 
% y     NxK shifted matrix. 
% 
% REMARKS: 
% The shift is performed within each row of the matrix, horizontally. If the shift value is positive,
% the elements in the row will be moved to the right. A negative shift value causes the elements to be
% moved to the left. The elements that are pushed off the end of the row are lost, and the fill value 
% will be used for the new elements on the other end. 
% 
% EXAMPLE:
% 
%      y = shiftr(x,s,f);
%                                                
%                 1    2             1             99 
%      If   x =             and s =     and  f =       
%                 3    4            -1             999 
%                                                
%  
%      Then  y =    99    1 
%                   4    999 
%   
% 
%                 1    2    3              0
%      If   x =   4    5    6   and  s =   1   and f=0
%                 7    8    9              2
%   
% 
% 
%      Then  y =    1    2    3 
%                   0    4    5 
%                   0    0    7 
%
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

if nargin == 2
    f=0;
elseif nargin == 1
    error('Wrong number of inputs')
elseif nargin > 3
    error('Wrong number of inputs')
end

if ndims(x)>2
    error('Not implemented for arrays')
end

rowsX = size(x,1);

if isscalar(s)
    s=s*ones(rowsX,1);
end
if isscalar(f)
    f=f*ones(rowsX,1);
end

y=zeros(size(x,1),size(x,2));

for i=1:rowsX
    B = circshift(x(i,:)',s(i));
    if s(i) < 0
        B(size(x,2) + s(i) + 1:size(x,2),1) = f(i);
    elseif s(i) > 0
        B(1:s(i)) = f(i);
    end
    y(i,:) = B';
end
