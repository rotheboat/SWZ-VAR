function [y] = trimr(x,t,b)
% PURPOSE:
% Trims rows from the top and/or bottom of a matrix. 
% 
% FORMAT:
% y = trimr( x, t, b); 
% 
% INPUT:
% x    NxK matrix from which rows are to be trimmed. 
% 
% t    scalar containing the number of rows which are to be removed from the top of x. 
% 
% b    scalar containing the number of rows which are to be removed from the bottom of x. 
% 
% OUTPUT: 
% y    RxK matrix where R=N-( t+ b), containing the rows left after the trim. 
% 
% REMARKS: 
% If either t or b is zero, then no rows will be trimmed from that end of the matrix. 
% 
% EXAMPLE:
% 
%      x = rndu(5,3);
%      y = trimr(x,2,1);
% 
%                                 x =    0.76042751    0.33841579    0.01844780 
%                                        0.05334503    0.38939785    0.65029973 
%                                        0.93077511    0.06961078    0.04207563 
%                                        0.53640701    0.06640062    0.07222560 
%                                        0.14084669    0.06033813    0.69449247 
%   
% 
%                                 y =    0.93077511    0.06961078    0.04207563 
%                                        0.53640701    0.06640062    0.07222560 
%   
% 
% See Also 
% submat, rotater, shiftr 
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

if nargin ~= 3
    error('Wrong number of arguments')
end

x([1:t size(x,1)-b+1:size(x,1)],:)=[];
y=x;