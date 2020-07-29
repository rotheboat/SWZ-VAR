function [c] = convgauss(b,x,f,l)
% NOTE: THIS FUNCTION IS EXACTLY THE FUNCTION 'conv' IN GAUSS. The name 'convgauss' is used to avoid conflicts
% with the MATLAB function 'conv'.
%
% PURPOSE:
% Computes the convolution of two vectors. 
% 
% FORMAT:
% c = convgauss( b, x, f, l);
% 
% INPUT: 
% b    Nx1 vector. 
% x    Lx1 vector. 
% f    scalar, the first convolution to compute. 
% l    scalAr, the last convolution to compute. 
% 
% OUTPUT: 
% c    Qx1 result, where Q = (l - f + 1) . 
% 
% If f is 0, the first to the l’th convolutions are computed. If l is 0,
% the f’th to the last convolutions are computed. If f and l are both zero,
% all the convolutions are computed. 
% 
% 
% REMARKS:
% If x and b are vectors of polynomial coefficients, this is the same as
% multiplying the two polynomials. 
%
% EXAMPLE:
% 
%      x = { 1,2,3,4 };
%      y = { 5,6,7,8 };
%      z1 = conv(x,y,0,0);
%      z2 = conv(x,y,2,5);
% 
%                                   5 
%                                  16 
%                                  34 
%                           z1 =   60 
%                                  61 
%                                  52 
%                                  32 
%   
% 
%         
%                                  16 
%                           z2 =   34 
%                                  60 
%                                  61 
%   
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde


if nargin == 2
    f=0;
    l=0;
end

w = conv(b,x);

if f~=0 && l~=0
    c = w(f:l);
elseif f==0 && l~=0
    c = w(1:l);
elseif f~=0 && l==0
    c = w(f:end);
elseif f==0 && l==0
    c=w;
end
