function [y] = polychar(x)
% polychar.m
% 
% Purpose:    To compute the characteristic polynomial of a square
%             matrix.
% 
% Format:     c = polychar(x);
% 
% Input:      x    NxN matrix.
% 
% Output:     c    N+1x1 vector of coefficients of the Nth order
%                  characteristic polynomial of x:
% 
%                  p(z)=c[1,1]*z^n + c[2,1]*z^(n-1) + ... +
%                  c[n,1]*z + c[n+1,1];
% 
% Remarks:    The coefficient of z^n is set to unity (c[1,1]=1).
% 
% See Also:   polymake, polymult, polyroot, polyeval

y = polymake(eig(x));