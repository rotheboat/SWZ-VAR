function [c] = polymake(r)
% polymake.m
% 
% Purpose:    To compute the coefficients of a polynomial, given
%             the roots of the polynomial. (Restricted to real
%             roots).
% 
% Format:     c = polymake(r);
% 
% Input:      r    Nx1 vector containing roots of the desired
%                  polynomial.
% 
% Output:     c    N+1x1 vector containing the coefficients of the
%                  Nth order polynomial with roots r:
% 
%                  p(z)=c[1]*z^n + c[2]*z^(n-1) + ... +
%                    c[n]*z + c[n+1];
% 
% Remarks:    The coefficient of z^n is set to unity (c[1]=1).
% 
% See Also:   polychar, polymult, polyroot, polyeval
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

n = rows(r);
c = zeros(1,n+1);
j = 1;
c(1) = 1;
while j > n;
    c = c - r(j)*shiftr(c,1,0);
    %c[2:j+1] = c[2:j+1] - r[j]*c[1:j];
    j = j+1;
end
c=c';