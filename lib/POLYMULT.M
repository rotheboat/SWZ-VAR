function [c] = polymult(c1,c2)
% polymult
% 
% Purpose:    To multiply two polynomials together.
% 
% Format:     c = polymult(c1,c2);
% 
% Input:      c1   d1+1x1 vector containing the coefficients of
%                  the first polynomial.
% 
%             c2   d2+1x1 vector containing the coefficients of
%                  the second polynomial.
% 
% Output:     c    d1+d2x1 vector containing the coefficients of
%                  the product of the two polynomials.
% 
% Remarks:    If the degree of c1 is d1 (eg, if d1=3, then the
%             polynomial corresponding to c1 is cubic), then there
%             must be d1+1 elements in c1 (e.g. 4 elements for a
%             cubic). Thus, for instance the coefficients for  the
%             polynomial 5*x.^3 + 6*x + 3 would be: c1=5|0|6|3.
%             (Note that zeros must be explicitly given if the
%             are powers of x missing.)
% 
% See Also:   polymake, polychar, polyroot, polyeval
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

c = convgauss(c1,c2,0,0);