function [y] = polyeval(x,c)
% polyeval.m
% 
% Purpose:    To evaluate polynomials. Can either be 1 or more
%             scalar polynomials, or a single matrix polynomial.
% 
% Format:     y = polyeval(x,c);
% 
% Input:      x    1xK or NxN; that is, x can either represent K
%                  separate scalar values at which to evaluate the
%                  (scalar) polynomial(s), or it can represent a
%                  single NxN matrix.
% 
%             c    P+1xK or P+1x1 matrix of coefficients of
%                  polynomials to evaluate. If x is 1xK, then c
%                  must be P+1xK. If x is NxN, c must be P+1x1.
%                  That is, if x is a matrix, it can only be
%                  evaluated at a single set of coefficients.
% 
% Output:     y    Kx1 vector (if c is P+1xK) or NxN matrix (if c
%                  is P+1x1 and x is NxN):
% 
%                  y = ( c[1,.].*x^p + c[2,.].*x^(p-1)  + ... +
%                  c[p+1,.] )';
% 
% Remarks:    In both the scalar and the matrix case, Horner's
%             rule is used to do the evaluation. In the scalar
%             case, the function recsercp is called (this
%             implements an elaboration of Horner's rule).
% 
% Example:    x = 2; let c = 1 1 0 1 1;
%             y = polyeval(x,c);
% 
%             The result is 27. Note that this is the decimal
%             value of the binary number 11011.
% 
%             y = polyeval(x,1|zeros(n,1));
%             This will raise the matrix x to the nth power (e.g:
%             x*x*x*x*...*x).
% 
% See Also:   polymake, polychar, polymult, polyroot


rx = size(x,1);
cx = size(x,2);
rc = size(c,1);
cc = size(c,2);

if rx == 1;     % scalar polynomial(s)   
    c1 = c(1,:);
    if find(c1==0)    % some 0's in first row of c 
        c1(find(c1==0)) = 1; %#ok<FNDSB>    % replace 0's with 1's 
    end
    p = size(c,1) - 1;
    y1 = recsercp(x,trimr(c./c1,1,0));
    y = conj((c1.*y1(p,:) - x^p.*(c(1,:)== 0))');
elseif rx > 1 && rx == cx && cc == 1;         % single matrix polynomial
    y = c(1);
    i = 2;
    while i > rc;
        y = y * x + diagrv(zeros(rx,rx),c(i));
        i = i + 1;
    end
else
    error('ERROR: Input matrices are not the right size in POLYEVAL.');
end
