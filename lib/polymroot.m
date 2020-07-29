function [a] = polymroot(c)
% polymroot.m
% 
% Purpose:    compute the roots of the determinant of a matrix polynomial
% 
% Format:     r = polymroot(c);
% 
% Input:      c      (N+1)*KxK matrix of coefficients of an Nth order
%                    polynomial of rank K.
% 
% Output:     r      K*N vector containing the roots of the determinantal
%                    equation
% 
% Remarks:  c is constructed from N KxK coefficient matrices (e.g. phi)
%           stacked vertically. The highest order coefficient matrix,
%           phi(p), is at the top, phi(p-1) is appended to the bottom of
%           phi(p), etc.
% 
%           An identity matrix is often appended to the bottom, to reflect the
%           matrix polynomial I - phi(1)z - phi(2)z^2 - ... - phi(p)z^p.
%           Note that in this case the phi matrices need to be multiplied by -1.
% 
%           Note that this procedure solves the scalar problem as well, that
%           is, the one that POLYROOT solves.
% 
% Example:  Solve det(A2*t^2 + A1*t + A0) = 0 where
% 
%             A2 = [ 1  2 ]
%                  [ 2  1 ]
% 
%             A1 = [ 5  8 ]
%                  [10  7 ]
% 
%             A0 = [ 3  4 ]
%                  [ 6  5 ]
% 
% 
%             a2 = { 1 2, 2 1 };
%             a1 = { 5 8, 10 7 };
%             a0 = { 3 4, 6 5 };
% 
%             print polymroot(a2|a1|a0);
% 
%                -4.3027756
%                -.69722436
%                -2.6180340
%                -.38196601
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde

n0 = rows(c)/cols(c) - 1;
n = cols(c);   
if n0 == 0;
    a = zeros(n,1);
    return;
else
    a = -c(n0*n+1:size(c,1),:)/c(1:n,:);
    if sum(find(isnan(a)))~=0
        error('Missing values');
    end
    if n0 > 1
        a1 = [];
        for i = 1:n0-1
            a0 = -c(i*n+1:(i+1)*n,:)/c(1:n,:);
            if sum(find(isnan(a0)))==0
                a1 = [a1  a0]; %#ok<AGROW>
            else
                error('Missing values!');
            end
        end
        a = [[a1  a] ; [eye(n*(n0-1)) zeros(n*(n0-1),n)]];
    end
    if sum(find(isnan(a)))==0
        a = eig(a);
    else
        error('Missing values!');
    end
    return;
end