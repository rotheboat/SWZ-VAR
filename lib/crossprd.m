function z = crossprd(x,y)
% Purpose:    Computes the cross products (vector products) of
%             sets of 3x1 vectors.
% 
% Format:     z = crossprd(x,y);
% 
% Input:      x    3xK matrix, each column is treated as a 3x1 vector.
% 
%             y    3xK matrix, each column is treated as a 3x1 vector.
% 
% Output:     z    3xK matrix, each column is the cross product
%                  (sometimes called vector product) of the
%                  corresponding columns of x and y.
% 
% Remarks:    The cross product vector (z) is orthogonal to both x and y.
%             sumc(x.*z) and sumc(y.*z) will be Kx1 vectors all of whose
%             elements are 0 (except for rounding error).

      
r1 = x(2,:).*y(3,:)-x(3,:).*y(2,:);
r2 = x(3,:).*y(1,:)-x(1,:).*y(3,:);   
r3 = x(1,:).*y(2,:)-x(2,:).*y(1,:);
z = [r1 ; r2 ; r3];

