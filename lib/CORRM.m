function [cx] = corrm(m)
  
% PURPOSE:
% Computes correlation matrix from moment (x'x) matrix.
% 
% FORMAT:
% cx = corrm(m);
% 
% INPUT:
% m     KxK moment (x'x) matrix.
%
% REMARKS: 
% Constant term MUST have been the first variable when moment matrix was computed.
% This version does NOT support arrays. m is only a 2-D matrix.
%
% OUTPUT:
% cx    KxK correlation matrix.
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde



% Check for complex input */
if isreal(vc) == 0;       
    error('ERROR: Not implemented for complex arguments.')
end

% if (type(m) == 6)      
    cc = seqa(2,1,cols(m)-1);
    xx = m(cc,cc);              % Pull out K-1xK-1 submatrix
    n = m(1,1);                 % Number of observations
    xbar = m(cc,1)/n;           % Vector of means
    vv = xbar*xbar';
    vc = (xx-n*vv)/(n-1);       % VC matrix
    std = sqrt(diag(vc));
    cx = vc./(std.*std');
% elseif (type(m) == 21)    
%     dims = getdims(m);
%     orders = getorders(m);    
%     if (dims > 2 )
%         index = ones(dims-2,1);
%         neworders = orders[1:dims-2];
%         neworders = neworders||orders[dims-1]-1||orders[dims]-1;
%         a = arrayalloc(neworders,0);
%         
%         loopni:
%         setarray a, index, corrm(getmatrix(m,index));
%         loopnextindex loopni, index, orders;
%     else
%         if (dims == 1);
%             neworders = orders-1;    
%         else
%             neworders = orders[dims-1]||orders[dims];
%             endif;
%             a = mattoarray(corrm(arraytomat(m)));
%             endif;
%             
%             retp(a);
%         else
%             error('ERROR: Type mismatch.');
%         end
%     end
% end


