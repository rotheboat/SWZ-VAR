function log_prior_density_at_point = GetPriorDensityAtParamVector( szprior, cbstar, cgstar )
% Given any parameter vectors (vbstar, vgstar), evaluate the prior density
% function parameterized by matrices stored in the szprior structure.
% Author: Roland Meeks
%
% @return:  [scalar] natural logarithm of the prior density at a point

cStld = szprior.Stld;
ciStld = szprior.Stld_inv;
ciHtld = szprior.Htld_inv;
cPtld = szprior.Ptld;

% 1. Density of the b coefficients
    % (a) The sum of the log absolute determinants of S
    % Use cellfun to apply the operation log(|det(.)|) to each element of
    % the cell array of prior covariance matrices Stld.  The output of 
    % cellfun in this case is an array of scalars, which are then summed.
    % * Note i: as S is a covariance matrix it is PD, and so has a
    % strictly positive determinant.
    % * Note ii: Uses the logdet()method to help avoid over/underflow. See:
    % http://tinyurl.com/zh6q47d
    sumLogDetS = sum( cellfun(@(S) logdet( S, 'chol' ), cStld, 'UniformOutput', true ) );
    
    % (b) The sum of the quadratics b'/S*b
    sumQuadB = sum( cellfun( @(b,iS) b'*iS*b, cbstar, ciStld, 'UniformOutput', true ) );
    
    % (c) Get the number of unrestricted entries in A in each equation, as 
    % measured by the number of columns in U, and sum them.
    sumq = sum( cellfun( @(U) size( U,2 ), szprior.cU, 'UniformOutput', true ) );
    
    % (d) The log conditional prior density evaluated at cbstar
    bdensity = -sumq*log(2*pi)/2 - sumLogDetS/2 - sumQuadB/2;

% 2. Density of the g coefficients
    % (a) See 1(a) above. The presence of a leading -1 on the method logdet
    % reflects that the prior struct contains a cell array of inv(H) 
    % matrices.  The determinant:
    %   det( inv(H) ) = 1/det(H)
    % therefore:
    %   log( det(H) ) = -log( det( inv(H) ) ).
    sumLogDetH = sum( cellfun(@(iH) -logdet( iH, 'chol' ), ciHtld, 'UniformOutput', true ) );
    
    % (b) The sum of the quadratics (g-P*b)'/H*(g-P*b)
    sumQuadG = sum( cellfun( @(b,g,iH,P) (g - P*b)'*iH*(g - P*b), ...
            cbstar, cgstar, ciHtld, cPtld, 'UniformOutput', true ) );
    
    % (c) Get the number of unrestricted coefficients in F in each 
    % equation, as measured by the number of columns in V, and sum them.
    sumr = sum( cellfun( @(V) size( V,2 ), szprior.cV, 'UniformOutput', true ) );    
    
    % (d) The log conditional prior density evaluated at (cbstar, cbstar)
    gdensity = -sumr*log(2*pi)/2 - sumLogDetH/2 - sumQuadG/2;

% Final return value
log_prior_density_at_point = bdensity + gdensity;

end