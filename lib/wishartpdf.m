% Evaluates the logarithm pdf for the Wishart distribution for a given PSD 
% matrix
%
% @input:  W = matrix to be evaluated
%          S = parameter - a PD matrix s.t. RV has expectation S*m
%          m = degrees of freedom \ge dim(S)
% @output: [scalar] = logarithm of the Wishart pdf at W
function pw = wishartpdf( W, Sigma, m )

%% Perform checks on method arguments
[ pp qq ] = size( W );
[ rr ss ] = size( Sigma );

if ~( pp == rr && qq == ss )
    error( 'wishartpdf::IncompatibleArguments', ...
        'The supplied matrices are not conformable.' );
else
    p = pp;
end;

try
    % Check if matrix is PD (TRUE if p == 0)
    [R, pd] = chol( Sigma );

    if ~( 0 == pd )
        warning( 'wishartpdf::SigmaMatrix', ...
                    'Non-positive definite matrix passed.' );
    end;

    % Check matrix is square, symmetric
    if ~all( Sigma == Sigma' )
        error( 'wishartpdf::SigmaMatrix', 'Non-symmetric matrix passed.' );
    end;
catch exception
    % make Sigma symmetric PD
    Sigma = R'*R; % R is upper triangular
end;

  

%% Compute logarithm of pdf; constant of integration is c
log_c = m*p*ln(2)/2 + m*ln(det(Sigma))/2 + p*(p-1)*ln(pi)/4 + ...
                sum(ln( gamma((m+1-[1:p])/2)));

pw = -log_c + (m-p-1)*ln( det( W ) )/2 - trace( Sigma\W )/2;
        
end        