%% Computes the log10 marginal likelihood of the VAR with an independent
%  normal-Wishart prior.
%
% @author:  Roland Meeks
% @date:    December, 2012
%
% @return:  [scalar] log10 marginal likelihood for the model
% @arguments:
%           [matrix] Y data on dependent variables
%           [matrix] X data on independent variables
%           [struct] prior coefficients
%           [struct] posterior coefficient draws
%           [struct] [optional] options settings for calculation
%
% The optional argument is to accomodate an options struct
% Options:
% * Parameter vector at which to evaluate the BMI
% * Compute Monte Carlo standard error for mlik
function [ mlik se_lik ] = GetVarLogMarginalLikelihood( Y, X, prior, post, varargin )

%% House keeping
num_iter = size( post.A_post, 1 );
[ T M ] = size( Y );
narg = size( varargin, 2 );
cL = 10; % default choice for number of lags at which to truncate the 
         % covariance function when computing the MC standard error
         
%% Default choices of A*, Sigma*
% Get the coarse maximum likelihood values A*, Sigma*
[ loglik, index ] = max( post.loglik );
A_star = squeeze( post.A_post( index, :, : ) );
Sigma_star = squeeze( post.SIGMA_post( index, :, : ) );
H_star = inv( Sigma_star );

%% Optional arguments, if present, are in a struct

if ( 0 < narg )
    
    options = varargin{1};
    
    % check whether optional approximation point is supplied - if so, apply
    % these user-defined values instead of the defaults
    if ~( isempty( options.A_star ) || isempty( options.Sigma_star ) )
        
        % Check the user-supplied matrices are the right size
        if( all(size( options.A_star ) == size( A_star) ) && ...
                     all(size( options.Sigma_star ) == size( Sigma_star )) )
                    
            A_star = options.A_star;
            Sigma_star = options.Sigma_star;

            % Log likelihood:
            % If optional A*, Sigma* were supplied, get the value of the log
            % likelihood
            params = {eye(M,M), A_star, Sigma_star};
            loglik = GetVarLogLikelihood( Y, X, params );
        else
            warning('GetVarLogMarginalLikelihood::options', ...
               'A* and Sigma* are incorrectly specified: Using defaults' );
        end;
    end;
end;
% end options

%% Log prior
% The value of the prior density for A at A*
pA_prior = ln( mvnpdf( ...
                vec( A_star )', vec( prior.A_prior )', prior.V_prior ...
             ));

% The value of the prior density for Sigma^{-1} (denoted H) at Sigma*
pH_prior = wishartpdf( H_star, inv(prior.S_prior), prior.v_prior );

% The log joint prior density at A*, Sigma*
log_prior = pA_prior + pH_prior;

%% Log posterior
% The conditional log posterior for A, evaluated at Sigma*
A_prec = inv( prior.V_prior ) + ...
    kron( eye(M), X )'*kron( H_star, eye(T) )*kron( eye(M), X );
A_mean = A_prec\( prior.V_prior\vec(prior.A_prior) + ...
    kron( eye(M), X )'*kron( H_star, eye(T) )*vec( Y ) );
try
    % Sometimes mvnpdf throws an exception because of numerical inaccuracy
    % when inverting the covariance matrix
    A_covar = inv( A_prec );
    log_pA_post = ln( mvnpdf( vec( A_star )', A_mean', A_covar ) );
catch exception
    [ UT, p ] = chol( A_covar );
    A_variance = UT'*UT;
    log_pA_post = ln( mvnpdf( vec( A_star )', A_mean', A_variance ) );
end;

% The marginal log posterior for Sigma, evaluated at Sigma^*
pH_post = zeros( num_iter, 1 );

for iter=1:num_iter,
    A_draw = squeeze( post.A_post( iter, :, : ) );
    S_post = prior.S_prior + (Y - X*A_draw)'*(Y - X*A_draw);
    pH_post(iter) = wishartpdf( H_star, inv( S_post ), ...
        prior.v_prior + T );
end;

% #BUGBUG for very small values of pH_post, exp could be numerically zero,
% so the mean returns zero, and the logarithm returns -Inf...
log_pH_post_mean = ln( mean( exp( pH_post ) ) );

% The estimated log posterior density
log_posterior = log_pA_post + log_pH_post_mean;

%% The log marginal likelihood, using the relation that 
% log_10(x) = ln(x)/ln(10)
mlik = ( loglik + log_prior - log_posterior )/log(10) ;

%% Get the numerical standard error of the estimated marginal likelihood
MM = num_iter - cL;

% There is one variance and L covariances to store
covpostH = zeros(cL+1,1);

% takes account of the fact that the running index for the covariances
% starts at 1 (for c_0, the variance) rather than zero, and so runs up to
% L+1
for j = 1:cL+1,
    % recall, these are the pH_post draws are for the *log* of the Wishart
    % density at H*
    mcov = cov( exp(pH_post( j:MM+j-1 )), exp(pH_post( 1:MM )), 1 );
    covpostH(j) = mcov(2,1); % pick out the off-diagonal
end;

% The MCSE for the posterior density of H evaluated at H*:
tau2 = covpostH(1) + 2*( 1 - [1:cL-1]/cL )*covpostH(2:cL);

% By the delta method
se_lik = sqrt( tau2/exp( 2*log_pH_post_mean ) )/log(10);

%% Optional display results: mlik (se) and 95% CI
try
    if (options.print_results)
        dispstr = [ 'Log10 marginal likelihood = ' num2str(mlik) ...
            ' (' num2str(se_lik) ')' ];
        fprintf( '\n%s', dispstr );
        dispstr = [ '95% CI = (' num2str(mlik - 1.96*se_lik) ', ' ... 
            num2str(mlik + 1.96*se_lik) ')' ];
        fprintf( '\n\t\t\t\t%s\n', dispstr );
    end;
catch err
    if strcmp( err.identifier, 'MATLAB:nonExistentField' )
        fprintf( '%s', '' );
    else
        rethrow( err );
    end; % end if
end; % end try-catch

end