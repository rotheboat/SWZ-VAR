function [ log10mdd ] = GetSimsZhaMDD( mY, mX, szprior, postest, options, varargin )
% Given the output of the Waggoner-Zha Gibbs sampler, returns an estimate
% of the marginal data density computed using the Chib (JASA, 1995)
% approach.
%
% @return:  log_10 marginal data density (marginal likelihood)

% Step 0: process options and set parameters

    % Step 0.1 Process options
    nVarargs = length(varargin);

    if 2 == nVarargs
        Astar = varargin{1};
        Fstar = varargin{2};
    end
    
    % Step 0.2 Set parameters...
    M = size( mY, 2 );
    K = size( mX, 2 );
    
    % Step 0.3 Check for other relevant settings
    if( ~isfield( options, 'reject_nonstationary' ) )
        reject_nonstationary = false;
        prob_stationary = 1;
    else
        reject_nonstationary = options.reject_nonstationary;
        
        if reject_nonstationary
            total_draws = postest.reject_count + postest.draws;
            prob_stationary = 1 - postest.reject_count/total_draws;
        else
            prob_stationary = 1;
        end
    end
% Step 1: find the posterior mode (A^*, F^*)

    % Step 1.1: There are three possibilities here: (a) a coarse 
    % approximation, that takes the highest posterior draw to be the mode; 
    % (b) a numerical optimum, that uses a numerical optimizer to find the 
    % mode; (c) we are already passed a value as an optional argument.

    % Change the (local) option to make GibbsBVAR initialize at Astar
    options.mA = Astar;
    % As we are supplying an impact matrix GibbsBVAR will set
    % number_independent_chains = 1.  We should therefore reset the number
    % of draws appropriately
    options.draws = options.number_independent_chains*options.draws;
    options.burn = options.number_independent_chains*options.burn;
    
    % Step 1.2: Convert the matrices of VAR parameters A^* and F^* into
    % cell arrays containing b^* and g^*, the vector parameters in each of
    % the M equations of the model under the restrictions.
    
        % Impact matrix A -> cell array b = U'*a
        ca = mat2cell(Astar,M,ones(1,M));
        cbstar = cellfun(@(a,U) U'*a, ca, szprior.cU, 'UniformOutput', false );
        % Lag matrix F -> cell array g = V'*f
        cf = mat2cell(Fstar,K,ones(1,M));
        cgstar = cellfun(@(f,V) V'*f, cf, szprior.cV, 'UniformOutput', false );

% Step 2: evaluate the ln posterior density at A^*, F^*
    % create storage
    log_marginal_posterior_b = NaN(postest.draws,M);

    % Step 2.1: p( b_1* | Y, X )
    % Evaluate the ln marginal posterior density for the first column
    % of A, a_1 = U_1*b_1, at the modal value a^*_1 using the output of the 
    % full Gibbs sampler.
    for d = 1:postest.draws
        equation = 1;
        postest.A_post(d,:,1:equation) = Astar(:,1:equation);
        mA = squeeze(postest.A_post(d,:,:));
        % Impact matrix A -> cell array a
        ca = mat2cell(mA,M,ones(1,M));
        % cell array columns a -> cell array b = U'*a
        cb = cellfun(@(a,U) U'*a, ca, szprior.cU, 'UniformOutput', false );        
        % evaluate marginal posterior density at b
        log_marginal_posterior_b(d,1) = ...
            GetMarginalPosteriorDensityAtParamVector( ...
                szprior, postest, equation, mA, cb );
    end

    %#BUGBUG ONLY IF M > 2#
    
    % Step 2.2: p( b_2* | b_1*, Y, X )
    % Evaluate the ln marginal posterior density for the 1 < i < m
    % column of A, a_i = U_i*b_i, at the modal value a^*_i using the output of
    % the reduced Gibbs sampler.
    for equation = 2:M-1
        % The reducedGibbs option tells GibbsBVAR where to start sampling.
        % Also the presence of this option tells GibbsBVAR to produce draws
        % from the reduced Gibbs sampler.
        options.reduced_gibbs = equation;
        
        % Execute a reduced Gibbs run
        postest = GibbsBVAR( mY, mX, szprior, options );
        
        % For each draw in the reduced run, evaluate the posterior density
        % for equation i as:
        %   p( b*_i | b*_1,...,b*_{i-1}, b(d)_{i+1},...,b(d)_M ) 
        % and store in column i of the array "log_marginal_posterior_b"
        for d = 1:postest.draws
            postest.A_post(d,:,1:equation) = Astar(:,1:equation);
            mA = squeeze(postest.A_post(d,:,:));
            % Impact matrix A -> cell array a
            ca = mat2cell(mA,M,ones(1,M));
            % cell array columns a -> cell array b = U'*a
            cb = cellfun(@(a,U) U'*a, ca, szprior.cU, 'UniformOutput', false );
            % evaluate marginal posterior density at b
            log_marginal_posterior_b(d, equation) = ...
                GetMarginalPosteriorDensityAtParamVector( ...
                    szprior, postest, equation, mA, cb );
        end
    end

    % Step 2.3: p( b_M* | b_1*,...,b_{M-1}*, Y, X )
    % Analytically evaluate the ln marginal posterior density for
    % column m of A, a_m = U_m*b_m.
    log_marginal_posterior_b(:,M) = ...
        GetMarginalPosteriorDensityAtParamVector( szprior, postest, M, Astar, cbstar );

    % Step 2.4: compute the mean across draws, then sum
    mean_log_marginal_posterior_b = mean( log_marginal_posterior_b );
    log_marginal_posterior_b_star = sum( mean_log_marginal_posterior_b );
    
    % Step 2.5: evaluate the ln marginal posterior density for the right hand 
    % side coefficients g_i at the posterior mode values (b^*, g^*) in 
    % equations i = 1,...,m.
    log_conditional_posterior_at_point = ...
        GetConditionalPosteriorDensityAtParamVector( postest, cbstar, cgstar );
    

% Step 3: evaluate the ln prior density at the posterior mode (b^*, g^*)
log_prior_density_at_point = GetPriorDensityAtParamVector( szprior, cbstar, cgstar );

    % Step 3.1: if there is an option to reject non-stationary draws, form
    % an estimate of the probability of being in that region of the
    % parameter space from the posterior draws
    if reject_nonstationary
        log_prior_density_at_point = ...
            log_prior_density_at_point - log( prob_stationary );
    end

% Step 4: evaluate the ln data density at the posterior mode (b^*, g^*)
log_data_density_at_point = GetDataDensityAtParamVector(mY, mX, Astar, Fstar);


% Step 5: using the outputs of the computations above, compute the Chib
% approximation to the ln marginal data density, according to Chib's BMI
% (Chib 1995, eq. 6).  Convert to base 10 logarithms by using the relation
% log_10(x) = ln(x)/ln(10).

% BUGBUG TEST
% Return an estimate of the log marginal data density (marginal likelihood)
% at the point (A*,F*).
log10mdd = ( log_prior_density_at_point ...         % prior density
           + log_data_density_at_point ...          % data density
           - log_conditional_posterior_at_point ... % posterior g density
           - log_marginal_posterior_b_star ...      % posterior b density
            )/ln(10);                               % convert to base 10

end


function rpostest = GetPosteriorOrdinateEstimate()
% Simulates draws from the posterior density of the Wagonner-Zha model when
% some blocks of parameters take on fixed values.  For example, if we have
% a full posterior p(b_1, b_2, b_3|y) for which all the conditional 
% distributions are known, and we wish to simulate draws from 
% p(b_1|b_2^*,y), we run a `reduced' Gibbs sampler on p(b_1|b_2^*,b_3,y) 
% and p(b_3|b_1,b_2^*,y) to obtain these.

end


function log_marginal_posterior_b = GetMarginalPosteriorDensityAtParamVector( szprior, postest, equation, mUb, cb )
% Evaluates the marginal posterior density p(b_i|b_{-i},Y,X) at a given point b^,
% returning the log.

% Number of variables in the VAR equals the number of equations, equals the
% number of restriction matrices in the cell array cU
M = length( szprior.cU );

% Catch input errors
assert( isequal( [M,M], size( mUb ) ), 'Expecting an M x M matrix in argument 4' );
assert( isequal( M, length( cb ) ), 'Expecting an M dimensional cell array in argument 5' );
assert( isequal( 'cell', class( cb ) ), 'Expecting a cell array in argument 5' );

% Step 0: set up constants
q_i = size( cb{equation}, 1 );
Tobs = postest.Tobs;
cH = postest.cH;
cT = postest.cT;
cU = szprior.cU;
cS_inv = postest.cS_inv;

% Step 1: compute the log normalizing constant
% To avoid over/underflow the function gammaln is used to compute the
% logarithm of the gamma function.
const_density = ...
        (Tobs-1+q_i)*log(Tobs/2)/2 - q_i*log(pi)/2 - gammaln( (Tobs+1)/2 ) ...
                + q_i*log( Tobs/(2*pi) )/2 - logdet( cH{equation}, 'chol' );

% Step 2: compute the log determinant of T_i
sumLogDetT = sum( cellfun( @(T) logdet(T), cT, 'UniformOutput', true ) );
    
% Step 3: compute the log of w'*U*b
mJ = mUb;
mJ(:,equation) = 0;
sw = null( mJ' );
logwUb = log( abs( sw'*cU{equation}*cb{equation} ) );

% Step 4: compute the quadratic b'*S*b
quadB = cb{equation}'*cS_inv{equation}*cb{equation};

% Return the log marginal posterior
log_marginal_posterior_b = const_density - sumLogDetT + Tobs*logwUb - Tobs*quadB/2;

end


function log_conditional_posterior_g = GetConditionalPosteriorDensityAtParamVector( postest, cbstar, cgstar )
% Evaluates the conditional posterior density p(g|b,Y,X) at a given point
% (g^, b^), returning the log.

% Step 1: compute the number of unrestricted coefficients in F in each 
    % equation, as measured by the number of elements in the g(i) vectors 
    % stored in the cell array cg.
    sumr = sum( cellfun( @(g) size( g,1 ), cgstar, 'UniformOutput', true ) );
    
% Step 2: compute the sum of the log absolute determinants of H, the
    % conditional posterior covariance matrix
    %
    % Use cellfun to apply the operation log(|det(.)|) to each element of
    % the cell array of prior covariance matrices Stld.  The output of 
    % cellfun in this case is an array of scalars, which are then summed.
    %
    % * Note i: as H is a covariance matrix it is PD, and so has a
    % strictly positive determinant.
    % * Note ii: Uses the logdet()method to help avoid over/underflow. See:
    % http://tinyurl.com/zh6q47d
    sumLogDetH = sum( cellfun(@(H) logdet( H, 'chol' ), postest.cH, 'UniformOutput', true ) );

% Step 3: compute the sum of the quadratics (g-P*b)'/H*(g-P*b)
    sumQuadG = sum( cellfun( @(b,g,H,P) (g - P*b)'/H*(g - P*b), ...
            cbstar, cgstar, postest.cH, postest.cP, 'UniformOutput', true ) );
    
% Return log conditional posterior
log_conditional_posterior_g = -sumr*log(2*pi)/2 - sumLogDetH/2 - sumQuadG/2;

end
