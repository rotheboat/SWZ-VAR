%% 
% Roland Meeks (2010-11)
%
% Given information on independent normal-Wishart prior distributions, OLS 
% parameter estimates, and data, uses Gibbs sampling to draw from the 
% posterior distributions B|Sigma and Sigma|B. Returns the draws, and
% functions of them such as impulse-responses.
function postest = GibbsINW( priorparams, strols, p, post_options )
%% ====================== Initialization & Storage ========================

% Data & important constants of dimension
Y = strols.Y;
X = strols.X;
[ T M ] = size( Y );
K = size( X, 2 );
Z = kron( eye( M ), X );

% Check if BVAR includes an intercept
if ( K > M*p )
    constant = 1;
elseif ( K == M*p )
    constant = 0;
end;
%% BUGBUG ******************
zero_select = ones( K, M );
zero_select( 11:12, 1:7 ) = 0;
zero_select( 22:23, 1:7 ) = 0;
zero_select( [2:10 12], 10 ) = 0;
zero_select( [13:21 23], 10 ) = 0;
zero_select = ~logical( zero_select );

%% END BUGBUG

% coefficient matrices from prior
f_prior = vec( priorparams.A_prior );
V_prior = priorparams.V_prior;
H_prior = inv( V_prior );

% The naming of the fields preserves notation used in GetINWishartPosterior
% but note that the A matrix is denoted F here, following Waggoner-Zha.
postest = struct( 'A_post', [], 'A_mean', [], 'V_post', [], 'v_post', [], ...
    'SIGMA_post', [], 'SIGMA_mean', [], 'IRF', [], 'loglik', [], 'is_restricted', [] );

% ---- Set restrictions ----------------------------------
% Default restrictions:
% For mA, the default is a triangular normalization
cU = cell(1,M);
for i=1:M
    cU{i} = eye(M,i);
end
% For mF, the default is no restrictions
mycell = cell(1,M); 
cV = cellfun(@(n)eye(K,K),mycell,'UniformOutput', false);

% My standard approach is that mRestrict is a 1x2 cell array.
% 	mRestrict{1} is a M-by-M matrix with ones where mA has a non-zero coefficient; zeros elsewhere
% 	mRestrict{2} is a K-by-M matrix with ones where mF has a non-zero coefficient; zeros elsewhere
% Either can be empty, meaning that the defaults defined above are used (e.g. 
% restrictions on A, but no restrictions on the lags in any equation). The call to 
% processRestrictions() takes the basic input in mRestrict and converts it to U and  
% V matrices that have the property a(i) = U(i)b(i) and f(i)=V(i)g(i), where b(i)
% and g(i) are the restricted parameter vectors.
if isfield( priorparams, 'restrictions' )
    mRestrict = priorparams.restrictions;
    
    % check if restrictions are empty
    postest.is_restricted = ...
                ~all(cellfun( @(n)isempty(n), mRestrict, 'UniformOutput', true ));    
    if postest.is_restricted
        [cU, cV] = processRestrictions( mRestrict, cU, cV );
        mRestrict = {cU, cV};        
        fprintf( '\nApplying parameter restrictions...' );
    end;
end;


% Initialize Bayesian posterior parameters using OLS values...
% Note that as we first draw from the conditional for f, we don't need
% to get an 'initial' f (if we were first drawing SIGMA, we would).
SIGMA = strols.SIGMA_OLS*(T-K)/(T-K+1); % n.b. check df?
H_temp = inv(SIGMA);

% Can drop strols
% clear strols;

% Storage space for posterior draws
F_draws = zeros(post_options.draws,K,M);
V_draws = zeros(post_options.draws,M*K,M*K);
SIGMA_draws = zeros(post_options.draws,M,M);
loglik_draws = zeros( post_options.draws, 1 );
imp = zeros( post_options.draws, M, M, post_options.irf_horizon );

% Gibbs-related preliminaries
ntot = post_options.draws + post_options.burn;  % Total number of draws

%% ========================= Start Sampling ===============================
tic;

disp('Number of iterations');

drawcount = 0 ;
rejectcount = 0 ;
flag = 1;

% Gibbs sampler
while ( drawcount < ntot+1 )
    
    if ( drawcount + rejectcount > post_options.max_draws )
        warning( 'GetINWishartPosterior:GibbsLoop', ...
            'Maximum number of draws exceeded.' );
        break;
    end;
    
    if ( (1 == flag) && (0 == mod(drawcount,post_options.it_print) ) )
        fprintf( '\nIteration = %d\n', drawcount );
        toc;
    end
    % Display a dot every it_print/10 draws...
    if ( 0 == mod(drawcount,ceil(post_options.it_print/10)) )
        fprintf( '.' );
    end
    
    %% (1) The F|SIGMA block -------------------------------------------
    VARIANCE = kron(H_temp,eye(T)); %Note that this is inconsistent
    % with equation 2.19 in Koop and Korobilis, since they adopt a
    % different notation for the vector of stacked observations y than the
    % one used in the code, and in Geweke (2005, p. 164)
    
    % Try-catch block will pick up if V_post is non-PD; if it is not then
    % we set the flag to zero and draw again. Most likely there is a
    % numerical problem with inv(). Note that if a matrix is non-PD you can
    % still get chol() to return without an error, but the resulting UT
    % matrix will not be of the expected dimension (or of course, rank).
    try
        % These H's are the inverse covariance of the parameters, denoted
        % H_beta by Geweke.  The H_temp variable, by contrast, is the
        % inverse of the covariance matrix of the errors
        H_post = H_prior + Z'*VARIANCE*Z;
        f_post = H_post\( (V_prior\f_prior) + Z'*VARIANCE*Y(:)); % RM replaced inv(V)*a
        V_post = H_post\eye( size( H_post ) ); % Replaced inv()
        
        % if check_pd is 0, then the V_post matrix is indeed PD
        [ R check_pd ] = chol( V_post, 'lower' );
        % if not PD, generate an exception which will be handled by the
        % first part of the catch block
        assert( 0 == check_pd, 'MATLAB:posdef', 'Parameter covariance not PD' );
        
        f = f_post + R*randn(K*M,1); % Draw of f     
        F = reshape(f,K,M); % Draw of F
        F( zero_select ) = 0; %% BUGBUG SHORTCUT

        % (2) Reject unstable draws ... see Cogley & Sargent ("Evolving
        % Post-WWII US Inflation Dynamics") for a formal justification of the
        % procedure below in terms of rejection sampling where the 'full'
        % posterior (what they call the 'virtual' posterior) density is the
        % proposal distribution.

        if ( post_options.reject_nonstationary )
            % The companion form of the VAR AR coefficient matrices...
            COMP = [ F( constant+1:end, : )'; eye( M*(p-1) ), zeros( M*(p-1), M ) ];
    
            % Check the roots of the companion matrix, noting draws that violate
            % stationarity condition
            if ( any( abs( eig( COMP ) ) > 1.001 ) )
                flag = 0;
                rejectcount = rejectcount+1;
            else
                flag = 1;
            end;
        else
            flag = 1;
        end;

    catch err
        if ( strcmp( err.identifier, 'MATLAB:posdef' ) )
            flag = 0;
            rejectcount = rejectcount+1; % RM 2014-03-25 Was missing-a bug
                                         % as failed to quit loop if many
                                         % non-PD matrices drawn
            
            % Possible solution: re-draw H with same S_post; notice that
            % S_post must be based on a stationary draw of F (if
            % post_options.reject_nonstationary is set to true); this is
            % merely an attempt to draw a valid H_post
            H_temp = wish(inv(S_post),v_post);
        else
            rethrow( err ); % There's some sort of serious problem
        end; % end if-else
    end; % end try-catch
    

    % --------------------------------------------------------------------
    % Draws that are non-stationary are rejected, as are draws for which
    % SIGMA is not positive definite. 
    % Previously code skipped the new draw of SIGMA whenever flag = 0, but
    % this allows the algorithm to become stuck 
    
    %% (3) Posterior of SIGMA|F,Data ~ iW(inv(S_post),v_post)
    if ( 1 == flag )
        v_post = T + priorparams.v_prior; % note: does not change on iteration
        S_post = priorparams.S_prior + (Y - X*F)'*(Y - X*F);
        % Draw SIGMA
        % RM 2014-03-25 replaced one line draw of SIGMA with two lines so
        % that H_temp is available above, removing the need to invert sigma
        % when computing VARIANCE
        H_temp = wish(inv(S_post),v_post);
        SIGMA = inv(H_temp);

        % (4) Store results for every iteration in excess of drawcount
        if drawcount > post_options.burn               
            % =========IMPULSE RESPONSES:
            %Btemp = reshape( F, constant+M*p, M );
            Btemp = F(constant+1:end,:);
            % note: chol( SIGMA ) delivers an UT matrix, thus SIGMA=UT'*UT
            imp_resp = impulse( reshape( Btemp', M, M, p ), ...
                                 chol( SIGMA ), post_options.irf_horizon );
            imp(drawcount-post_options.burn,:,:,:) = imp_resp ;

            % ----- Save draws of the parameters
            F_draws(drawcount-post_options.burn,:,:) = F;
            V_draws(drawcount-post_options.burn,:,:) = V_post;
            SIGMA_draws(drawcount-post_options.burn,:,:) = SIGMA;

            loglik_draws(drawcount-post_options.burn) = ...
                    GetVarLogLikelihood( Y, X, {eye(M,M), F, SIGMA} );
        end % end saving results
    else
        % We have a new F which is non-stationary, or a new SIGMA which
        % is not positive definite
        
    end; % end if non-stationary flag is zero
    
    drawcount = drawcount + flag;
end; %end the main Gibbs for loop

fprintf( '\nDone with Gibbs sampler\n' );
disp( ['Rejected ' num2str( 100*rejectcount/(rejectcount+drawcount) ) ...
                                        '% of draws.'] )

toc;
%% ===================== End Sampling Posteriors ==========================

% If drawcount < post_options.draws + post_options.burn because the number
% of rejected draws is high, need to truncate the arrays of posterior
% draws (which are initialized with zeros)
ndraws_actual = drawcount - post_options.burn - 1;

SS = zeros( M );
for j = 1:(drawcount-post_options.burn-1)
    SS = SS + squeeze( SIGMA_draws(j,:,:) );
end;


postest.A_post = F_draws(1:ndraws_actual,:,:);
postest.A_mean = squeeze( mean( F_draws(1:ndraws_actual,:,:) ) );
postest.SIGMA_post = SIGMA_draws(1:ndraws_actual,:,:);
postest.SIGMA_mean = SS/(drawcount-post_options.burn);
postest.IRF = imp(1:ndraws_actual,:,:,:);
postest.loglik = loglik_draws( 1:ndraws_actual );

end


% Subfunction to process restrictions.  This is called if a parameter
% vector has been passed that must be correctly mapped into structural form
% parameter matrices.  
function [cU, cV] = processRestrictions( mRestrict, cU, cV )
    
    A_zeros = mRestrict{1};
    F_zeros = mRestrict{2};
    M = size( cU{1}, 1 );
    K = size( cV{1}, 1 );
    
    if ~isempty( A_zeros )
        
        for i = 1:M,
            % The more obvious way to do create Q produces a logical matrix
            % that needs to be cast to an uint16 before passing to svd
            mQ = ones( M, 1 );
            mQ( logical( A_zeros(:,i) ) ) = 0;
            mQ = diag( mQ );
            cU{i} = null( mQ );
        end;
        
    end;
    
    if ~isempty( F_zeros )

        for i = 1:M,
            mR = ones( K, 1 );
            mR( logical( F_zeros(:,i) ) ) = 0;
            mR = diag( mR );
            cV{i} = null( mR );
        end;

    end;
    
% end function
end