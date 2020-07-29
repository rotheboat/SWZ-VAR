% Roland Meeks (2010-11)
%
% Returns posterior parameter estimates, given prior information, OLS
% quantities, and information on various VAR dimensions. Returns empty
% matrices if elements in (A, V, v, S) are not relevant for the chosen 
% prior.
%
% NB a generic problem of numerical accuracy arises when using operations 
% such as inv(), resulting in non-PD covariance matrices. I solve this by
% taking lower-triangular Choleski factors of the thing with error, then
% multiplying them up. The difference seems to be O(1e-14).
%
% BUGBUG it may be better to simply return the LT Cholesky factor as often
% this is required ultimately, rather than the full matrix
function s = GetAnalyticPosterior( strpr, strols, prior, p, post_options )

% Data & important constants of dimension
Y = strols.Y;
X = strols.X;
[ T M ] = size( Y );
K = size( X, 2 );

% Check if BVAR includes an intercept
if ( K > M*p )
    constant = 1;
elseif ( K == M*p )
    constant = 0;
end;

%--------- Output of function ---------------------------------------------
s = struct( 'A_post', [], 'V_post', [], 'v_post', [], 'S_post', [], ...
    'IRF', [] );

%--------- Posterior hyperparameters of ALPHA and SIGMA with Diffuse Prior
if prior == 1
    % Posterior of alpha|Data ~ Multi-T(kron(SSE,inv(X'X)),alpha_OLS,T-K)
    s.V_post = inv(strols.XpX);
    % Need to guarantee positive definiteness
    B = chol( s.V_post, 'lower' );
    s.V_post = B*B';
    a_post = strols.A_OLS(:);
    s.A_post = reshape(a_post,K,M);
    
    % posterior of SIGMA|Data ~ inv-Wishart(SSE,T-K)
    s.S_post = strols.SSE;
    % Need to guarantee positive definiteness
    B = chol( s.S_post, 'lower' );
    s.S_post = B*B';
    s.v_post = T-K;
        
%--------- Posterior hyperparameters of ALPHA and SIGMA with Minnesota Prior
elseif prior == 2
    % Retrieve coefficients from prior and OLS structs
    A_prior = strpr.A_prior;
    V_prior = strpr.V_prior;
    XpX = strols.XpX;
    A_OLS = strols.A_OLS;
    SIGMA_OLS = strols.SIGMA_OLS;
    
    % Note SIGMA is equal to the OLS quantity (and is fixed), but it will
    % be useful to set S_post equal to SIGMA on the understanding that it
    % represents a different quantity from the other priors.
    SIGMA = chol( SIGMA_OLS, 'lower' );
    SIGMA = SIGMA*SIGMA'; % guarantee positive definiteness
    s.S_post = SIGMA; % caution on intepretation
    
    % posterior of alpha|Data ~ N( a_post, V_post )
    s.V_post = inv( inv(V_prior) + kron(inv(SIGMA),XpX) );
    
    s.a_post = s.V_post * ( V_prior\A_prior(:) + ...
                        kron(inv(SIGMA),XpX)*A_OLS(:) );
    s.A_post = reshape(s.a_post,K,M);
    
    %=========IMPULSE RESPONSES:
    drawcount = 1;
    rejectcount = 0;
    flag = 1;
    
    imp = zeros( post_options.draws, M, M, post_options.irf_horizon );
    
    while ( drawcount < post_options.draws + 1 )
        
        if ( 0==(mod(drawcount,post_options.it_print)) && (1==flag) )
            disp( [ 'Iteration = ' num2str( drawcount ) ]);
            toc;
        end
    
        a_draw = s.a_post + chol(s.V_post)'*randn(K*M,1); % Draw of alpha
        Btemp = reshape(a_draw,K,M);
        Btemp = Btemp(constant+1:end,:); % changed from Btemp(2:end,:)
        % Form the companion matrix
        COMP = [ Btemp'; eye( M*(p-1) ), zeros( M*(p-1), M ) ];
       % Check the roots of the companion matrix, noting draws that violate
        % stationarity condition
        if ( (post_options.reject_nonstationary) && ...
                                ( any( abs( eig( COMP ) ) > 1.001 ) ) )
            flag = 0;
            rejectcount = rejectcount+1;
        else
            flag = 1;
            imp_resp = impulse( reshape( Btemp', M, M, p ), ...
                                     chol( SIGMA ), post_options.irf_horizon );
            imp(drawcount,:,:,:) = imp_resp ;
        end;        
        
        drawcount = drawcount + flag;
    end; % end while

    disp( 'Done with computing IRFs.' );
    disp( ['Rejected ' num2str( 100*rejectcount/(rejectcount+post_options.draws) ) ...
                                        '% of draws as non-stationary.'] )

    s.IRF = imp;

    %--------- Posterior hyperparameters of ALPHA and SIGMA with Normal-Wishart Prior
elseif prior == 3
    % Retrieve coefficients from prior and OLS structs
    A_prior = strpr.A_prior;
    V_prior = strpr.V_prior;
    v_prior = strpr.v_prior;
    S_prior = strpr.S_prior;
    XpX = strols.XpX;
    A_OLS = strols.A_OLS;
    SSE = strols.SSE;
    
    % posterior for alpha|Sigma,Data ~ N( a_post, Sigma**V_post )
    s.V_post = inv( inv(V_prior) + XpX );
    % Need to guarantee positive definiteness
    B = chol( s.V_post, 'lower' );
    s.V_post = B*B';

    s.A_post = s.V_post * ( V_prior\A_prior + (XpX)*A_OLS );
    
    % posterior for inv(Sigma)|Data ~ Wishart( v_post, inv(S_post) )
    s.S_post = SSE + S_prior + A_OLS'*(XpX)*A_OLS + ...
                    A_prior'*(V_prior\A_prior) - ...
                            s.A_post'*( inv(V_prior) + XpX )*s.A_post;
    % need to guarantee positive definiteness
    B = chol(s.S_post, 'lower' );
    s.S_post = B*B';
    s.v_post = T + v_prior;
end
