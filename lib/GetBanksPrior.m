% @author       Roland Meeks
% @date         2012-08-03
% @see          GetAnalyticPrior.m
% This function sets prior information for a subset of the variables in a
% Bayesian VAR. Priors for equations corresponding to other variables are
% assumed to have been set (the fields in the `priorparams' structure are 
% assumed not to be empty, e.g. GetAnalyticPrior may have been called).
%
% Suppose macro and banking variables were stacked as follows:
%   y = [y_m' y_b' ]'
% Then partition the AR matrix as
%
%          | A_mm    A_mb |
%   A =    |              |
%          | A_bm    A_bb |
%
% The prior information on AR coefficients imposed below is on the A_bm
% and A_bb matrices.
%
% This function deals with the possibility that the macro and banking
% variables are not in the order given above.
%
% The prior precision of A is given by V. This function allows the overall
% precision of the information in A_bm and A_bb to be set. When precision
% is zero, the prior is diffuse. For 'large' values of the precision
% parameter, the prior is very concentrated around the micro estimates.
% 
% The function is passed the following struct from GetAnalyticPrior():
%
% s = struct( 'A_prior', [], 'V_prior', [], ...
%                                      'v_prior', [], 'S_prior', [])
%
% panelprior is expected to be a struct with the fields:
%
% index_prior_parameters : tells you where each variable in the VAR
%           specification appears in the panel specification, so you can
%           map from A_prior into the VAR coefficient matrix
% A_prior : the estimated coefficient matrix, with values for A_bm and A_bb
%           and zeros elsewhere
% V_prior : the estimated precision matrix for A
% use_Abm_prior : if TRUE then estimated A_bm matrix (see above) is imposed
% use_V_prior : if TRUE then the estimated precision matrix for A is used
% rho : scale factor for inverse of prior precision
%
function [s, panelprior] = GetBanksPrior( panelprior, s, index, rho, prior, M, K, p )
    
%% A.0 Preliminaries
    % Check if BVAR includes an intercept
    if ( K > M*p )
        constant = 1;
    elseif ( K == M*p )
        constant = 0;
    end;
    
    % Ensure that index is a row vector
    [ ir, ic ] = size( index );
    if ( ir > ic )
        index = index';
    end;
    
    % Ensure that the panelprior coefficient matrix has right dimension
    if ~all( size( s.A_prior(1+constant:end,:) ) == size( panelprior.A_prior ) )
        error( 'GetBanksPrior:MatrixDimMismatch', ...
            ['Something is wrong: Perhaps the number of lags in the ' ...
             'panel estimates is different from the number of lags in ' ...
             'the VAR'] );
    end;

    % Prior AR matrices for vector y (ignore constant) 
    A = s.A_prior(1+constant:end,:) ;
    
    % Set dimensions of the banking block (b) and the macro block (m)
    n_b = numel( index );
    n_m = size( A, 2 ) - n_b;
    
    % reorder the rows & cols of panelprior coefficient matrices to match
    % the ordering in the BVAR as defined by the user in `names'.
    Q = eye( M, M );
    Q = Q( 1:M, panelprior.index_prior_parameters );
    panelprior.A_prior = kron( eye( p, p ), Q )*panelprior.A_prior*Q';
    % same for matrix of diagonal elements of covariance matrix of
    % parameter estimates from panel (if used)
    if panelprior.use_V_prior
        assert( ~isempty( panelprior.V_prior ), ...
                        'use_V_prior is true but V_prior is empty.' );
        panelprior.V_prior = kron( eye( p, p ), Q )*panelprior.V_prior*Q';
    end
    
%% A.1 Need to reorder the rows and columns of certain matrices so that they
    % corresponds to the block structure given above (pre and post multiply
    % by an orthogonal matrix). Variables are not reordered within the 
    % macro and banking  blocks. So an initial ordering of { m1, b1, m2, 
    % b2, b3, m3 } winds up as { m1, m2, m3, b1, b2, b3 }. 
    %
    % Make the vector of reordered indexes:
    v0 = 1:M;
    v0( index ) = [];
    v0 = [v0 index];
    % make the permutation matrix; note that *pre* multiplication permutes 
    % rows, *post* multiplication permutes columns.     
    P = eye( M, M );
    P = P( v0, 1:M ); % think of the arguments here as coordinates
        
%% A.2. Matrix A is (Mp x  M), so AP correctly exchanges columns, but to do 
    % the same for rows need to pre-multiply by kron( I_p, P ). The 
    % dimensions of the product are (pM x pM) x (pM x M) x (M x M ) = 
    % (pM x M) as with the initial A matrix and as required.
    A = kron( eye( p, p ), P )*A*P' ;
    panelprior.A_prior = kron( eye( p, p ), P )*panelprior.A_prior*P' ;
    
%% A.3. The expected form of banksprior is a matrix [ A_bm  A_bb ] x p. The
    % arrangements of elements in A_prior is equations in columns and
    % variables in rows. So column m corresponds to the equation for
    % variable m in the reordered VAR. In rows are variables x lags, so we
    % first get coefficients on each endogenous variable for lag 1, then
    % lag 2 and so on.
    if ( 0 == panelprior.use_Abm_prior )
        row_select = kron( ones(p,1), [(n_m+1):M]' ) + ...
                    M*kron( tril(ones(p,p),-1)*ones(p,1), ones(n_b,1) );
    % n.b. the length of row_select is n_b*p            
    else
        row_select = 1:M*p;
    end;
    
    % only select columns corresponding to the banking equations
    A( row_select, n_m+1:end ) = panelprior.A_prior( row_select, n_m+1:end );

    
%% A.4. Reorder the rows and columns of A to return it to its *original* 
    % form
    A = kron( eye( p, p ), P' )*A*P ;
 
%% A.5. Set the new value of A_prior in the priorparams structure, taking
    % account the presence/absence of a constant in the model
    s.A_prior(1+constant:end,:) = A;
    
 % --- end AR coefficients prior --- %
 
%% *** Set prior covariance matrix *** %%
 if ( 2 == prior || 4 == prior ) % Minnesota or Indep N-Wi
 
    V = s.V_prior;
    % this is a M*(constant+M*p) x M*(constant+M*p) dimensional object
    v1 = 1:M*K; % index vector for all coefficients
    
%% V.1. Where a constant is present, get the rows corresponding only to the 
    % autoregressive coefficients
    if ( 1 == constant )
        % drop the rows corresponding to the constant (1, 2+M*p, 3+2*M*p,...)
        const_select = tril( ones(M,M) )*ones(M,1) + ...
                            tril( ones(M,M), -1)*ones( M,1)*M*p;
        v1( const_select ) = []; % drop constants
        V = V( v1, v1 ) ;   % select only rows/cols
    end;        

    
%% V.2. Pick out coordinates corresponding to variables in the banking 
    % block of the vectorized coefficient matrix vec(A)
    bank_row_select = ...
        kron( ones(p,1), index' ) + ...
            M*kron( tril(ones(p,p),-1)*ones(p,1), ones(n_b,1) );
    % gives all the rows of A where a banking variable appears in a banking
    % equation, as specified in 'index' function argument
    index_select = kron( index'-1, ones( n_b*p,1) )*M*p + ...
        kron( ones(n_b,1), bank_row_select ) ; 
%         kron( tril(ones(M,M),-1)*ones(M,1), ones( n_b*p,1) )*M*p + ...
%             kron( ones(M,1), bank_row_select );
    % gives all the rows of vec(A) where a banking equation appears, as
    % specified in 'index' function argument
    
  
%% V.3. Impose the banksprior matrix on the submatrix of V that corresponds
    % to the row and column indexes of the banking block in vec(A)

    % if we are not to use prior information on variances from all the
    % coefficient in the banking equations, just use the restricted
    % coefficient list
    if ( 0 == panelprior.use_Abm_prior )
        % create a temporary vector which has non-zero elements in the
        % position of banking variables in banking equations
        D1 = zeros( M*p, 1 );
        D1( index_select ) = index_select ;        
        % create a temporary matrix which has D1 on its diagonal
        D2 = diag( D1 );
    else
        % use all the prior information on variances, by just adopting all
        % the parameters that actually get read in...
        
        % create a temporary vector which contains the estimated parameter
        % variances from the panel, stacked by equation
        D1 = vec( 0 ~= panelprior.V_prior );
        % create a temporary matrix with the estimated parameter variances
        % along its diagonal, zeros elsewhere
        D2 = diag( D1 );
    end;
    
    % if the flag use_V_prior is set, then replace the relevant elements of
    % the V matrix with the estimated variances in D1
    if ( panelprior.use_V_prior )
        V1 = diag( vec( panelprior.V_prior ) );
        V( 0 ~= D2 ) = V1( 0 ~= D2 );
    end;
    
    % case 1: no Abm prior
    % scale the prior variances of banking variables in banking equations
    % by a factor rho
    %
    % case 2: Abm prior
    % scale the prior variance of banking and macro variables in banking
    % equations by a factor rho
    V( 0 ~= D2 ) = rho*V( 0 ~= D2 ) ; 
    
    % Ensure the new V prior is not singular
    if ( rank( V ) ~= size( V, 1 ) )
        warning( 'GetBanksPrior:CovarianceMatrix', ...
            'The prior covariance matrix V may be singular.' );
    end;
    
    
%% V.4. Set the new value of V_prior in the priorparams structure, where we
    % again ignore rows corresponding to constant (where present) since we
    % have done nothing to those rows above.
    s.V_prior( v1, v1 ) = V;
    
end; % end if Minnesota prior

disp( 'Done setting micro prior' );

end % end mfile