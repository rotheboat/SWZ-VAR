% Roland Meeks (2010-11)
%
% Given basic data and information about a BVAR, supplies prior parameter
% values for Minnesota or Normal-Wishart priors.
%
% Note that the expected form of the design matrix Ylag excludes any 
% constants that may be present in the VAR.
%
% Adapted from BVAR_ANALYT.m
%
% 2010-11-18 RM
function s = GetAnalyticPrior( prior_options, Y1, Ylag, T, varargin )


K = prior_options.K;
M = prior_options.M;
p = prior_options.p;

s = struct( 'A_prior', [], 'V_prior', [], ...
                      'v_prior', [], 'S_prior', []);
                  
% check if the optional arg giving diffs/levels prior is supplied
if ( 5 == nargin )
    % get the optional argument 'user difference vs level' prior
    usrdifflev = varargin{1};
    if ( M ~= length (usrdifflev) )
        warnmsg = [ 'GetAnalyticPrior::The vector `usrdifflev` is the ' ...
                 'wrong length. Using defaults.' ];
        warning( warnmsg );
        usrdifflev = 0;
    end;
else
    % supply the same default as in KK - that the VAR is specified in
    % differences
    usrdifflev = 0;
end;

%------------------Check if BVAR includes an intercept-----------------
constant = 0;

if ( K > M*p ) %( 1 == mod( K, M*p ) )
    constant = 1;
elseif ( K < M*p ) %( 1 < mod( K, M*p ) )
    % throw an exception
end;
%-----------------Prior hyperparameters for bvar model-----------------
% Define hyperparameters
% If Normal-Diffuse
% I guess there is nothing to specify in this case! Just return the set of
% empty matrices.
% Posteriors depend on OLS quantities
    
if (2 == prior_options.prior || 4 == prior_options.prior) % Minnesota or IN-Wishart
    % init the A matrix
    s.A_prior = zeros(K,M);
    
    % is there an intercept
    if ( 1 == constant )
        B = s.A_prior( 2:M+1, 1:M ); % select A_1 ...
        % ... then replace diagonal with optional levels/diffs opt vector
        B( linspace(1, numel(B), length(B)) ) = usrdifflev;
        s.A_prior( 2:M+1, 1:M ) = B;    
    else % there is not an intercept
        B = s.A_prior( 1:M, 1:M ); % select A_1 ...
        % ... then replace diagonal with optional levels/diffs opt vector
        B( linspace(1, numel(B), length(B)) ) = usrdifflev;
        s.A_prior( 1:M, 1:M ) = B;
    end;
    
        
    % Hyperparameters on the Minnesota variance of alpha
    % NOTE RM: in line with Banbura et al. parameterize in terms of lambda,
    % which is the square of the Koop-Korobilis a_bar_* parameters.
    a_bar_1 = prior_options.lambda;
    a_bar_2 = prior_options.lambda;
    a_bar_3 = 10^4; % RM Changed from 10^2, 2014-03-07

    % Now get residual variances of univariate p-lag autoregressions. Here
    % we just run the AR(p) model on each equation, ignoring the constant
    % and exogenous variables (if they have been specified for the original
    % VAR model). When we are doing h-step forecasts, these are going to be
    % regressions of x(t+h) on x(t), x(t-1),..., x(t-p).
    sigma_sq = zeros(M,1); % vector to store residual variances
    for i = 1:M
        % Dependent variable in i-th equation
 		Y_i = Y1(:,i);
        
        % Independent variable in i-th equation - select the relevant (own
        % lag) columns from the Ylag matrix
        Ylag_i = Ylag(:,(0:p-1)*M+i);
        
		% Because we only want to use in-sample data, drop the last obs.
		Ylag_i = Ylag_i(1:T,:);
		Y_i = Y_i(1:T,:);
        
        % OLS estimates of i-th equation [checked against PcGive 2010-11-18]
        alpha_i = (Ylag_i'*Ylag_i)\(Ylag_i'*Y_i);
        sigma_sq(i,1) = (1./(T-p))*(Y_i - Ylag_i*alpha_i)'*(Y_i - Ylag_i*alpha_i);
    end
    
    % Now define prior hyperparameters.
    % Create an array of dimensions K x M, which will contain the K diagonal
    % elements of the covariance matrix, in each of the M equations.
    V_i = zeros(K,M);
    
    % index in each equation which are the own lags
    ind = zeros(M,p);
    for i=1:M
        ind(i,:) = constant+i:M:K;
    end
    for i = 1:M  % for each i-th equation
        for j = 1:K   % for each j-th RHS variable
            if (1 == constant)
                if (1==j)
                    V_i(j,i) = a_bar_3*sigma_sq(i,1); % variance on constant                
                elseif find(j==ind(i,:))>0
                    V_i(j,i) = (a_bar_1^2)./(p); % variance on own lags (RM Changed to p from p^2)          
                else
                    for kj=1:M
                        if find(j==ind(kj,:))>0
                            ll = kj;                   
                        end
                    end
                    V_i(j,i) = ( (a_bar_2^2)*sigma_sq(i,1))./((p^2)*sigma_sq(ll,1));           
                end
            else
                if find(j==ind(i,:))>0
                    V_i(j,i) = (a_bar_1^2)./(p); % variance on own lags (RM Changed to p from p^2)
                else
                    for kj=1:M
                        if find(j==ind(kj,:))>0
                            ll = kj;
                        end                        
                    end
                    V_i(j,i) = ((a_bar_2^2)*sigma_sq(i,1))./((p^2)*sigma_sq(ll,1));            
                end
            end
        end
    end
    
    % Now V is a diagonal matrix with diagonal elements the V_i
    s.V_prior = diag(V_i(:));  % this is the prior variance of the vector a  
    
    if ( 2 == prior_options.prior )
        % Not relevant for pure MN prior (with fixed covariance)
        s.v_prior = 0;
        s.S_prior = 0;
    else
        % Independent normal-Wishart; see BVAR_GIBBS
        s.v_prior = M;
        % Changed [2013-02-18] from 1e-4*eye(M). The proposed covariance
        % matrix has the variances of the OLS residuals from a set of
        % autoregressions for each variable (see above)
        s.S_prior = 1e-3*diag( sigma_sq )/M; 
                                     % #BUGBUG Can one be sure this is 
                                     % small enough not to be influential
    end;
    
% ~~end MN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif (3 == prior_options.prior) % Normal-Wishart
    % Hyperparameters on a ~ N(a_prior, SIGMA x V_prior)

    % init the A matrix
    s.A_prior = zeros(K,M);
    
    % is there an intercept
    if ( 1 == constant )
        B = s.A_prior( 2:M+1, 1:M ); % select A_1 ...
        % ... then replace diagonal with optional levels/diffs opt vector
        B( linspace(1, numel(B), length(B)) ) = usrdifflev;
        s.A_prior( 2:M+1, 1:M ) = B;    
    else % there is not an intercept
        B = s.A_prior( 1:M, 1:M ); % select A_1 ...
        % ... then replace diagonal with optional levels/diffs opt vector
        B( linspace(1, numel(B), length(B)) ) = usrdifflev;
        s.A_prior( 1:M, 1:M ) = B;
    end;

    s.V_prior = 1e2*eye(K);
    % Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior)) - the prior S
    % = I and v = M ensure the distribution is non-singular...
    s.v_prior = M;
    s.S_prior = eye(M); % #BUGBUG This is a horrible default choice
    
% ~~end N-Wi~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end