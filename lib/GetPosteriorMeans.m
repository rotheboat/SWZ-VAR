% Roland Meeks (2010-11)
%
% Given posterior parameter estimates, returns posterior means. Require OLS
% estimates for the MN prior.
% @args prior = {1 | 2 | 3}
%       olsparams = struct
%       postest = struct
%       M = int
% @ returns struct( 'SIGMA', 'alpha_mean', 'alpha_var', 'v_post' )
%
% Roland Meeks 2010-11-23
function s = GetPosteriorMeans( prior, olsparams, postest, M )

% The reason for including v_post is to avoid having to include both
% postest and postmean data structures for methods that require both
% posterior degrees of freedom and posterior means.
s = struct( 'SIGMA', [], 'alpha_mean', [], 'alpha_var', [], 'v_post', [] );

if ( 1==prior )
    % mean of SIGMA (GCSR 574-5)
    s.SIGMA = (olsparams.SSE)./(postest.v_post - M - 1); 
    % Now get the mean and variance of the Multi-t marginal posterior 
    % of alpha
    s.alpha_mean = postest.A_post(:);
    s.alpha_var = ...
      (1/(postest.v_post - M - 1))*kron(postest.S_post,postest.V_post);
    s.v_post = postest.v_post;
    
elseif ( 2==prior )
    % In the MN case, SIGMA is fixed
    s.SIGMA = olsparams.SIGMA_OLS;
    % In the MN case, the mean is a_post and the variance is V_post
    s.alpha_mean = postest.A_post(:);
    s.alpha_var = postest.V_post;
        
elseif ( 3==prior )
    % mean of SIGMA (GCSR 574-5)
    % Get non-PDness due to numerical error unless you do this step...
    UT = chol(postest.S_post); % ensure S_post is symetric PD
    s.SIGMA = (UT'*UT)./(postest.v_post - M - 1);
    % Now get the mean and variance of the Multi-t marginal posterior 
    % of alpha
    s.alpha_mean = postest.A_post(:);
    s.alpha_var = ...
      (1/(postest.v_post - M - 1))*kron(postest.S_post,postest.V_post);
    s.v_post = postest.v_post;
  
end;