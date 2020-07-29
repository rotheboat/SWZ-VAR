% Sample from the posterior predictive distribution of Y using simulation
% methods.
% 
% This function is required whenever we have analytic posteriors for the
% parameters, but lack an analytic posterior predictive for Y (as in the
% case of any multi-step forecast). It assumes we're interested in a 
% particular forecast horizon, rather than the whole forecast path, i.e in
% the distribution of Y(T+h) not the distribution of the entire path
% Y(T+1),...,Y(T+h). Thus the matrix returned has one row per simulation,
% with variables in columns.
%
% Roland Meeks 2010-11-24
function Ypred = DoForecastSimulation( Y, X, prior, postest, ndraws, h )

M = size(Y, 2); % dim VAR
K = size(X, 2); % number of params per equation
p = floor( size(X, 2)/M ); % num lags
constant = [];
if ( K > M*p )
    constant = 1;
end;
% The object of the following is to produce a posterior predictive
% distribution for the variables Y. To accomplish this, we simulate from
% the posterior of the parameters, then simulate Y forward given the
% parameters and initial Y(T). The resulting distribution of Y(T+h) is the
% required predictive distribution, i.e. p[ Y(T+h) | Y(T) ] having
% integrated out the parameters.

% Set the initial value of Y, i.e. X(T) including appropriate deterministic
% terms and lags of Y

% Set storage in memory for Ypred. It is dimension ndraws x dimvar: each
% posterior draw of Y(T+h) is in a single row, and variables are in the
% columns (i.e. same orientation as Y itself).
Ypred = zeros( ndraws, M );

% Loop over draws
for draw = 1:ndraws,
    % Because of the slightly dumb-ass way I'm just gonna iterate to form
    % forecasts below, I'll append the Y and X matrices with zeros. And while
    % I'm at it, I'll just truncate them so position 1 is the forecast origin,
    % position 2 is T+1, and so forth.
    Yf = [ Y(end,:); zeros(h, M) ];
    Xf = [ X(end,:); zeros(h, K) ];
    Vf = zeros( h+1, M );
    
    % Draw parameter matrices from posterior distribution. First draw the
    % covariance/scale param. Note that for the MN prior, it is fixed.
    % Conditional on Sigma, draw alpha. We then have a valid draw from the
    % joint posterior of all the parameters.
    
    % normal-diffuse
    if ( 1 == prior )
        % random draw from Wishart distribution
        SIGMA = wish( postest.S_post, postest.v_post );
        % the variance of alpha
        V_A = kron(SIGMA, postest.V_post);

    % MN
    elseif( 2 == prior )
        % variances are fixed
        SIGMA = postest.S_post;
        V_A = postest.V_post;
        
    % natural conjugate
    elseif( 3 == prior )
        % BUGBUG is this the correct way to draw from the inverse-Wishart?
        % I'm using the inverse scale in the regular wish() function. I'm
        % inverting the whole thing since its Sigma^-1 that is ~IW.
        SIGMA = inv( wish( inv(postest.S_post), postest.v_post ) );
        % the variance of alpha
        V_A = kron( SIGMA, postest.V_post );
        
    else
        error( 'Option `prior` must be one of {1,2,3}' );
    end;

    % Draw a (1xK) vector of params
    alpha = mvnrnd( postest.A_post(:)', V_A, 1 );
    % Make alpha conformable with Y
    A = reshape( alpha, K, M );
    % draw a (hxM) matrix of errors
    Vf( 2:end, : ) = mvnrnd( zeros(1, M), SIGMA, h );
    
    % Now we have a parameter draw and some errors, we can generate a path
    % for Y(s), s = T+1,...,T+h. As usual, keep time in the rows, and
    % variables in the columns. This piece shd be quick, as it's called
    % ndraws times, which could be lots.
    for j = 1:h,
        
        if ( p > 1)
            Xf(j+1,:) = [ constant Yf(j,:), Xf(j,end-p*M+1:end-M) ];
        else
            Xf(j+1,:) = [ constant Yf(j,:) ];
        end;
        
        Yf(j+1,:) = Xf(j+1,:)*A + Vf(j+1,:);
    end;

    % Store Y(T+h) in the Ypred matrix
    Ypred(draw,:) = Yf(h+1,:);
end;