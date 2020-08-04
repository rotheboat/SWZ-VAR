function [ priorparams, options ] = SetSimsZhaPrior( mY, mX, prior_name, varargin )
% Given the data (Y, X), an option determining the nature of the prior, and
% some optional arguments specifiying parameters of the prior, returns a
% structure containing the matrices describing the prior distributions of 
% the parameters in the A and F matrices of the VAR:
%   y(t)'A = x(t)'F + e(t)', e(t) ~ N(0,I)
% where x(t)' = [ 1', y(t-1)',...,y(t-p)' ].  Notation follows Waggoner and 
% Zha (JEDC 2003).
%
% Where restrictions are present, the code proceeds as
% 
%   GetNamedPrior -> SetRestrictions -> GetDerivedPriorUnderRestrictions
%
% 
%
% @author: Roland Meeks
%
% A vector of hyperparameters can be passed using options.lambda, with
% lambda a 6-by-1 vector of hyperparameters (the following numbers for 
% Atlanta Fed's monthly forecasting model)
%   lambda(1): overall tightness and also for A0;  (0.57)
%   lambda(2): relative tightness for A+;  (0.13)
%   lambda(3): tightness on lag decay;  (1)
%   lambda(4): relative tightness for exogenous variables;  (0.1)
%   lambda(5): weight on nvar sums of coeffs dummy observations (unit roots);  (5)
%   lambda(6): weight on single dummy initial observation including constant
%               (cointegration, unit roots, and stationarity);  (5)
%
% @update: June, 2020
% Amended the prior to deal with the case where there are exogenous
% variables (including deterministic variables). This is allowed for in the
% WZ set-up, but the literature does not provide specific guidance on how
% to set priors on such variables (aside from the constant term). I choose
% to treat all exogenous variables symmetrically, and recommend a diffuse
% prior. There may be N exogenous variables (including the constant) in the
% M equations of the VAR, with N typically a small number. Where a single
% exogenous variable appears at various lags, I do not apply shrinkage to
% the lag coefficients.
%
% @update: July, 2020
%   
%
    % init: create an empty options struct.  This object is returned by the
    % function
    options = struct([]);

    % If an optional argument has been passed, it is an options struct
    nargin = length( varargin );
    if ( 1 == nargin )
        options = varargin{1};
    end

    % init: create an empty priorparams struct.  This object is returned by the
    % function.
    priorparams = struct( ...
        'Hbar', [], 'Pbar', [], 'Sbar', [], 'Stld', [], ...
        'Stld_inv', [], 'Ptld', [], 'Htld_inv', [], ...
        'cV', [], 'cU', [], 'mYd', [], 'mXd', [] );

    % Get the parameter matrices in equation (6) of WZ, describing the
    % distribution of the unrestricted prior density
    priorparams = GetNamedPrior( mY, mX, prior_name, priorparams, options );
    
    % Set the normalizing or overidentifying restrictions, which are
    % summarized in the matrices inside the cU and cV structures
    priorparams = SetRestrictions( mY, mX, priorparams, options );
    
    % Take the prior without restrictions, and derive the prior with them.
    % Get the parameter matrices in equation (10) of WZ, describing the
    % distribution of the restricted prior density
    priorparams = GetDerivedPriorUnderRestrictions( priorparams );

end


%% GetNamedPrior
function priorparams = GetNamedPrior( mY, mX, prior_name, priorparams, options )
% This method returns the matrices of the unrestricted prior distribution
% as per the original Sims and Zha prior.  For the prior distribution under
% the restrictions summarized by Q_i and R_i, see method
% GetDerivedPriorUnderRestrictions(.)

% Given the data matrices mY and mX, derives and returns:
%   M : the dimension of the VAR
%   K : the number of RHS variables per equation (including exog/deter)
%   N : the number of exog RHS variables, incl. constant (0 if none)
%   p : the order of the VAR
%   constant : whether the VAR has a constant (1) or not (0)
[T, M, K, N, p, constant] = GetVARDimensions( mY, mX, options );

% Check for the existence of a vector of hyperparameters
if isfield( options, 'lambda' )
    lambda = options.lambda;
    
    if ( 0==lambda(4) && constant )
        error( 'Where a constant is present, its prior variance may not be set to zero.' );
    end
else            
    % Use defaults: On monthly data Tao Zha uses:
    %   lambda = [ 0.57, 0.13, 1, 0.1, 5, 5 ]; 
    % On quarterly data it's suggested to use fairly 'diffuse'
    % values:
    options.lambda = [ 1, 1, 0.1, 1, 1, 1, 0 ];
    lambda = options.lambda;
    % Brandt (UTD) suggests the following "tight" prior:
    %   lambda = [ 0.6, 0.1, 1, 0.1, 5, 5 ]
    % "loose" prior:
    %   lambda = [ 0.6, 0.15, 1, 0.15, 2, 2 ]
    % and "diffuse" prior:
    %   lambda = [ 1, 10, 0, 10, 0, 0 ]
end


% convert prior name to upper case
prior_name = upper( prior_name );

% deal with the possible cases
switch prior_name
    
    case 'LITTERMAN'
        % ---- prior P matrix ------------------
        % For the Litterman-style prior, the P matrix takes a particularly
        % simple form embodying the idea that the reduced form variables
        % follow a random walk...

        % For each equation, prior mean is E[ f(i)|a(i) ] = P(i)a(i).  This
        % is more general than E[ f|a ] = kron( I, P )a in that it allows
        % different equations to be different functions of a.  Naturally
        % it's not fully general since prior independence between equations
        % is maintained.  The baseline Sims-Zha interpretation of Litterman
        % has the same form for each equation however.
        
        % Note that if certain options have been passed, they will be
        % ignored where LITTERMAN is the prior_name, and a warning message
        % will appear
        custom_prior_warning_message = ...
            cat( 2, 'SetSimsZhaPrior::PriorOptions::To set prior ', ...
            'parameters to custom values,\nset prior_name to "CUSTOM". ', ...
            'Continuing with "LITTERMAN" prior.' );
        
        % RM 2020-06-08
        % Amended to make allowance for exogenous variables. I assume that
        % the prior for all exogenous variables (including deterministic
        % ones) is the same by default, as the literature provides no
        % guidance on this point. Thus e.g. M = 2 with a constant:
        %   mP = [ 1 0;
        %          0 1;
        %          0 0 ]; 
        mP = cat(1, eye( M*p, M ), zeros( N, M ) );
        
        mycell = cell(1,M);
        P_prior = cellfun( @(X)mP, mycell, 'UniformOutput', false );

        % ---- prior S matrix ------------------
        % Now get residual variances of univariate p-lag autoregressions. 
        % Here we just run the AR(p) model on each variable (no intercept)
        mS = zeros(M,1); % vector to store residual variances

        if isfield( options, 'sigma_sq' )
            warning( custom_prior_warning_message, '%s' );
        end

        for i = 1:M
            % Dependent variable in i-th equation
            Y_i = mY(:,i);
            %Y_i = Y_i - mean( Y_i );

            % Lagged variables in i-th equation
            % #BUGBUG define N as #exog vars
            Ylag_i = mX(:,i:M:K-N );
            %Ylag_i = Ylag_i - ones(T,1)*mean( Ylag_i );

            % OLS estimates of i-th equation
            alpha_i = (Ylag_i'*Ylag_i)\(Ylag_i'*Y_i);
            mS(i) = sumsqr(Y_i - Ylag_i*alpha_i)/(T-p);
        end % end for
        
        % SZ (p. 961): "The individual elements [of the non-zero parts of 
        % A] assumed independent, with prior standard deviations of all 
        % elements in the i'th row set to lambda(1)/hat(sigma)_i
        mS = lambda(1)^2*diag( 1./mS );

        % The prior is described by a set of M identical mV matrices.  The
        % implicit assumption is that the prior is unrestricted.  The derived
        % prior distribution under linear restrictions is a function of S_prior.
        mycell = cell(1,M);
        S_prior = cellfun( @(X)mS, mycell, 'UniformOutput', false );

        
        % ---- prior H matrix ------------------
        % Check if a prior variance matrix has been passed
        if isfield( options, 'mV' )
            warning( custom_prior_warning_message, '%s' );
        end
        
        % The default prior for H is very similar to that for S; see 
        % SZ (p. 955). 
        mH = repmat( lambda(2)^2*diag( mS ), p, 1 );
        mH = mH./power( kron( 1:p, ones( 1, M ) )', 2*lambda(3) );
        % all exogenous variables including the deterministic ones 
        % (e.g. constant term) have s.d. lambda(1)*lambda(4)
        if 0 < N
            % Exogenous variables, then constant, appear after lags of the
            % endogenous variables. Exogenous variables should be scaled by
            % their variances to avoid allowing arbitrary re-scalings of
            % data to effect the prior.
            
            % Get the 1 x N vector of variances of exogenous variables
            var_exog = var( mX( :, M*p+1:K ) )';
            
            if constant
            % Where a constant is present, set normalization to unity
                var_exog( N ) = 1;
            end
            
            mH = diag( [ mH; (1./var_exog)*(lambda(1)*lambda(4))^2 ] );  
        else
            % No exogenous variables, no constant
            mH = diag( mH );
        end
            
        mycell = cell(1,M);
        H_prior = cellfun( @(X)mH, mycell, 'UniformOutput', false );
    
    case 'CUSTOM'
        % The custom case is one where the user can supply the key input
        % parameters to construct the Sims-Zha prior. These include:
        %   options.mB
        %   options.mV
        %   options.sigma_sq
            
        % ---- prior P matrix ------------------
        % For the Litterman-style prior, the P matrix takes a particularly
        % simple form embodying the idea that the reduced form variables
        % follow a random walk...
        cP = cell(1,M);
        
        % Check whether we've been passed a reduced form parameter matrix
        % mB, which will stand in for mP below
        if isfield( options, 'mB' )
            mB = options.mB;
        end
        
        
        % For each equation, prior mean is E[ f(i)|a(i) ] = P(i)a(i).  This
        % is more general than E[ f|a ] = kron( I, P )a in that it allows
        % different equations to be different functions of a.  Naturally
        % it's not fully general since prior independence between equations
        % is maintained.
        for i = 1:M
            
            if isempty( mB )
                mP = eye( M*p, M );
            else
                mP = mB;
            end
            
            if constant
                mP = [ mP; ...
                       zeros(1,M) ];
            end
            
            cP{i} = mP;           
        end % end for 
        
        P_prior = cP;

        % ---- prior S matrix ------------------
        % Now get residual variances of univariate p-lag autoregressions. 
        % Here we just run the AR(p) model on each variable (no intercept)
        mS = zeros(M,1); % vector to store residual variances

        % check if the user has supplied a vector containing the residual
        % variances
        if isfield( options, 'sigma_sq' )
            assert( eq( M, size( options.sigma_sq, 1 ) ), ...
                'SetSimsZhaPrior:PriorOptions', ...
                    'Option `sigma_sq` must be a size(Y,2)-by-1 vector' );
            mS = options.sigma_sq;
        else
            for i = 1:M
                % Dependent variable in i-th equation
                Y_i = mY(:,i);
                %Y_i = Y_i - mean( Y_i );

                % Lagged variables in i-th equation
                % #BUGBUG define N
                Ylag_i = mX(:,i:M:K-N );
                %Ylag_i = Ylag_i - ones(T,1)*mean( Ylag_i );

                % OLS estimates of i-th equation residual variance
                alpha_i = (Ylag_i'*Ylag_i)\(Ylag_i'*Y_i);
                mS(i) = sumsqr(Y_i - Ylag_i*alpha_i)/(T-p);
            end
        end
        
        % SZ (p. 961): "The individual elements [of the non-zero parts of 
        % A] assumed independent, with prior standard deviations of all 
        % elements in the i'th row set to lambda(1)/hat(sigma)_i
        mS = diag( 1./mS );

        % The prior is described by a set of M identical mV matrices.  The
        % implicit assumption is that the prior is unrestricted.  The derived
        % prior distribution under linear restrictions is a function of S_prior.
        mycell = cell(1,M);
        for i=1:M
            mycell{i} = lambda(1)^2*mS;
            
            if isfield( options, 'bank_var_index' )
                if any( options.bank_var_index == i )
                    mycell{i} = options.bank_lambda(1)^2*mS; % Mult. by factor
                end
            end
        end
        S_prior = mycell;
        %cellfun( @(X)mS, mycell, 'UniformOutput', false );

        
        % ---- prior H matrix ------------------
        mycell = cell(1,M);
        
        % BUGBUG here I seem to be allowing a matrix mV to be passed, which
        % directly sets mH
        if isfield( options, 'mV' )
            if ~isempty( options.mV )
                % Process
                mV = diag( options.mV );
                for i=1:M
                    mH = mV;%( (i-1)*K+1:i*K );

                    % ****BUGBUG****
                    if false %any( 8:11 == i )
                        mH = lambda(1)^2*10*mH; % Mult. by factor
                    else
                        mH = lambda(1)^2*lambda(2)^2*mH;                
                    end
                    
                    if constant % use own settings for constant
                        % #BUGBUG changed order of constant
                        mH = [ mH(2:end); lambda(1)^2*lambda(4)^2 ];  
                    end
                    mycell{i} = diag( mH );
                end
            else
                error( 'SetSimsZhaPrior:PriorOptions', ...
                       'Prior option `mV` is an empty matrix' );
            end
        else
            for i = 1:M
                % The prior for H is very similar to that for S; SZ (p. 955)            
                mH = repmat( diag( mS ), p, 1 );
                mH = lambda(1)^2*lambda(2)^2*mH;
                    
                mH = mH./power( kron( 1:p, ones( 1, M ) )', 4*lambda(3) );
                % constant term has s.d. lambda(1)*lambda(4)
                if constant
                    mH = diag( [ lambda(1)^2*lambda(4)^2; mH ] );  
                else
                    mH = diag( mH );
                end

                mycell{i} = mH;
            end % end for
        end % end if

        H_prior = mycell;% cellfun( @(X)mH, mycell, 'UniformOutput', false );
                
    otherwise
end

priorparams.Hbar = H_prior; 
priorparams.Pbar = P_prior;
priorparams.Sbar = S_prior;

%% ---- dummy variable priors ----------
if isfield( options, 'dummy_data' )
    priorparams.mYd = options.dummy_data{1};
    priorparams.mXd = options.dummy_data{2};   
end

% Compute the average initial condition
temp = mX(1, 1:K-N);
Ymean = mean( reshape( temp, M, p ), 2 )';
%         Ymean = mean( mY, 1 );
%         is_stationary = zeros(1,6); %% ***** BUGBUG *****
%         is_stationary( [4:8] ) = 1; %% ***** BUGBUG *****
%         Ymean( ~logical( is_stationary ) ) = Ymean( ~logical( is_stationary ) )/1e3;  %% ***** BUGBUG *****
        
%% ---- dummy inital observation ----------
mYd = NaN(1,M);
mXd = NaN(1,K);

mYd = lambda(6)*Ymean;

% #BUGBUG
mXd(1, M*p+1:K ) = lambda(6);
mXd(1, 1:K-N) = lambda(6)*repmat( Ymean, 1, p );

% if hyperparameter lambda(6) = 0, don't add dummy obs
if 0 < lambda(6)
    priorparams.mYd = mYd;
    priorparams.mXd = mXd;
end
        
%% ---- sum of coefficients ----------        
mYd = NaN(M,M);
mXd = NaN(M,K);

%         Ymean = mean( reshape( temp, M, p ), 2 )'; %% ***** BUGBUG *****
%         Ymean( logical( is_stationary ) ) = Ymean( logical( is_stationary ) )/100; %% ***** BUGBUG *****
        
mYd = lambda(5)*diag( Ymean );
mXd = cat( 2, repmat( mYd, 1, p ), zeros(M,N) );
        
% if hyperparameter lambda(5) = 0, don't add dummy obs        
if 0 < lambda(5)
    priorparams.mYd = [ priorparams.mYd; mYd ];
    priorparams.mXd = [ priorparams.mXd; mXd ];
end
        
%% ---- M unit roots with constant ----------
mYd = NaN(M,M);
mXd = NaN(M,K);

for i = 1:M
    Ytemp = zeros(1, M);
    Ytemp(i) = Ymean(i);
    mYd(i,:) = lambda(7)*Ytemp;
    mXd(i, M*p+1:K ) = lambda(7);
    mXd(i, 1:M*p) = lambda(7)*repmat( Ytemp, 1, p );
end
        
%         % ***** BUGBUG *****
%         Ytemp = zeros(1, M); % ***** BUGBUG *****
%         Ytemp(3) = Ymean(3); % ***** BUGBUG *****
%         mYd(3,:) = 10*lambda(7)*Ytemp; % ***** BUGBUG *****
%         mXd(3, logical(constant) ) = 10*lambda(7); % ***** BUGBUG *****
%         mXd(3, constant+1:end) = 10*lambda(7)*repmat( Ytemp, 1, p );
%         % ***** BUGBUG *****
        
        
% if hyperparameter lambda(6) = 0, don't add dummy obs
if 0 < lambda(7)
    priorparams.mYd = [ priorparams.mYd; mYd ];
    priorparams.mXd = [ priorparams.mXd; mXd ];
end        


    
end



%% GetDerivedPriorUnderRestrictions
function priorparams = GetDerivedPriorUnderRestrictions( priorparams )
% This method computes parameters of the derived prior under the 
% restrictions.  Cell arrays of 'tilde' (tld) matrices, WZ (10).  Note that 
% the dimensions of the tilde matrices are not the same as of the prior 
% `bar' matrices above; their dimensions depend on the number of free 
% parameters in equation i; this is always greater than 0.

    cV = priorparams.cV;
    cU = priorparams.cU;
    cSbar = priorparams.Sbar;
    cPbar = priorparams.Pbar;
    cHbar = priorparams.Hbar;

    cHtld_inv = cellfun( @(V,Hbar) V'*(Hbar\V), ...
                       cV, cHbar, ...
                       'UniformOutput', false );

    cPtld = cellfun( @(V,Htld_inv,Hbar,P,U) Htld_inv\(V'/Hbar)*P*U, ...
                       cV, cHtld_inv, cHbar, cPbar, cU, ...
                       'UniformOutput', false );

    % Just compute the inverse, since this is all that is required
    cStld_inv = cellfun( @(U,S,P,Ptld,H,Htld_inv) ...
                        U'*(S\U) + U'*P'*(H\P)*U - Ptld'*Htld_inv*Ptld, ...
                        cU, cSbar, cPbar, cPtld, cHbar, cHtld_inv, ...
                        'UniformOutput', false );
                    
    cStld = cellfun( @(Stld_inv) Stld_inv\eye(size(Stld_inv)), ...
                    cStld_inv, 'UniformOutput', false );
                    
     
    % Assign parameters for the derived prior for a to priorparams struct:
    % To avoid numerical errors making the covariance/inverse covariance
    % non-symmetric or non-PD, use Cholesky factorization.
    cCholStld_inv = GetCellArrayCholesky( cStld_inv );
    priorparams.Stld_inv = cellfun( @(S) S'*S, cCholStld_inv, ...
        'UniformOutput', false );
    
    cCholStld = GetCellArrayCholesky( cStld );
    priorparams.Stld = cellfun( @(S) S'*S, cCholStld, ...
        'UniformOutput', false );
    
    % Assign parameters for the derived prior for b and g to priorparams 
    % struct:
    priorparams.Ptld = cPtld;
    priorparams.Htld_inv = cHtld_inv;
end


%% SetRestrictions
function priorparams = SetRestrictions( mY, mX, priorparams, options )
% Sets restriction matrices U and V, either to "default" values or
% according to zero restrictions passed in the options.
% The default "restrictions" are for A to be normalized to an upper
% triangular matrix, and F to be unrestricted.  However, any restriction
% matrix can be passed. The call to processZeroRestrictions() is triggered
% by the presence of the options.restrictions takes the 
% basic input in mRestrict and converts it to U and  
% V matrices that have the property a(i) = U(i)b(i) and f(i)=V(i)g(i), 
% where b(i) and g(i) are the restricted parameter vectors.

% My standard approach is that mRestrict is a 1x2 cell array.

    [~, M, K, ~, ~, ~] = GetVARDimensions( mY, mX, options );

    % Default "restrictions": 
    % * A is normalized to be triangular; F is
    %   unrestricted.
    cU = cell(1,M);
    for i=1:M
        cU{i} = eye(M,i);
    end
    
    % * For mF, the default is no restrictions
    mycell = cell(1,M); 
    cV = cellfun(@(n)eye(K,K),mycell,'UniformOutput', false);

    % Check if any restrictions are present (if not the defaults are used).
    % There can never be "no restrictions".
    if isfield( options, 'restrictions' )
        % ** Zero restrictions
        % if there is a 'zero_restrictions' field, it should contain a cell
        % array
        if isfield( options.restrictions, 'zero_restrictions' )

            % check if the cell array of restrictions are completely empty.
            % It is expected behavior if one is empty and the other is not
            % (so restrictions apply only contemporaneously, or only at
            % lags and/or on exogenous variables).
            priorparams.is_restricted = ...
                        ~all( ...
                            cellfun( @(n)isempty(n), ...
                                options.restrictions.zero_restrictions, ...
                                    'UniformOutput', true ) ...
                            );

            % If there are some zero restrictions, process them
            if priorparams.is_restricted
                
                [cU, cV] = processZeroRestrictions( ...
                            options.restrictions.zero_restrictions, ...
                                M, K, cU, cV );
                
                % Signal that restrictions are present
                fprintf( '\nApplying parameter restrictions...' );

            end
        end

        % ** General restrictions
        % Where general linear restrictions are to be imposed, these are
        % already passed in the options struct:
        %   options.restrictions.linear_restrictions
        % However they must be validated here
        if isfield( options.restrictions, 'linear_restrictions' )
            
            % For each equation, retrieve the Q matrix from the options
            % structure, and store the nullity of Q in U. We expect Q to be
            % a 3-D numerical array, with the third dimension corresponing
            % to equation numbers with ordering exactly as in
            % options.names; same comments for R -> V.
            for i=1:M
                cU{i} = ...
                    null( options.restrictions.linear_restrictions.Q(:,:,i) );
                cV{i} = ...
                    null( options.restrictions.linear_restrictions.R(:,:,i) );
            end
        end
    end

    priorparams.cU = cU;
    priorparams.cV = cV;
    
end



function [cU, cV] = processZeroRestrictions( mRestrict, M, K, cU, cV )
% This method processes the restrictions summarized in mRestrict, a cell
% array of restriction matrices.  At present, only zero restrictions are
% handled.  There is a simple pattern-matching approach whereby:
% 	* mRestrict{1} is a M-by-M matrix with ones where mA has a non-zero 
%           coefficient; zeros elsewhere
% 	* mRestrict{2} is a K-by-M matrix with ones where mF has a non-zero 
%           coefficient; zeros elsewhere
% Either can be empty, meaning that the defaults defined above are used 
% (i.e. a triangular normalization on A, but no restrictions on the lags in 
% any equation). 
    if ~isempty( mRestrict{1} )

        A_zeros = mRestrict{1};
        assert( all( eq( [M,M], size( A_zeros ) ) ), ...
            'PriorRestrictions::Size of A restriction matrix invalid' );
        
        cU = cell( 1, M );

        for i = 1:M
            % The more obvious way to do create Q produces a logical matrix
            % that needs to be cast to an uint16 before passing to svd
            mQ = ones( M, 1 );
            mQ( logical( A_zeros(:,i) ) ) = 0;
            mQ = diag( mQ );
            cU{i} = null( mQ );
        end

    end

    if ~isempty( mRestrict{2} )

        F_zeros = mRestrict{2};
        assert( all( eq( [K,M], size( F_zeros ) ) ), ...
            'PriorRestrictions::Size of F restriction matrix invalid' );
        cV = cell( 1, M );

        for i = 1:M
            mR = ones( K, 1 );
            mR( logical( F_zeros(:,i) ) ) = 0;
            mR = diag( mR );
            cV{i} = null( mR );
        end

    end
    
% end function
end


function cCholFac = GetCellArrayCholesky( pdMat )
% Given a cell array of positive definite matrices, return a cell array of
% their Cholesky factors.  Returns the upper triangular matrix R such that
% R'*R = A.
    cCholFac = cellfun( @(S) chol(S), pdMat, 'UniformOutput', false );
end