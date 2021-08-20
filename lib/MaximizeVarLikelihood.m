function [params, maxlik, mA0, mA1, mSigma, exitflag, minoutput] = ...
                                    MaximizeVarLikelihood( Y, X, varargin )
% Gaussian maximum likelihood estimation of a multivariate time series
% model.  Computes MLEs for the structural VAR:
%           y(t)'A = x(t)'F + e(t)'
% Note that variables are in rows, equations in columns. By default, i.e. 
% when no options are specified, estimates subject to the normalization 
% that the A matrix is upper triangular (recursive form). Linear 
% restrictions of the form:
%   Q(i)*a(i) = 0
%   R(i)*f(i) = 0  ... for i = 1,...,M
% are permitted on the columns of A and F.  
%
% @author:  Roland Meeks
% @date:    2014-06-10
%
%
% The optimization routine minFunc is used. Default starting values for the
% maximization routine are the OLS estimates.
% - Optionally use prior parameters as starting values (makes sense in the
%   context of estimating a Bayesian VAR).
% 
%
% @args:
%   Required
%       Y = data matrix, dependent variables
%       X = data matrix, explanatory variables
%   Optional
%       theta = vector containing starting values for numerical
%           optimization
%       mRestrict = cell array.  
%           mRestrict{1} is a matrix of size(A) with zeros wherever a 
%               restriction is placed on A.
%           mRestrict{2} is a matrix of size(F) with zeros wherever a
%               restriction is placed on F
% @returns:
%
% @requires:
% minFunc package by Mark Schmidt, available at:
% https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html
% If using minFunc for the first time, be sure to run mexAll.m before
% running this code.
%
%

% Establish dimensions of the VAR
[ T1 M ] = size( Y );
[ T2 K ] = size( X );
theta = []; % starting parameter values

model = [];
model.var_type = 'reduced_form';
model.is_restricted = false;

% *** create 1-by-M cell arrays of restriction matrices cU, cV

% For mA, the default is a triangular normalization
cU = cell(1,M);
for i=1:M
    cU{i} = eye(M,i);
end

% For mF, the default is no restrictions
mycell = cell(1,M); 
cV = cellfun(@(n)eye(K,K),mycell,'UniformOutput', false);
mRestrict = cell(1,2);


% Inputs:
nargin = length( varargin );

switch nargin
    case 0 % We must use the defaults
        error( 'Right now the defaults are not setup #BUGBUG' );
    
    case 1 % We are passed some starting values or some restrictions
        input = varargin{1};
        if iscell( input )
            [cU cV] = processRestrictions( input, cU, cV );
            mRestrict = {cU, cV};
        else
            theta = checkParameterVector( input, M, K );
        end;

    case 2 % We are passed start values and a string indicating var_type
        
        input = varargin{2};
        
        % unrestricted reduced/structural form
        if ischar( input )
            model.var_type = input;
            theta = checkParameterVector( varargin{1}, M, K );            
            
        % restricted structural form (if restrictions not empty) else
        % unrestricted reduced form
        else
            mRestrict = input;
            assert( iscell(mRestrict), 'MaximizeVarLikelihood:InputFormat', ...
                'Second input must be char or cell.' );
            
            % check if restrictions are empty
            model.is_restricted = ...
                ~all(cellfun( @(n)isempty(n), mRestrict, 'UniformOutput', true ));
            
            % if there are some restrictions, process them
            if model.is_restricted
                [cU cV] = processRestrictions( input, cU, cV );
                mRestrict = {cU, cV};
                model.var_type = 'structural_form';
                % Because of the restrictions, 
                theta = varargin{1};
                
            % if the model is in fact unrestricted (empty restriction
            % matrices were passed) then we'll go with the default cU
            % and cV.  
            else
                warning( 'MaximizeVarLikelihood:RestrictionMatrices', ...
                            'Restriction matrices are empty.'  );
                theta = checkParameterVector( varargin{1}, M, K );                
            end;
        end;
        
    otherwise % We are passed stuff we don't need
        % ignore for now
end;

fprintf( ['\nStarting maximum likelihood estimation of ' model.var_type ' VAR' ] );

% *** Options for minFunc ***
%
%     - 'sd': Steepest Descent
%         (no previous information used, not recommended)
%     - 'csd': Cyclic Steepest Descent
%         (uses previous step length for a fixed length cycle)
%     - 'bb': Barzilai and Borwein Gradient
%         (uses only previous step)
%     - 'cg': Non-Linear Conjugate Gradient
%         (uses only previous step and a vector beta)
%     - 'scg': Scaled Non-Linear Conjugate Gradient
%         (uses previous step and a vector beta, 
%             and Hessian-vector products to initialize line search)
%     - 'pcg': Preconditionined Non-Linear Conjugate Gradient
%         (uses only previous step and a vector beta, preconditioned version)
%     - 'lbfgs': Quasi-Newton with Limited-Memory BFGS Updating
%         (default: uses a predetermined nunber of previous steps to form a 
%             low-rank Hessian approximation)
%     - 'newton0': Hessian-Free Newton
%         (numerically computes Hessian-Vector products)
%     - 'pnewton0': Preconditioned Hessian-Free Newton 
%         (numerically computes Hessian-Vector products, preconditioned
%         version)
%     - 'qnewton': Quasi-Newton Hessian approximation
%         (uses dense Hessian approximation)
options = [];
options.method = 'lbfgs';    % see algorithm options above
options.numDiff = 2;         % 'central differencing' for num. derivatives
options.MaxIter = 5e4;       % default 500
options.MaxFunEvals = 5e7;   % default 1000
options.Display = 'iter';
% note: have been using optTol = 1e-9 and progTol = 1e-12
options.optTol = 1e-6;       % default 1e-5 Term. tol. 1st-order optimality
options.progTol = 1e-9;     % default 1e-9 Term. tol. func/param chgs

% *** minFunc ***
if model.is_restricted
    funObj = @(params)(-1)*GetVarLogLikelihood( Y, X, params, model.var_type, mRestrict );
else
    funObj = @(params)(-1)*GetVarLogLikelihood( Y, X, params, model.var_type );
end;

[params lik exitflag minoutput] = minFunc( funObj, theta, options );

% Post-estimation output
fprintf( '\nIterations = %6.0f\n', minoutput.iterations );
fprintf( 'Function evals = %6.0f\n\n', minoutput.funcCount );

%
if model.is_restricted
    [ maxlik mA0 mA1 mSigma ] = GetVarLogLikelihood( Y, X, params, model.var_type, mRestrict );
else
    [ maxlik mA0 mA1 mSigma ] = GetVarLogLikelihood( Y, X, params, model.var_type );
end;


end % end function


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

% If a vector is passed, check that it is of the right dimensions for a
% reduced form VAR
function theta = checkParameterVector( input, M, K )

    theta = input;
    b_input_ok = isnumeric( theta ) && ...
                any( eq( 1, size( theta ) ) ) && ...
                    any( eq( M*K+M*(M+1)/2, size( theta ) ) );
    err_msg = 'Method requires vector size M*K+M*(M+1)/2';
    assert( b_input_ok, 'MaximizeVarLikelihood:ParamVector', err_msg );

end