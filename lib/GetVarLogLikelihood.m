function [ lik mA0 mA1 mSigma ] = GetVarLogLikelihood( mY, mX, varargin )
% Computes the log likelihood for any Gaussian VAR in reduced or
% structural form, with or without linear restrictions on the parameters. 
% Returns the natural logarithm; divide the result by log( a ) for 
% base( a ) logarithm.
%
% @author:  Roland Meeks
% @date:  2012-12-18
%         2014-06-05 Modified to allow for a parameter *vector* to be
%         passed directly
%
% The function can handle four cases:
%   1. Unrestricted reduced form (default)
%   2. Unrestricted structural form (triangular normalization)
%   3. Restricted reduced form (probably not frequently used)
%   4. Restricted structural form
%
% The general form of the VAR in SF/RF is:
% 
%   SF: y(t)'A = x(t)'F + e(t)',  e(t) ~ N( 0, I )
%   RF:  y(t)' = x(t)'B + u(t)',  u(t) ~ N( 0, Sigma )
%
% where mSigma = inv( mA*mA' ) and mB = mF*inv( mA ).  Knowing which form
% we have avoids repeated matrix inversions when optimizing.  If a
% structural form model is specified without restrictions, a triangular
% normalization is assumed.
% 
% Restrictions on A and F are summarized by the conditions:
%   a(i) = Q(i)*b(i)
%   f(i) = R(i)*g(i)
% In reduced form, restrictions are allowed only on B and take a similar
% form.
%
% Input formats
% =============
% Multiple input formats are required.  
%
% ** One input
%   1. vB - assumes unrestricted reduced form by default
%   2. Cell array only:
%        SF: { mA, mF, I }          unrestricted/restricted
%        RF: { I, mB, mSigma }      unrestricted/restricted
%
% ** Two inputs
%   1. vB, char   - SF/RF     unrestricted
%
% ** Three inputs
%   1. vB, char, cell - SF    restricted
%      The cell contains two sets of arrays - the 'already processed'
%      restrictions cU{i} and cV{i} on mA and mF resp.
%
%   Note: don't support restricted reduced form for now
% 
% Optimizers such as minFunc pass a parameter vector, as in case 1.
% Although by default the code maps theta->{mB,mSigma}, we also need to
% handle the case where theta->{mA,mF,I}.  Optionally, can pass a string
% 'reduced_form' or 'structural_form' as argument 4.
%
%
%
% @args:
%   mY = a T-by-M matrix of observables
%   mX = a T-by-K matrix of lagged observables and deterministic terms
%
% @varargin:
%      
%
% @returns:
%   lik = scalar natural logarithm of the likelihood function evaluated at
%       the observations mY, conditional on mX and parameters mB, mSigma
%   mA0 = matrix of contemporaneous coefficients
%   mA1 = matrix of lag coefficient
%   mSigma = error covariance matrix 
%
% @see: GetVarLogLikelihood.tex for notation etc.

% Establish dimensions of the VAR
[ T1 M ] = size( mY );
[ T2 K ] = size( mX );

% Test for data matrix error
assert( eq( T1,T2 ), 'GetVarLogLikelihood:DataMatrix', ...
            'Y and X have an unequal number of rows' );

% init --------------------------------------------------------------------
nargin = length( varargin );
model = [];
model.var_type = 'reduced_form';
model.is_restricted = false;
model.input_type = 'vector';
mA = eye(M,M);
mB = [];
mF = [];
mSigma = eye(M,M);
vP = [];
cU = cell( 1, M );
cV = cell( 1, M );


%% check inputs -----------------------------------------------------------
try
    switch nargin 
        
        case 1
            input = varargin{1};

            % unrestricted reduced form
            if isnumeric( input ) % parameters, must be as a vector
                vP = input;
                assert( any(eq(1,size(vP))), ...
                    'GetVarLogLikelihood:InputFormat1', '' );
                [mB, mSigma] = mapParameterVectorToReducedFormMatrices( M, K, vP );
                
            % unrestricted structural/reduced form
            else 
                assert( iscell( input ), ...
                    'GetVarLogLikelihood:InputFormat1', '' ) % cell array of parameters
                mSigma = input{3};
                if eye(M,M) == mSigma, % structural form
                    mA = input{1};
                    mF = input{2};
                    model.var_type = 'structural_form';
                else                   % reduced form (default)
                    mB = input{2};
                    mSigma = input{3};
                end;
            end;
            
        case 2
            input = varargin{2};
            assert( ischar( input ), 'GetVarLogLikelihood:InputFormat2', '' );

            model.var_type = input;
            vP = varargin{1};
            assert( isnumeric( vP ), 'GetVarLogLikelihood:InputFormat2', '' );
            assert( any(eq(1, size(vP))), 'GetVarLogLikelihood:InputFormat2', '' );
                
            if strcmp( 'reduced_form', model.var_type )
                [mB, mSigma] = mapParameterVectorToReducedFormMatrices( M, K, vP );
            else
                [mA, mF] = mapParameterVectorToStructuralFormMatrices( M, K, vP );
            end;
            
        case 3
            % restricted structural/reduced form
            vP = varargin{1};
            model.var_type = varargin{2};
            cRestrict = varargin{3};
            
            b_inputs_ok = isnumeric( vP ) && any(eq(1, size(vP))) && ...
                ischar( model.var_type ) && iscell( cRestrict );
            assert( b_inputs_ok, 'GetVarLogLikelihood:InputFormat3', '' );

            model.is_restricted = ~all(cellfun( @(n)isempty(n), cRestrict, 'UniformOutput', true ));
            
        otherwise
            throw( MException( 'GetVarLogLikelihood:InputArguments', ...
                '' ) );
    end;
catch err% the inputs are wrong
    switch err.identifier
        case 'MATLAB:badsubscript'
            err_msg = ['Cell array must be in the format: ' ...
                       '{I,mB,mSigma} or {mA,mF,I}.' ];
                   
        case 'GetVarLogLikelihood:InputFormat1'
            err_msg = ['Call with a single input must supply either a ' ...
                       'parameter vector or a cell array with the format: ' ...
                       '{I,mB,mSigma} or {mA,mF,I}.' ];
                   
        case 'GetVarLogLikelihood:InputFormat2'
            err_msg = ['Call with a two inputs must supply either a ' ...
                       'parameter vector with a string `reduced_form`/'...
                       '`structural form`, or a pair of cell arrays'...
                       ' containing parameters and restrictions' ];
                   
        case 'GetVarLogLikelihood:InputFormat3'
            err_msg = ['Call with three inputs must supply a' ...
                       'parameter vector with a string `reduced_form`/'...
                       '`structural form`, and a cell array containing'...
                       ' restrictions'];

        case 'GetVarLogLikelihood:InputArguments'
            if 0 == nargin
                err_msg = 'Too few input arguments';
            else
                err_msg = 'Too many input arguments';
            end;
            
        case 'GetVarLogLikelihood:ParamVector'
            err_msg = ['The parameter vector is not conformable with ' ... 
                'the data matrix given model option: ' model.var_type];
        otherwise
            err_msg = 'Unknown error';
    end;
        
    error( 'GetVarLogLikelihood:InputFormat', err_msg );
end; % end try-catch        


%% end of input format checks ---------------------------------------------


% check if any restrictions have been passed ------------------------------
if model.is_restricted
    
    cU = cRestrict{1};
    cV = cRestrict{2};

    try
        [mA mF] = mapRestrictedParameterVectorToCoefficientMatrices( vP, cU, cV );
    catch err
        if strcmp( err.identifier, 'GetVarLogLikelihood:ParamVector' )
            err_msg = [ 'The parameter vector passed is not compatible' ...
                ' with the restrictions on the model as specified. ' ...
                'Check restrictions matrices?' ];
            error( 'GetVarLogLikelihood:ParamVector', err_msg );
        else
            rethrow(err);
        end
        
    end
end


% -------------------------------------------------------------------------
% Invoke library file mvnpdf; first argument is the data (variables in
% rows, observations in columns); second argument is the conditional mean
% at each observation date; third argument is conditional variance. If we
% have a structural VAR, adjust p(z) by the Jacobian of the transformation 
% z(t)' = y(t)'A to obtain p(y).
switch model.var_type
    case 'structural_form'
        p_Y = abs( det( mA ) )*mvnpdf( mY*mA, mX*mF, eye(M, M) );
        mA0 = mA;
        mA1 = mF;
    case 'reduced_form'
        p_Y = mvnpdf( mY, mX*mB, mSigma );
        mA0 = mA;
        mA1 = mB;
    otherwise
        error( 'GetVarLogLikelihood:ModelType', ...
            'model_type must be either `structural_form` or `reduced_form`' );
end;

% Obtain log likelihood
lik = sum( ln( p_Y ) );
% -------------------------------------------------------------------------

% end function
end


% Subfunction to allow for a parameter vector of length M*M + M*K to be
% passed to the GetVarLogLikelihood function, instead of two parameter
% matrices.  This is the form which, for example, a numerical optimizer
% will typically require.
function [mB, mSigma] = mapParameterVectorToReducedFormMatrices( M, K, theta )
    
    assert( eq( length(theta), M*(K+(M+1)/2) ), ...
                        'GetVarLogLikelihood:ParamVector', '' );
    
    % autoregressive coefficient matrix
    mB = theta( 1:M*K );
    mB = reshape( mB, K, M );
    
    % covariance matrix
    C = theta( M*K+1:end );
    C = unvech( C, 1 ); % unvech with zeros in lower triangle
    mSigma = C'*C; % ensure positive definiteness
    
% end function    
end


% SVAR counterpart to mapParameterVectorToReducedFormMatrices()
function [mA, mF] = mapParameterVectorToStructuralFormMatrices( M, K, theta )
    
    assert( eq( length(theta), M*(K+(M+1)/2) ), ...
                        'GetVarLogLikelihood:ParamVector', '' );
    
    % autoregressive coefficient matrix
    mF = theta( M*(M+1)/2+1:end );
    mF = reshape( mF, K, M );
    
    % covariance matrix
    C = theta( 1:M*(M+1)/2 );
    mA = unvech( C, 1 ); % unvech with zeros in lower triangle
    
% end function    
end




% Subfunction to map restricted parameter vector to coefficient matrices
function [mA mF] = mapRestrictedParameterVectorToCoefficientMatrices( vP, cU, cV )

    M = size( cU{1}, 1 );
    K = size( cV{1}, 1 );

    mA = zeros( M, M );
    mF = zeros( K, M );

    nnz_A = nnz( cell2mat( cU ) );
    nnz_F = nnz( cell2mat( cV ) );
    
    % There should be the same number of non-zero elements indicated by the
    % restrictions as parameters passed...
    assert( eq( nnz_A + nnz_F, length( vP ) ), ...
                    'GetVarLogLikelihood:ParamVector', '' );

    % loop over equations
    index_A = 1:M + 1;
    index_F = nnz_A + 1: nnz_A + M + 1;

    for i = 1:M
        index_A( i+1 ) = index_A( i ) + nnz( cU{i} ); % # of non-zero elements in U{i}
        mA(:,i) = cU{i}*vP( index_A(i):index_A(i+1)-1 );
        
        index_F( i+1 ) = index_F( i ) + nnz( cV{i} ); % # of non-zero elements in V{i}
        mF(:,i) = cV{i}*vP( index_F(i):index_F(i+1)-1 );
    end;
        
% end function
end

