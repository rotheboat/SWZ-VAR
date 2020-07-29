function [T, M, K, N, p, constant] = GetVARDimensions( mY, mX, options )

    % Given the data matrices mY and mX, derives and returns:
    %   M : the dimension of the VAR
    %   K : the number of RHS variables per equation (including exog/deter)
    %   N : the number of exog RHS variables, incl. constant (0 if none)
    %   p : the order of the VAR
    %   constant : whether the VAR has a constant (1) or not (0)
    [ T, M ] = size( mY );
    [ T2, K ] = size( mX );
    N = 0; % init to no exogenous variables
    
    % Test for data matrix error
    assert( eq( T,T2 ), 'GetVarLogLikelihood:DataMatrix', ...
                'Y and X have an unequal number of rows' );
    
    % Check for the existence of exogenous variables Z. The number of
    % exogenous variables (not including the constant) is N. That includes
    % any lags of endogenous varibles too
    if ( isfield( options, 'has_exo_vars' ) )
        if options.has_exo_vars
            min_exo_lag = options.exo_lags(1);
            max_exo_lag = options.exo_lags(2);
            % number of exog vars x number of lags
            N = options.num_exo_vars*(max_exo_lag-min_exo_lag+1);
        end
        
    else
        options.has_exo_vars = 0;

    end
    
    % Check for existence of constant. Could be a field in options
    % struct
    if isfield( options, 'constant' )
        constant = options.constant;
        N = N + constant;
        p = (K-N)/M;
        
    else
    % If not in options struct, infer its presence from K and M
        if eq(0, rem( K-N, M ))
            constant = 0;
            p = (K-N)/M;
        else
            constant = 1;
            p = ( (K-N)-rem((K-N),M))/M;
        end

    end

end