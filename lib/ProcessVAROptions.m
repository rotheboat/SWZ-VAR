function [mY, mZ] = ProcessVAROptions( data_table, options, Ynames, shock_to )
% Checks that the options set for the VAR make sense, that the necessary
% ones appear, and that the variables selected for estimation exist in the
% database. Returns the variables selected for estimation in the order
% present in options.names, options.exo_names, as matrices of length T.

    constant_names = {'excluded', 'included'};
    % Code ....................................................................
    assert( ( 0 == options.constant || 1 == options.constant ), ...
        'Option "constant" must be set to 0 (exclude) or 1 (include)' );

    % The height of the data table is the number of observations
    T = height( data_table );
    
    % M is the dimension of the VAR
    M = length(options.names);
    mY = nan( T, M );

    % Loop over the variables in names (the endogenous variables to be
    % modeled in the VAR)
    for i = 1:M
        % endogenous variables
        if ismember( options.names{i}, Ynames )
            mY(:,i) = data_table.( options.names{i} );
        else
            error( 'Data:UserVariableList', ...
                cat(2, 'Selected variable ', options.names{i}, ...
                ' not in data set. Check "names" list.') );
        end
    end

    % K is the number of exogenous variables (not including the constant)
    K = length( options.exo_names );
    mZ = nan( T, K );
    
    % Loop over the variables in exo_names (the exogenous variables to be 
    % included in the VAR, aside from the constant) (may be empty)
    for i = 1:K
        if ismember( options.exo_names{i}, Ynames )
            mZ(:,i) = data_table.( options.exo_names{i} );
        else
            error( 'Data:UserExoVariableList', ...
                'Selected variable not in data set. Check "names" list.' );
        end
    end

    % Process shock_to option
    for j = 1:numel( shock_to )
        if ~ismember( shock_to{j}, options.names )
            error( 'UserOptions:shock_to', ...
                ['Impulse to ' shock_to{j} ' equation: Selected variable not in data set.'] );
        end
    end
    % End process shock_to option

    
    fprintf( '\n%s\n\n', ['VAR order p = ' num2str(options.nlags) ...
            '; constant ' constant_names{options.constant+1} ] );

end