function [mY, mZ] = ProcessVAROptions( Yraw, options, Ynames, shock_to )
% Checks that the options set for the VAR make sense, that the necessary
% ones appear, and that the variables selected for estimation exist in the
% database. Returns the variables selected for estimation in the order
% present in options.names, options.exo_names, as matrices of length T.

    constant_names = {'excluded', 'included'};
    % Code ....................................................................
    assert( ( 0 == options.constant || 1 == options.constant ), ...
        'Option "constant" must be set to 0 (exclude) or 1 (include)' );

    % M is the dimension of the VAR
    M = length(options.names);

    % The columns of the dataset in which the variables in names fall is:
    index_selected_data = zeros(1,M);
    for i = 1:M
        % endogenous variables
        if ismember( options.names{i}, Ynames )
            index_selected_data(i) = find( strcmpi( options.names{i}, Ynames ) );
        else
            error( 'Data:UserVariableList', ...
                cat(2, 'Selected variable ', options.names{i}, ...
                ' not in data set. Check "names" list.') );
        end
    end

    % K is the number of exogenous variables (not including the constant)
    K = length( options.exo_names );

    index_selected_exo = zeros(1,K);
    for i = 1:K
        if ismember( options.exo_names{i}, Ynames )
            index_selected_exo(i) = find( strcmpi( options.exo_names{i}, Ynames ) );
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

    
    % The exogenous variables are retained in the Zraw matrix.
    mZ = Yraw( :, index_selected_exo );
    % only the selected data is retained in the Yraw matrix.
    mY = Yraw( :, index_selected_data );
    
    fprintf( '\n%s\n\n', ['VAR order p = ' num2str(options.nlags) ...
            '; constant ' constant_names{options.constant+1} ] );

end