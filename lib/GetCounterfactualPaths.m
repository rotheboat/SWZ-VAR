% Computes the paths of endogenous variables under alternative,
% counterfactual, parameter values. 
%
%
function data = GetCounterfactualPaths( counterfactual, Y, X )

    disp('Retrieving counterfactual paths.');
    
    % init
    [T M] = size( Y );
    K = size( X, 2 );
    constant = rem( K, M );
    p = floor( K/M );
    
    if (constant)
        inter = 1;
    else
        inter = [];
    end
    
    Y_new = zeros( T, M );
    X_new = zeros( T, constant+M*p );
    X_new(1,:) = X(1,:);
    
    % Not all the posterior draws may have been kept; check how many
    % counterfactual paths are available
    ndraws = counterfactual.num_draws;
    
    % Create storage for counterfactual Y paths
    data = zeros( ndraws, T, M );
    
    % Loop over draws
    for draw = 1:ndraws
        A_new = squeeze( counterfactual.A_new( draw, :, : ) );
        B_new = squeeze( counterfactual.B_new( draw, :, : ) );
        V = squeeze( counterfactual.V( draw, :, : ) );
        
        for t = 1:T
            Y_new(t,:) = ( X_new(t,:)*B_new + V(t,:) )/A_new;
            X_new(t+1,:) = [ inter, Y_new(t,:), X_new(t,constant+1:constant+M*(p-1)) ];
        end % end for t time periods
        
        data( draw, :, : ) = Y_new;
    end % end for draw
    
end