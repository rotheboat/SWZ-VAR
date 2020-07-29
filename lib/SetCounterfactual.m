%% SetCounterfactual - set counterfactual policy rule
% Takes as arguments parameter draws from posterior simulator, returns new
% 'counterfactual' parameter matrices.  
%
%
% ctrfact.theta     the vector of counterfactual policy parameters is
%                   called theta; there may be one policy rule for each
%                   draw, in which case this is a matrix
%
% @author   Roland Meeks
% @date     March 2014
function counterfactual = SetCounterfactual( names, postest, ctrfact, Y, X  )

% init
counterfactual = struct( 'A_new', [], 'B_new', [], 'Phi_new', [], ...
    'V', [], 'num_draws', [] );
[T M] = size(Y);
[T K] = size(X);

if eq(0, rem( K, M ))
    constant = 0;
    p = K/M;
else
    constant = 1;
    p = (K-rem(K,M))/M;
end

% Create storage
ndraws = size( postest.A_post, 1 );
A_new = zeros( ndraws, M, M );
B_new = zeros( ndraws, K, M );
Phi_new = zeros( ndraws, K, M );
% The structural errors are recalculated on each draw
V_d = zeros( ndraws, T, M );

% Find out where the policy instrument is in the data set
index_policy = strmatch( ctrfact.instrument_name, names, 'exact' );

if isempty( index_policy )
    warning('GetPolicyRule: No such policy instrument in model.');
    return;
end;

% Check what type of policy rule is being supplied:
if strcmp( 'single', ctrfact.policy_type )
    % supply a single policy rule
    % Replicate the rule ndraws times
    ctrfact.theta = kron( ones( ndraws, 1 ), ctrfact.theta(1,:) ); 
elseif strcmp( 'multiple', ctrfact.policy_type )
    % supply an ndraws x M+K matrix of policy rule
    % reject if the matrix of rules is the wrong size
    if ~( size( ctrfact.theta, 1 ) == ndraws )
        warning( 'Option `policy_type` set to `multiple`.' );
        return;
    end
elseif strcmp( 'random', ctrfact.policy_type )
    % create an ndraws x M+K matrix of random policy rules
else
    warning( ['No counterfactual produced: option `policy_type` must' ...
        'be one of `single`, `multiple`, `random`.' ]);
    return;
end

if isempty( postest.SIGMA_post )
    model_type = 'structural_form';
else
    model_type = 'reduced_form';
end;

% Start looping over draws...
counter = 1; % keep track of number of draws that are kept

for draw = 1:ndraws

    switch model_type
        case 'structural_form'
            A_cf = squeeze( postest.A_post( draw, :, : ) );
            B_cf = squeeze( postest.F_post( draw, :, : ) );
            % (1) If there is a constant, will need the structural shocks for the
            % initial period, retrived at the original parameter values
            V = Y*A_cf - X*B_cf;
            
        case 'reduced_form'
            % Retrieve parameter matrices for this draw
            C = chol( squeeze( postest.SIGMA_post( draw,:,: ) ), 'upper' );
            Phi = squeeze( postest.A_post( draw,:,: ) );

            % (1) If there is a constant, will need the structural shocks for the
            % initial period, retrived at the original parameter values
            V = (Y - X*Phi)/C;
            
            % Counterfactual A and B
            A_cf = eye(M)/C;
            B_cf = Phi/C;
        otherwise
            % do nothing
    end;
    % (2) Retrieve, or generate, the vector of policy parameters for this
    % draw
    theta = ctrfact.theta( draw, : );
    
    % (3) Map counterfactual parameters into A, B matrices
    A_cf(:, index_policy) = theta(1:M);
    B_cf(:, index_policy) = -theta(M+1:end);

       
    % (4) If there's a constant in policy rule, ensure that it is set so 
    %     that Y_new(1,i)=Y(1,i)
    if (constant)
        mu = Y(1,:)*A_cf - X(1,2:end)*B_cf(2:end,:) - V(1,:);
        B_cf( 1, : ) = mu;
    end
    % end counterfactual

    % (5) Set the reduced form parameter matrix 
    Phi_cf = B_cf/A_cf;
    
    % (6) Check the model is stable at the new parameter values by checking
    % the roots of the counterfactual companion matrix
    if ( ctrfact.check_stability )
        keep = ~any( ge( abs( eig( [Phi_cf(constant+1:end,:), eye(M*p,M) ] ) )', 1) );
    else
        keep = 1;
    end
    
    % store the draws that are OK in the relevant arrays
    if ( keep )
        A_new(counter,:,:) = A_cf;
        B_new(counter,:,:) = B_cf;
        Phi_new(counter,:,:) = Phi_cf;
        V_d(counter,:,:) = V;
        
        counter = counter+1;
    end % end if keep
    
end % end draw

fprintf( '\nCounterfactual: %d draws retained out of %d.\n', counter, ndraws );
counterfactual.num_draws = counter-1;

% Drop any zero rows
A_new( counter:end,:,:) = [];
B_new( counter:end,:,:) = [];
Phi_new( counter:end,:,:) = [];
V_d( counter:end,:,:) = [];

% Set the values of the struct for return
counterfactual.A_new = A_new;
counterfactual.B_new = B_new;
counterfactual.Phi_new = Phi_new;
counterfactual.V = V_d;

end


%%
% When policy is random, get draws according to some distribution...
function drawPolicyRule( param_mean, param_covar )

end
