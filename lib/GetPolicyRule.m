% Designate one variable in a VAR model the policy instrument, then get the
% contemporaneous and lagged coefficients in the corresponding policy rule
% for every posterior draw.
%
% @author   Roland Meeks
% @date     March 2014
% @return   policy_params : an ndraws x M matrix, with the parameters of 
%               the policy function in rows.  Model assumed to be of the
%               form Y*A + X*B = V, so be aware of signs
%           policy_params_norm : as policy_params, but normalized so the
%               coefficient on the policy instrument is 1
% @arguments
% prior_options = a structure containing a wealth of info on the model
% postest = a structure returned from GetINWishartPosterior
% instrument_name = something like 'BANKRATE'
function [ policy_params policy_params_norm ] = GetPolicyRule( prior_options, postest, instrument_name, varargin )

fprintf('\n');

% init
names = prior_options.names;
M = prior_options.M;
p = prior_options.p;
K = prior_options.K;

if ( K > M*p )
    constant = 1;
else
    constant = 0;
end;


% Which element of Y is the policy variable of interest?
index_policy = strmatch( instrument_name, names, 'exact' );

if isempty( index_policy )
    warning('GetPolicyRule: No such policy instrument in model.');
    return;
end;

ndraws = size( postest.A_post, 1 );

% Matrix to hold the policy parameter vector corresponding to each
% posterior draw.  Draws in rows.  Variables in columns.
policy_params = zeros( ndraws, M*(p+1) + constant );

% vector of structural intercepts (if constant present)
mu = [];

for draw = 1:ndraws

    if ( isempty( postest.SIGMA_post ) )
        A0 = squeeze( postest.A_post( draw,:,: ) );
        B = squeeze( postest.F_post( draw,constant+1:end,: ) );

        if (constant)
            mu = squeeze( postest.F_post( draw,1,: ) )';
        end;

    % Backwards compatibility    
    else 
        C = chol( squeeze( postest.SIGMA_post( draw,:,: ) ), 'upper' );
        Phi = squeeze( postest.A_post( draw,constant+1:end,: ) );
        A0 = inv( C );
        B = Phi/C;

        if (constant)
            mu = squeeze( postest.A_post( draw,1,: ) )'/C;
        end;
    end;
    
    if constant
        policy_params( draw, : ) = [ A0(:,index_policy); -mu( index_policy ); -B(:,index_policy) ]';
    else
        policy_params( draw, : ) = [ A0(:,index_policy); -B(:,index_policy) ]';
    end
%     % Outputs the policy parameter vector for every draw.
%     if( 0 )
%         norm_const = A0(index_policy,index_policy);
%         
%         if (constant)
%             fprintf( '%d:  0 = %12.3f * \t %s \n', draw, mu( index_policy ), 'constant' );
%         else        
%             disp( '%d:  0 = ', draw );
%         end;
%         for lag = 0:p,
%             for j = 1:M,
%                 fprintf( '%d: %8.3f * \t %s (L^%d) \n', M*lag + j, ...
%                     policy_params( draw, M*lag + j )/norm_const , ...
%                     names{j}, lag );
%             end;
%         end;
%         fprintf( '\n' );
%     end;

end

% Normalize each row of the policy parameter matrix with respect to the 
% coefficient on the policy variable.
policy_params_norm = diag(1./policy_params(:,index_policy))*policy_params;

if (4 == nargin)
    if ( 1 == varargin{1} )
        PlotParameterCDF( policy_params_norm );
    end
end


% Take the average of the policy parameters from each draw.  
mean_policy_params = mean( policy_params_norm );

if ( constant )
    fprintf( 'Estimated policy rule (mean):  0 = %6.3f * \t %s \n', ...
        mean_policy_params( M+1 ), 'constant' );
else
    fprintf( 'Estimated policy rule (mean):  0 = \n' );
end;
for lag = 0:p,
    for j = 1:M,
        if (0 == lag)
            fprintf( '%d: %8.3f * \t %s (L^%d) \n', j, ...
                mean_policy_params( :, j ), ...
                names{j}, lag );
        else
            fprintf( '%d: %8.3f * \t %s (L^%d) \n', constant+M*lag+j, ...
                mean_policy_params( :, constant+M*lag+j ), ...
                names{j}, lag );            
        end
    end
end




end
%% end of function --------------------------------------------------------

function PlotParameterCDF( policy_params_norm )

num_vars_to_plot = size( policy_params_norm, 2 ); 
plot_rows = floor( sqrt(num_vars_to_plot) );
plot_cols = ceil( sqrt(num_vars_to_plot) );
if ( num_vars_to_plot > plot_rows*plot_cols )
    plot_rows = plot_cols;
end;

% new figure
figure; 

for i=1:num_vars_to_plot
        [f, xi] = ksdensity( policy_params_norm( :, i ), 'function', 'cdf' );
        subplot( plot_rows, plot_cols, i );
        plot( xi, f, 'b' );
end % end for 

end % end function
