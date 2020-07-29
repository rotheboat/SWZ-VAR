%% Impose 'soft' zero restrictions on specified banking variables in the
% macro block of the VAR. Specify tightness with which the restrictions
% are imposed.
%
% @arguments [struct] priorparams - information on the VAR prior parameters
%            [vector] variable_index - position of variables to be excluded
%            [scalar] tightness - larger numbers mean a tighter prior
%
% @author    Roland Meeks
% @date      2013-01-23
function s = GetExcludeBankVarFromMacroBlockPrior( priorparams, ...
                                            prior_options, banks_index )

% The return object (copy the priorparams struct)
s = priorparams;

%---------------------Get dimensions of VAR----------------------------                               
  
K = prior_options.K;
M = prior_options.M;
p = prior_options.p;                    

%------------------Check if BVAR includes an intercept-----------------
constant = 0;


if ( K > M*p ) %( 1 == mod( K, M*p ) )
    constant = 1;
elseif ( K < M*p ) %( 1 < mod( K, M*p ) )
    % throw an exception
end;                               
                               
% The variables which are not in banks_index are in the macro block, by
% definition. Get index of macro block (q):
q = [1:M]';
q( banks_index ) = [];

% The positions of the restricted variables (r)
r = prior_options.soft_exclude;

% The dimensions of r and q
n_r = length( r );
n_q = length( q );

%-----------------------Set the prior A matrix-----------------------------
A = priorparams.A_prior(constant+1:end,:);
alpha = A(:);

% See documentation for construction of this expression, which selects the
% elements of 
index_select = kron( p*M*(q-ones(n_q,1)), ones(p*n_r,1) ) + ...
                kron( ones(n_q,1), ( M*kron( tril(ones(p,p),-1)*ones(p,1), ones(n_r,1) ) + ...
                 kron( ones(p,1),r) ) );

alpha( index_select ) = 0;
A = reshape( alpha, M*p, M );
% Assign the modified prior A matrix to the struct
s.A_prior( constant+1:end,: ) = A;


%----------------------Set prior precision matrix V------------------------

V = priorparams.V_prior;

% this is a M*(constant+M*p) x M*(constant+M*p) dimensional object
v1 = 1:M*K; % index vector for all coefficients
    
% V.1. Where a constant is present, get the rows corresponding only to the 
% autoregressive coefficients
if ( 1 == constant )
    % drop the rows corresponding to the constant (1, 2+M*p, 3+2*M*p,...)
    const_select = tril( ones(M,M) )*ones(M,1) + ...
                            tril( ones(M,M), -1)*ones(M,1)*M*p;
    v1( const_select ) = []; % drop constants
    V = V( v1, v1 ) ;   % select only rows/cols
end;

% To apply restriction to diagonal elements of V, need to create a
% selection matrix the same dimension as V...
D1 = zeros( M*M*p, 1 );
D1( index_select ) = index_select ;        
% create a temporary matrix which has D1 on its diagonal
D2 = diag( D1 );
V( 0 ~= D2 ) = V( 0 ~= D2 )*1/(prior_options.kappa);

% Assign the modified prior V matrix to the struct
s.V_prior( v1, v1 ) = V;

            
end