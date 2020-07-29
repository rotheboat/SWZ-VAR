function [ Bhat, Ghat ] = FindPosteriorMode( priorparams, postest )
% Find the peak of the joint posterior density of the Sims/Waggoner/Zha
% type Bayesian VAR.  

% Extract the restriction matrices
cU = priorparams.cU;
cV = priorparams.cV;
cS = priorparams.cS;
cP = priorparams.cP;
cH = priorparams.cH;


% Get parameters at the highest posterior draw
[~, index] = max( postest.loglik );
mA_mode = squeeze( postest.A_post( index, :, : ) );

% Break the mA matrix into a cell array of column vectors
cA = mat2cell( mA_mode, M, ones(M,1) );

% Perform elementwise operation cU{k}'*ca{k} to transform column vector
% elements ca into restricted parameter space cb
cb = cellfun( @(U,a), U'*a, cU, cA, 'UniformOutput', false );
 
Bhat = FindFunctionMaximum( ...
    @(mB) GetMarginalPosteriorB( mB, X, Y ), cb );

% Conditional on Bhat, look for maximimum for 
Ghat = FindFunctionMaximum( ...
    @(mG) GetMarginalPosteriorB( mG, Bhat, X, Y ), init_mG );

end


function param_vector = FindFunctionMaximum( func_to_maximize, init_val )
% This function sets up and executes a numerical optimization using various
% options

param_vector = fminsearch( func_to_maximize, init_val );

end


function lik = GetMarginalPosteriorB( mB, X, Y )
% Given data X and Y, return the marginal posterior density for the LHS
% restricted coefficient matrix B evaluated at mB

end

function lik = GetMarginalPosteriorG( mG, mB, X, Y )
% Conditional on Bhat, and given the data, return the marginal posterior
% density for the RHS possibly restricted coefficient matrix G evaluated at
% mG

end

