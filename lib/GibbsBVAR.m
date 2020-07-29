function postest = GibbsBVAR( mY, mX, prior, options)
% Given data and a Sims-Zha type prior, produces posterior draws of the
% coefficients of the VAR:
%       y(t)'A = x(t)'F + e(t)'
% subject to general linear restrictions of the form:
%       a(i) = Q(i)b(i)
%       f(i) = R(i)g(i)
% as in:
% Waggoner and Zha (2003): "A Gibbs sampler for structural vector 
%   autoregressions", J. Econ. Dynamics & Control (28), 349-366.
%
% @author   Roland Meeks
% @date     2014-06-16
%
% March 2016
% ==========
% 1. Implemented Tao's estimation tips which reverse the process of finding
%   the posterior mode and drawing from the posterior distribution given in 
%   algorithm 1 (see footnote 17).
% -- Initialization from a random draw for A0, followed by 100 burn-in
%    draws.  Removes the need for an expensive initial maximization. To
%    implement, introduced an dual loop in the Gibbs sampler: the outer
%    loop for different initial A0 draws; the inner loop as usual.
% -- Record values of likelihood on each iteration, as start-point for
%    finding the posterior mode.
%   These changes also require a modification to the WZ normalization rule,
%   which needs the MLEs at one of the likelihood modes (there are many
%   reflections of the likelihood across the parameter space).  Instead of
%   applying it within-iteration, I take a second pass in which I pick the
%   value of A corresponding to the highest marginal posterior density for
%   the b coefficients (recalling A = [U{1}*b{1}|...|U{M}*b{M}])
%
%
% 2. Reorganized the processing of restrictions on the prior distribution
%   to lie in SetSimsZhaPrior().  That method now returns the restricted
%   prior density q( a_i, f_i | Q_i a_i = 0; R_i f_i = 0 ), see WZ (352).
%
% October 2016
% ============
% 1. The full log posterior kernel is now correctly computed and stored. I
%   compute log likelihood plus log prior using the methods 
%   GetPriorDensityAtParamVector() and GetDataDensityAtParamVector().
%
% 2. Added to options struct an additional, optional, setting to compute 
%   reduced Gibbs samples for the b_i (vectors of the restricted impact
%   matrix).  The idea is to produce samples for the b_i conditional on the
%   b_j (j<i) taking specified values b*_j.  The reduced Gibbs sampler
%   therefore requires an options.mA to be supplied, giving the values on
%   which to condition.
% -- The reduced Gibbs sampler will run M-2 new simulations.  The two
%    simulations it is not required to produce correspond to p(b*_1|Y,X)
%    which can be obtained from the full Gibbs sampler, and to 
%    p(b*_M|b*_1,...,b*_{M-1},Y,X) which can be obtained analytically.
% -- The reduced Gibbs sampler returns a 1 x M-2 dimensional cell array of 
%    draws.  In each element of the cell, there is a 3D array of the form
%    mA( num_draws, M, M ).
% -- Because we assumed to be starting each reduced Gibbs step at a high
%    posterior density point, number_independent_chains = 1.  Options.draws
%    and options.burn should be modified accordingly by the user.
% -- We do not require IRFs (indeed do not draw any F matrices), so
%    options.irf_horizon is set to 0.
%   
% Storage for outputs
postest = struct('A_post',[],'F_post',[],'SIGMA_post',[],'B_post',[],...
    'cH',[],'cP',[],'cS_inv',[],'cT',[],...
    'Tobs',[],'Tstar',[],...
    'IRF',[],'loglik',[],'is_restricted',false,'reject_count',0);

% Set up basic information about the VAR and check if a exogenous variables 
% are present 
[Tobs, M, K, ~, p, ~] = GetVARDimensions( mY, mX, options );

%% Dummy observations
% Dummy observations are to be used by default, unless options.dummy_obs is
% set to false
dummy_obs_flag = true;
% check if option is specified
if isfield( options, 'dummy_obs' )
    assert( islogical( options.dummy_obs ), 'GibbsBVAR:options', ...
        'dummy_obs may be true or false' );
    dummy_obs_flag = options.dummy_obs;
else
    % check the dummy observations exist (hyperparams <> 0)
    if isempty( [ prior.mYd, prior.mXd ] )
        dummy_obs_flag = false;
    end
end
if dummy_obs_flag
    [mY, mX, Tstar] = processDummyObservations( mY, mX, prior, options );
else
    Tstar = Tobs;
end

postest.Tobs = Tobs;
postest.Tstar = Tstar;

%% Statistics required for posterior computation
XpX = mX'*mX;
XpY = mX'*mY;
YpY = mY'*mY;

%% Possible reduced Gibbs sampler
do_reduced_gibbs = false; % init
% Checks for an optional argument 'reduced_gibbs' indicating that M-2
% reduced Gibbs simulations will be run instead of one full Gibbs
% simulation. Used in the Chib MDD routines.
if isfield( options, 'reduced_gibbs' )
    assert( isfield( options, 'mA' ), 'The reduced Gibbs sampler requires an options.mA' );
    % If the dimension of the model is greater than 2, the special cases
    if M <= 2
        warning( ['The reduced Gibbs sampler is not required for a model of dimension ' num2str(M)] )
        % In this instance no reduced Gibbs step is required. Return
        % control to the calling script
        return
    else
        % In this instance the reduced Gibbs routine will be run. Set index
        % limits for loop, set flag, and initialize storage.
        index_first_equation_in_loop = options.reduced_gibbs;
        do_reduced_gibbs = true;
        reducedGibbsDraws = cell(1,M-2);
    end
    options.irf_horizon = 0; % note local modify
    
else
    % We have a regular Gibbs loop, in which we wish to draw sequences of
    % b_i coefficients from all equations i = 1,...,M
    index_first_equation_in_loop = 1;
end

fprintf( ['\n--- Waggoner-Zha Sampler for equations ' num2str(index_first_equation_in_loop) ...
    ' thru ' num2str(M) ' ---'] );

%% Info
if options.reject_nonstationary
    msg_1 = 'rejecting non-stationary draws ';
else
    msg_1 = 'keeping all posterior draws ';
end

if 0 < options.irf_horizon
    msg_2 = 'computing IRFs';
else
    msg_2 = '';
end

% If no burn-in option is present set to 100 burn-in draws (can be made
% zero by choosing options.burn = 0).
if ~isfield( options, 'burn' )
    options.burn = 100;
end

if ~isfield( options, 'number_independent_chains' )
    number_independent_chains = 5;
else
    number_independent_chains = floor( options.number_independent_chains );
    assert( number_independent_chains >= 1, ...
        'Error: Option "number_independent_chains" must be at least one' );
end

if isfield( options, 'quit_multiple' )
    max_draws = number_independent_chains*options.draws*options.quit_multiple;
else
    max_draws = number_independent_chains*options.draws*1e5;
end

if ~isfield( options, 'it_print' )
    it_print = 5e2;
else
    assert( isnumeric( options.it_print ), 'GibbsBVAR:Option it_print should be an integer' );
    it_print = options.it_print;
end


% ---- Get prior parameters ----------------------------------
cHtld_inv = prior.Htld_inv;
cPtld = prior.Ptld;
cStld_inv = prior.Stld_inv;
cStld = prior.Stld;
cV = prior.cV;
cU = prior.cU;

% ---- Compute parameters of the posterior ------------------
% Cell arrays of posterior matrices, WZ	(12) and (13)
cH = cellfun( @(Htld_inv,V) inv( V'*XpX*V + Htld_inv ), ...
              cHtld_inv, cV, ...
              'UniformOutput', false );
          
cP = cellfun( @(H,Htld_inv,V,U,Ptld) H*(V'*XpY*U + Htld_inv*Ptld ), ...
                cH, cHtld_inv, cV, cU, cPtld, ...
                'UniformOutput', false );

% Use mldivide to compute the inversion required for S matrix
cS_inv = cellfun( @(U, Stld_inv, P, Ptld, H, Htld_inv) (1/Tstar)*( U'*YpY*U + Stld_inv + Ptld'*Htld_inv*Ptld - P'*(H\P) ), ...
				   cU,cStld_inv,cP,cPtld,cH,cHtld_inv, ...
                'UniformOutput', false );
% BUGBUG why is this line here:
% cellfun( @(S_inv) chol( S_inv ), cS_inv, 'UniformOutput', false );
% as nothing is assigned to the output returned by the call?
% BUGBUG U'*YpY*U + Stld_inv + Ptld'*Htld_inv*Ptld - P'*(H\P) is numerically unstable
% Could be that Ptld'*Htld_inv*Ptld and P'*(H\P) are about the same size?
			
% There are two equivalent ways of inverting S_inv to obtain S...
%cS = cellfun( @(S_inv) S_inv\eye(size(S_inv)), cS_inv, 'UniformOutput', false );
cS = cellfun( @(S_inv) inv( S_inv ), cS_inv, 'UniformOutput', false );

% Defined by WZ (p. 355)
cT = cellfun( @(S) chol(S,'lower'), cS, 'UniformOutput', false );

% Store posterior parameters in struct
postest.cH = cH;
postest.cP = cP;
postest.cS_inv = cS_inv;
postest.cT = cT;

% #BUGBUG WHICH WAY OF COMPUTING cT IS CORRECT?#
% mT is the cholesky factor of mS, i.e. mS=mT*mT'. cT is a cell array 
    % of mT matrices, corresponding to each of the equations in the model
%    cT = cellfun( @(S_inv) inv( chol( S_inv, 'lower' ) ), ...
%       cS_inv, 'UniformOutput', false );

						
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Gibbs sampler ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Gibbs-related preliminaries

% Preallocate storage for columns of the impact matrix A, denoted vb, 
% that are computed each iteration of the Gibbs loop. There will be M 
% of these per iteration. Each loop, the contents of the array will be
% overwritten with the new draw.  In the case of a reduced Gibbs run, the
% first options.reduced_gibbs columns of the impact matrix are not going to
% change across draws.  In this case, pre-populate the cell array.
if do_reduced_gibbs
    % the reduced Gibbs run case
    ca = mat2cell(options.mA,M,ones(1,M));
    % cell array columns a -> cell array b = U'*a
    cb = cellfun(@(a,U) U'*a, ca, prior.cU, 'UniformOutput', false );
else
    % otherwise just initialize with an empty cell array
    cb = cell(1,M);
end
    
% Ditto for cg.
cg = cell(1,M);

nswitch = 0; % for nmlzvar.m
% set simulation options

%% ---- start the loop to sample A matrices ----
fprintf( '\rStarting to draw A matrices...\r' );
tic;

% Set the initial 'A' or 'A0' matrix to the MLE.  Note that I call this
% matrix mUb, as it is comprised of the entries [ U(1)*b(1) ... U(M)*b(M) ]
if isfield( options, 'mA' )
    mUb = options.mA;
    % In the case where a single, fixed initial A0 is used (e.g. the peak
    % of the posterior density), we will run only one chain.
    number_independent_chains = 1;
end

%% Initialize storage
% Total number of draws
ntot = number_independent_chains*( options.draws + options.burn );  

% Initialize storage
% Structural impact matrices
mmA = NaN( number_independent_chains*options.draws, M, M );
% Structural lag matrices
mmF = NaN( number_independent_chains*options.draws,  K, M );
% Reduced form parameter matrices
mmB = NaN( number_independent_chains*options.draws,  K, M );
% Kernel of the log marginal posterior density for the impact matrix
mmPDFb = NaN( number_independent_chains*options.draws, 1 );
% Kernel of the log posterior density
mmPostKer = NaN(number_independent_chains*options.draws,1);
% Impulse response functions
imp_resp = NaN( number_independent_chains*options.draws, M, M, options.irf_horizon );


drawcount = 1;
savecount = 1;

% Outer loop: case of a random initial A matrix (options.mA is absent)
%  > if options.mA is set, number_independent_chains = 1;
for outer_loop = 1:number_independent_chains

    % If no MLE is present (the default) then take a draw from the prior
    % distribution for A.  Columns a_i = U_i b_i and b_i ~ N( 0, Stld_i ).
    % The mvnrnd function checks the covariance argument is symmetric and
    % positive definite.  The method SetSimsZhaPrior::GetDerivedPriorUnderRestrictions
    % therefore ensures that this is the case by using the product of
    % Cholesky factors to construct the matrices in Stld_i.
    if ~isfield( options, 'mA' )
        ca_prior_draw = cellfun( @(U,S) U*mvnrnd( zeros(size(U,2),1), S )', ...
            cU, cStld, 'UniformOutput', false);
        mUb = cell2mat( ca_prior_draw );
    end
      
    % Inner loop: conditional on the randomly drawn A matrix, sample from
    % posterior
    while drawcount <= outer_loop*( options.draws + options.burn )
        % Show progress...
        if eq(0, mod(drawcount,it_print) )
            fprintf( '\nIteration = %d\n', drawcount );
            toc;
        end
        % Display a dot every it_print/100 draws...
        if ( 0 == mod(drawcount,ceil(it_print/50)) )
            fprintf( '.' );
        end


        % #NOTE Moved outside the loop RM 06/10/16
        % cT = cellfun( @(S)chol(S,'lower'), cS, 'UniformOutput', false );

        % Loop over equations.  The Gibbs sampler gives us draws of b(i*) conditional
        % on { b(1),...,b(i*-1),b(i*+1),...,b(n) }.  The vector b(i*), i.e. the 
        % impact vector in the i*th equation, is seen to have q(i*) free parameters.
        % As a(i) = U(i)*b(i), lenk = size( U{i}, 2 );
        for k = index_first_equation_in_loop:M
            % Get the vectors w(i) needed to construct b(i) as per WZ (15)
            %mJ = cT{k}*mUb;
            mJ = mUb;
            mJ(:,k) = 0;
            % Let sw be a non-zero Mx1 vector perpendicular to each vector
            % in {mU(i)b(i)|i\neq i*}, for 1 <= i* <= M.
            sw = null( mJ' );
            mW = cT{k}'*cU{k}'*sw/norm( cT{k}'*cU{k}'*sw ); % w(1), WZ (14)

			% Now compute the remaining w_j vectors, j=2,...,q_i_*, to produce a 
			% set of vectors {w_1,...,w_q_i_*} that are an orthonormal basis for R^{q_i_*}
            lenk = size( cU{k}, 2 );
            for j=1:lenk
                [mW, ~] = qr( [mW unit(lenk,j)], 0 );  % Economy size qr...
            end
            % Produces w(j) = mW(:,j) as required in WZ (p. 357)

            % Now obtain the beta's needed to construct b(i) as per WZ (15)
            % TZ: draw independent beta's that combine uR to form a's (columns of A0)
            % Notice that this gives you a WZ \beta_j, 2 <= j <= q_i_*, as per Theorem 1 (b)
            % i.e. normally distributed: mean 0 variance 1/T.  In TZ code, the indexes of
            % the \beta's are 1 <= j <= q_i_* - 1, and beta(q_i_*) is computed
            % next (p. 355).
            %
            % ** I found that it is essential to keep the indexing of the betas
            % exactly as in the WZ paper -- beta(1) must be weighted by 
            % w(1) = mW(:,1).  If instead I follow TZ's code, I get some
            % obviously silly draws for A0 -- far from the MLEs with bimodality
            % in the off-diagonals. ** RM (2014-06-13)
            beta = zeros(lenk,1);
            beta(2:lenk) = sqrt(1/Tstar)*randn(lenk-1,1);

            % gamma or 1-d Wishart draw
            % I have assumed that TZ kdf=Tobs+Tdummy
            % Draw x from p(x) by letting r = x^2
            r = sqrt(1/Tstar)*randn(Tstar+1,1);

            % beta(q_i_*):
            % For each draw, assign x=+/-sqrt(r) each with prob 1/2 to obtain a draw from WZ 
            % (16). Notice we assign the result to the last element of beta although in the 
            % WZ paper it is \beta_1 that has this non-standard distribution
            if rand(1)<0.5
               beta(1) = sqrt(r'*r);
            else
               beta(1) = -sqrt(r'*r);
            end

            % Form k^th equation of A0 = [U(1)b(1),..., U(M)b(M)] matrix;  
            vb = cT{k}*mW*beta;
            % store vector of b coefficients in a cell array
            cb{k} = vb;
            % fill in the k^th column of the impact matrix
            mUb(:,k) = cU{k}*vb;

        % end loop over equations
        end

        
        if drawcount > outer_loop*options.burn + (outer_loop-1)*options.draws
            % Time to start saving the draws; to understand the condition
            % for when to start saving, think as follows:
            %
            % outer loop    inner loop      total draws
            %       1       burn            burn
            %       2       burn            2*burn + draws
            %       .
            %       n       burn            n*burn + (n-1)*draws
            %
            % The counter savecount is only clicked once the burn-in period
            % has been passed inside each inner loop
            mmA(savecount,:,:) = mUb;

            % We also need to save the kernel of the marginal posterior for the
            % matrix [b(i)]_i (in logarithms)
            % BUGBUG Tobs or Tstar?
            sumQuadB = sum( cellfun( @(b,iS) b'*iS*b, cb, ...
                               prior.Stld_inv, 'UniformOutput', true ) );
            mmPDFb(savecount,:) = ...
                Tobs*log( abs( det( mUb ) ) ) - Tobs*( sumQuadB )/2;
            
            savecount = savecount+1;
        end

        % We're done drawing for this iteration on the inner loop
        drawcount = drawcount+1;
    end  % end inner loop
end % end outer loop
    
%% Normalization step

% Find the largest value of the marginal posterior density of (12). If an
% mA matrix has been supplied in options, assume that this is the posterior
% mode and skip this step.
if ~isfield( options, 'mA' )
    [~,max_index] = max( mmPDFb );
    % Select the draw that maximizes (12) as the "peak" likelihood value
    % against which to normalize
    mA_mode = squeeze( mmA(max_index,:,:) );
else
    mA_mode = options.mA;
end
    
% Call to Waggoner-Zha (JEcotx 2003) normalization routine.  Option 1 
% is the ML distance rule (Algorithm 1).  Requires nmlzvar.m from Tao
% Zha's matlab library.
for good_draw = 1:number_independent_chains*options.draws
    mUb = squeeze( mmA(good_draw,:,:) );
	[mUb,nswitch] = nmlzvar(mUb,mA_mode,[],1,nswitch,[]);
	mUb(abs(mUb)<eps) = 0;  % Avoid cumulative round-off error
    % Replace the draw for A with its normalized counterpart
    mmA(good_draw,:,:) = mUb;
end

fprintf( '\nDone drawing A matrices... starting F draws:\n\t>%s \n\t>%s\n', msg_1, msg_2 );
toc;

% If we are running a reduced Gibbs sampler, then we can end here
if do_reduced_gibbs
    postest.reducedGibbsDraws = reducedGibbsDraws;
    postest.draws = savecount-1;
    postest.A_post = mmA(1:savecount-1,:,:);
    return
end


%% ---- given draws of A, start sampling F matrices ----
% Assuming that the function has not already returned control to the
% caller, because a reduced Gibbs sampler is being run, go ahead and draw
% from the density p(g|b).

% set the largest eigenvalue of the companion matrix that will be
% considered 'stationary' when sampling F matrices:
max_stationary_eigenvalue = 1 - eps;

% initialize
mF = NaN(K,M);
drawcount = 1;
flag = 1; % if options.reject_nonstationary = false, flag is aways 1
rejectcount = 0;

while drawcount <= number_independent_chains*options.draws
    % Show progress...
    if ( (1 == flag) && (0 == mod(drawcount,it_print) ) )
        fprintf( '\nIteration = %d\n', drawcount );
        toc;
    end
    % Display a dot every it_print/100 draws...
    if ( 0 == mod(drawcount,ceil(it_print/50)) )
        fprintf( '.' );
    end
    % Quit if we're going nowhere
    if ( drawcount + rejectcount > max_draws )
        fprintf( '\r' );
        warning( 'GibbsBVAR:Fdraws', ...
            'Maximum number of draws exceeded.' );
        break;
    end
    
    % Retrieve the A draw
    mA = squeeze( mmA(drawcount,:,:) );
    % make a cell array of post-normalization impulse vecs (cols of A)
    ca = mat2cell(mA, M, ones(M,1) ); 
    % transform back to b's, using b(i) = U(i)'a(i)
    cb = cellfun(@(U,a) U'*a, cU, ca, 'UniformOutput', false );
 
    % Loop over structural equations
    for k = 1:M
        % Conditional on the draw for b(i) we can draw g(i) as per WZ (13)
        % Draw g(i)|b(i) ~ normal (WZ 13)
        mHsqrt = chol( cH{k} );
        %vg = cP{k}*cb{k} + mHsqrt*randn( size(cV{k},2), 1 );
        vg = mvnrnd( cP{k}*cb{k}, cH{k} )';
        % store vector of b coefficients in a cell array
        cg{k} = vg;
        % transform g vector to f vector using f(i) = V(i)g(i)
        mF(:,k) = cV{k}*vg;
    
    end % end loop over equations
    
    % Check the model is stable at the newly drawn parameter values if the
    % option reject_nonstationary is set to true, otherwise just go ahead
    % and save the structural and reduced form matrices every time.
    %   flag = 1 if the draw is accepted as stationary
    %   flag = 0 if the draw is rejected as non-stationary
    % The test for non-stationarity is based on the reduced form AR matrix
    mB = mF/mA;
    if ( options.reject_nonstationary )
        flag = ~any( abs( eig( GetCompanionForm( mB, options ) )' ) > ...
            max_stationary_eigenvalue );
    end
    % *** POSSIBLY NO NEED TO COMPUTE mB WITHIN THE LOOP? *** 
    % IS y' = x'*B + u' STATIONARY IFF y'*A = x'*F + e' STATIONARY
    
    if flag
        % retain the structural RHS matrix
        mmF( drawcount,:,: ) = mF;
        % retain the reduced form RHS matrix
        mmB( drawcount,:,: ) = mB;
        % retain the estimate of the posterior kernel
        %   posterior kernel = ln(likelihood) + ln( prior )
        log_data_density_at_point = GetDataDensityAtParamVector(mY,mX,mA,mF);
        log_prior_density_at_point = GetPriorDensityAtParamVector( prior, cb, cg );
        mmPostKer(drawcount) = log_data_density_at_point + log_prior_density_at_point;
    else
        % See "Quit if we're going nowhere", above
        rejectcount = rejectcount+1;
    end        

    % If irf_horizon is > 0 then compute irf
    if ( flag && 0 < options.irf_horizon )
        imp_resp(drawcount,:,:,:) = ...
            impulse( reshape( mB(1:M*p,:)', M, M, p ), ...
                inv( mA ), options.irf_horizon );
    end


    % We keep on re-drawing mF until flag = 1
    drawcount = drawcount + flag;
    
end % end sampling F matrices
fprintf( '\nDone with F draws!\n' );

toc;

%% Outputs
postest.A_post = mmA(1:drawcount-1,:,:);
postest.F_post = mmF(1:drawcount-1,:,:);
postest.B_post = mmB(1:drawcount-1,:,:);
postest.loglik = mmPostKer(1:drawcount-1,:);
postest.nswitch = nswitch;
postest.IRF = imp_resp(1:drawcount-1,:,:,:);
postest.reject_count = rejectcount;
postest.draws = drawcount-1;

% function end	
end



% Utility function returns a column vector of length n with a unit entry at the
% jth position, and zeros elsewhere.  Used in the main Gibbs loop as part
% of the QR decomposition.
function x = unit(n,j)
	x = eye( n );
	x = x(:,j);
end



function [mYstar, mXstar, Tstar] = processDummyObservations( mY, mX, prior, options )
% Helper function that appends the dummy observations passed as part of the
% prior structure (prior.mYd, prior.mXd) to the data matrices
% (mY, mX).
    mYd = prior.mYd;
    mXd = prior.mXd;
    
    mYstar = [mYd; mY];
    mXstar = [mXd; mX];
    
    [ Tstar, ~ ] = size( mYstar );

    fprintf( '\nUsing dummy prior(s):\r\t>lambda(5) = %d\r\t>lambda(6) = %d', ...
        options.lambda(5), options.lambda(6) );
end