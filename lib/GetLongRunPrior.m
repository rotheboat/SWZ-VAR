%% 
% @author       Roland Meeks
% @date         2012-08-03
% @see          GetAnalyticPrior.m
%
function s = GetLongRunPrior( s, data_file, names, p, constant, varargin )

%% INIT
rel_tight = 1e0; % default relative tightness
cal = DateHandler(); % instantiate DateHandler class

if strcmp( 'Windows_NT', getenv('OS') )
    [num text raw] = xlsread( data_file, 'longrun' );
else
    load( [ data_file '_longrun.mat' ] );    
end;


% edit calendar settings to match database (**not** estimation dates)
START_YEAR = 1963; START_Q = 2; % DO NOT EDIT
END_YEAR = 2012; END_Q = 2;     % UNLESS EXCEL CHGD
FREQ = 4;

% edit to set estimation dates
EST_START_YEAR = 1963; EST_START_Q = 3;
EST_END_YEAR = 1989; EST_END_Q = 3;


% varargin is an options struct that can override these settings, if
% present
if ( 0 < size( varargin, 2 ) )
    options = varargin{1};
    
    % Relative tightness option
    if ~isempty( options.long_run_tightness )
        rel_tight = options.long_run_tightness;
    end;
    
    % Estimation dates options
    if ~isempty( options.LR_START_YEAR )
        EST_START_YEAR = options.LR_START_YEAR;
    end;
    if ~isempty( options.LR_START_Q )
        EST_START_Q = options.LR_START_Q;
    end;
    if ~isempty( options.LR_END_YEAR )
        EST_END_YEAR = options.LR_END_YEAR;
    end;
    if ~isempty( options.LR_END_Q )
        EST_END_Q = options.LR_END_Q;
    end;
    if ~isempty( options.long_run_V )
        use_long_run_V = options.long_run_V;
    end;

end;


fprintf( '\n%s\n', 'Estimating long-run prior:' );

% The variables that are to be modelled in the VAR are denoted names. All
% named variables must match exactly entries in text for the selection to
% work. The usrdifflev variable determines whether the priors are centered
% persistent (1) or non-persistent (0) coeff. values.
Yraw = num;
Ynames = text(1,3:end); % note: first two 'name's are date-related; ignore


% Go through the variables included in the main VAR, in order, and find
% where (if at all) they appear in the long-run database. The result will
% be a trimmed database for estimation in this method, which includes the
% same variables as the main VAR in the same order.
index_selected_data = [];

for i = 1:numel( names )
    j = find( strcmp( names{i}, Ynames ) );
    index_selected_data = [ index_selected_data j ]; %#ok<AGROW>
end;

% only the selected data is retained in the Yraw matrix.
Yraw = Yraw( :, index_selected_data );

disp( 'Database information:' );
cal.setFrequency( FREQ ); 
cal.setEndYearMonth( END_YEAR, END_Q );
cal.setStartYearMonth( START_YEAR, START_Q );
cal.displayDateInfo();

%fprintf( '\n%s\n\n', ['VAR order p = ' num2str(p) ...
%            '; constant ' constant_names{constant+1} ] );

% get estimation indexes
t0 = cal.getIndexFromDate( EST_START_YEAR, EST_START_Q );
t1 = cal.getIndexFromDate( EST_END_YEAR, EST_END_Q );
% and reverse...
[ y0 q0 ] = cal.getDateFromIndex( t0 );
[ y1 q1 ] = cal.getDateFromIndex( t1 );
dispstr = [ num2str(y0) '-' num2str(q0) ' thru ' ...
                                   num2str(y1) '-' num2str(q1) ];
disp( 'Estimation sample:' );
fprintf( '%s\n\n', dispstr );

%% -------------------------DATA HANDLING----------------------------------
% Drop any data not used in estimation. Then reset the DateHandler object.
Yraw = Yraw( t0:t1, : );
cal.setEndYearMonth( y1, q1 );
cal.setStartYearMonth( y0, q0 );

% Outputs Y1 and X1 are the raw data, manipulated so that the full sample
% is available, after taking away observations to create lags and to
% account for possibly h-step ahead specifications such as y(t+h)=B y(t) +
[ Y, X, Tstar ] = DoRawDataManipulation( ...
                      Yraw, constant, p, ...
                      0, 0, 1, cal );

% First get Ordinary Least Squares (OLS) estimates [struct]
olsparams = GetOLSestimates( Y, X );
s.data = {Y, X};
%% -------------------------SETTING THE AR PRIOR---------------------------

% Find the positions of regression variables in the VAR
index_select = [];
% number of long-run variables
n_lr = numel( Ynames( index_selected_data ) ); % same as size( Y, 2 )

for i = 1:n_lr
    j = find( strcmp( Ynames( index_selected_data(i) ), names ) );
    index_select = [ index_select j ]; %#ok<AGROW>
end;

% Assign the priorparams AR matrix to a new variable A
A = s.A_prior(constant+1:end,:);

% note index_select needs to be a column vector
rows_select = kron( ones(p,1), index_select' ) + ...
    kron( tril( ones(p,p), -1 )*ones(p,1)*numel(names), ones( n_lr, 1 ) );

% Set the rows, column entries of the prior parameter matrix corresponding 
% to the long-run variables (ignoring the constant)
A(rows_select,index_select) = olsparams.A_OLS(constant+1:end,:);

% set the priorparams coefficient matrix to the new value
s.A_prior( constant+1:end, : ) = A;


%% -------------------SETTING THE COVARIANCE PRIOR-------------------------
if ( use_long_run_V )
    V = s.V_prior;
    [ K, M ] = size( s.A_prior );

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


    % V.2. Pick out coordinates corresponding to variables in the long run
    % macro block of the vectorized coefficient matrix vec(A).

    % gives all the rows of A where a long-run macro variable appears in the 
    % model equations, as specified in 'index_select' function argument
    lr_row_select = ...
          kron( ones(p,1), index_select' ) + ...
                M*kron( tril(ones(p,p),-1)*ones(p,1), ones(n_lr,1) );

    % gives all the rows of vec(A) where a long-run macro variable appears in a
    % long-run macro equation
    index_select_all = kron( index_select'-1, ones( n_lr*p,1) )*M*p + ...
          kron( ones(n_lr,1), lr_row_select ) ; % BUGBUG


    % V.3. Impose the OLS estimates of the long run covariance matrix of the 
    % parameters onto the submatrix of V that corresponds to the row and column 
    % indexes of the long-run macro variables block in vec(A). The OLS estimate
    % for the covariance of A is E[(A-A0)(A-A0)'] = inv(X'X)*X'SIGMA*X*inv(X'X)
    V_OLS = olsparams.V_OLS; 
    % Drop rows corresponding to constant
    if ( 1 == constant )
        v2 = 1:n_lr*(1+n_lr*p);
        % drop the rows corresponding to the constant (1, 2+M*p, 3+2*M*p,...)
        const_select = tril( ones(n_lr,n_lr) )*ones(n_lr,1) + ...
                               tril( ones(n_lr,n_lr), -1)*ones( n_lr,1)*n_lr*p;
        v2( const_select ) = []; % drop constants
        V_OLS = V_OLS( v2, v2 ) ;   % select only rows/cols
    end;

    % -------------------------------------------------------------------------
    % There are two ways one could impose prior information on the covariances
    % between parameters; one takes only variances into account, the other the
    % full covariance between parameters...

    % Option 1: use OLS estimates of parameter variances alone
    D1 = zeros( M*M*p, 1 );
    D1( index_select_all ) = index_select_all ;        
    % create a temporary matrix which has D1 on its diagonal
    D2 = diag( D1 );
    V( 0 ~= D2 ) = rel_tight*diag( V_OLS ); % one could adjust 'tightness' here

    % Option 2: use OLS estimates of parameter variances and covariances (n.b.
    % can cause problems in combination with variance estimates from panel)
    %V( index_select, index_select ) = V_OLS; % one could adjust 'tightness' here
    % -------------------------------------------------------------------------

    % Ensure the new V prior is not singular
    if ( rank( V ) ~= size( V, 1 ) )
        warning( 'GetLongRunPrior:CovarianceMatrix', ...
                 'The prior covariance matrix V may be singular.' );
    end;


    % V.4. Set the new value of V_prior in the priorparams structure, where we
    % again ignore rows corresponding to constant (where present) since we
    % have done nothing to those rows above.
    s.V_prior( v1, v1 ) = V;    
end;    
%% ---------------- SETTING THE PRIOR ERROR COVARIANCE --------------------
%
s.S_prior( index_select, index_select ) = olsparams.SIGMA_OLS; %BUGBUG divide by v_prior??

%% ----------------------------WE'RE DONE----------------------------------
fprintf( '\n%s\n', 'Done setting long-run prior.' );

end
