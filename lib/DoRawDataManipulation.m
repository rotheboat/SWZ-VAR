% Outputs Y1 and X1 are the raw data, manipulated so that the full sample
% is available, after taking away observations to create lags and to
% account for possibly h-step ahead specifications such as y(t+h)=B y(t) +
% ... The outputs Y and X are truncated such that a single forecast step
% can be undertaken, given the forecast_method.
function [ Y1, X1, Tstar, options ] = DoRawDataManipulation( Yraw, Zraw, options, cal )

% Get the VAR order, check it makes sense
if isfield( options, 'nlags' )
    p = options.nlags;
    assert( and( isint( p ), 0 < p ), 'Option `nlags` must be integer > 0' );
else
    % If there is no 'nlags' field in options struct, estimate a VAR(1)
    p = 1;
    options.nlags = p;
    warning( 'VAROptions::No lag length specified; continuing with VAR(1)' );
end
    
% Get initial dimensions of dependent variable
[Traw, ~] = size(Yraw);

% Get the number of exogenous variables, if any
[~, Nexog] = size(Zraw);

% Create any needed lags etc. where exogenous variables are present
if 0 < Nexog
    options.has_exo_vars = 1;
    options.num_exo_vars = Nexog;
    
    min_exo_lag = options.exo_lags(1);
    max_exo_lag = options.exo_lags(2);
    
    lag_cutoff = max( max_exo_lag, p );
    
else
    options.has_exo_vars = 0;
    options.num_exo_vars = 0;
    
    lag_cutoff = p;
end

% The model specification is different when implementing direct forecasts,
% compared to the specification when computing iterated forecasts.
if isfield( options, 'h_oos' )
    assert( options.h_oos >= 0, 'You must set h >= 0')
end    

% Now create VAR specification according to forecast method

% Direct forecasts require h-step-ahead VAR set-up
if ( 0 == options.forecast_method && 1 == options.forecasting )
    Y1 = Yraw(h+lag_cutoff:end,:);
    Y2 = Yraw(1:end-(h-1),:);
        
    % Generate lagged Y matrix. This will be part of the X matrix
    Ylag = mlag2(Y2,lag_cutoff); % Y2 is [T-(h-1) x M]. 
    Ylag = Ylag(1+lag_cutoff:end,:); % Ylag is [T-(h-1)-p x (Mp)]
        
    % The total number of available obs for case::direct
    Tstar = Traw - (h - 1) - lag_cutoff;
    
    % retrieve the date corresponding to the start index in the original data
    [ yr, qtr ] = cal.getDateFromIndex( h + lag_cutoff );
    % set the database dates to reflect the truncated data (after lags etc. accounted for)
    cal.setStartYearMonth( yr, qtr );
    disp('Done data manipulation.');
    cal.displayDateInfo();
        
% In the case of iterated forecasts, or no forecasts at all, the VAR is set
% up in one-step-ahead mode.
elseif ( 1 == options.forecast_method || 0 == options.forecasting ) 
    Y1 = Yraw(1+lag_cutoff:end,:);   % obviously agrees with abv for h==1 but for
    Y2 = Yraw(1:end,:);     % general h>1 does not
        
    % Generate lagged Y matrix. This will be part of the X matrix
    Ylag = mlag2(Y2,lag_cutoff); % Y2 is [T x M]. 
    Ylag = Ylag(1+lag_cutoff:end,:); % Ylag is [T-p x (Mp)]
        
    % The total number of available obs for case::iterated
    Tstar = Traw-lag_cutoff;

    % retrieve the date corresponding to the start index in the original data
    [ yr, qtr ] = cal.getDateFromIndex( 1 + lag_cutoff );
    % set the database dates to reflect the truncated data (after lags etc. accounted for)
    cal.setStartYearMonth( yr, qtr );    
    disp('Done data manipulation.');
    cal.displayDateInfo();
    
else
    error('Wrong choice of forecast_method');
end
    
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Create any needed lags etc. where exogenous variables are present
if options.has_exo_vars
    
    % If max_exo_lag = 0, then mlag2 returns an empty matrix.
    Zlag = mlag2( Zraw, max_exo_lag );
    
    % Start by creating a matrix [ Z(t), Z(t-1), ..., Z(t-max_exo_lag) ]
    Z1 = cat( 2, Zraw(1+lag_cutoff:end,:), Zlag( 1+lag_cutoff:end, : ) );
    
    % Truncate the full matrix of exogenous variables from the left if the
    % first required lag is 1 or above to obtain [Z(t-min_exo_lag), ... ]
    Z1 = Z1(:, 1+Nexog*min_exo_lag:end );
       
else
    options.has_exo_vars = 0;
    options.num_exo_vars = 0;
    
end


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Now define matrix X which has all the R.H.S. variables (constant, lags of
% the dependent variable and exogenous regressors/dummies)
X1 = Ylag;

% If there are exogenous vars, all their lags go after the lagged Ys: 
% X1 = [ Ylag, Z1 ]
if options.has_exo_vars
    X1 = cat(2, X1, Z1);
end

% If there is a constant, it goes after all the other variables:
% X1 = [ Ylag, Z1, 1 ]
if (0 < options.constant)
    X1 = cat(2, X1, ones(Tstar,1)); % #BUGBUG constant moved *after* lags
end

%% Check the selected data contains no nan values
% Note that it is possible that some exogenous variables have nan values
% that are truncated out of the sample once the Y data has been manipulated
% to create lags. This is not a bug.
err_msg = 'Data::Selected variables/data sample must not contain nan values.';
assert( all( all( ~isnan( X1 ) ) ), err_msg );
assert( all( all( ~isnan( Y1 ) ) ), err_msg );

end
% end function
