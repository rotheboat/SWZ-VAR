% Linear restrictions
i_r_L = find( strcmpi( 'LENDRATE', options.names ) );
i_r_D = find( strcmpi( 'DEP3M', options.names ) );
i_Q_L = find( strcmpi( 'LPRIVCRED', options.names ) );

M = numel(options.names);
P = options.nlags;
N = numel( options.exo_names );
if 0 < N
    K = M*P + N + N*(options.exo_lags(2) - options.exo_lags(1)) + options.constant;
else
    K = M*P + options.constant;
end

% Think of the rows of Q, R representing individual restrictions. The
% columns represent the variables. We initialize everything to be
% unrestricted in the first instance.

%% A matrix
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Equation: cost of debt
i_rst = 0;
Q( :, :, i_r_D ) = zeros( M, M );

% i_r_L = 0
i_rst = i_rst+1;
Q( i_rst, i_r_L, i_r_D ) = 1;

% i_Q_L = 0
i_rst = i_rst+1;
Q( i_rst, i_Q_L, i_r_D ) = 1;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Equation: markup
i_rst = 0;
Q( :, :, i_r_L ) = zeros( M, M );

% i_Q_L = 0
i_rst = i_rst+1;
Q( i_rst, i_Q_L, i_r_L ) = 1;

% i_r_D + i_r_L = 0
% i_rst = i_rst+1;
% Q( i_rst, i_r_L, i_r_L ) = 1;
% Q( i_rst, i_r_D, i_r_L ) = 1;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Equation: quantity of loans
i_rst = 0;
Q( :, :, i_Q_L ) = zeros( M, M );

% r_D = 0
i_rst = i_rst+1;
Q( i_rst, i_r_D, i_Q_L ) = 1;

% Store the Q matrix in the options structure
options.restrictions.linear_restrictions.Q = Q;

%% F matrix
i_CAR = M*P + find( strcmpi( 'CAR_WTSUM', options.exo_names ) );
i_CARMIN = M*P + find( strcmpi( 'CAR_MIN_WTSUM', options.exo_names ) );
i_NPL = M*P + find( strcmpi( 'NPLR', options.exo_names ) );
i_RFR = M*P + find( strcmpi( 'POLRT', options.exo_names ) );
i_GDP = M*P + find( strcmpi( 'LRGDP', options.exo_names ) );
i_PRICE = M*P + find( strcmpi( 'LPRICES', options.exo_names ) );

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Equation: cost of debt
i_rst = 0;
R( :, :, i_r_D ) = zeros( K, K );

% Q_L = 0 and i_L = 0 at all lags
for lag = 1:P
    i_rst = i_rst+1;
    R( i_rst, i_Q_L+(lag-1)*M, i_r_D ) = 1;
    i_rst = i_rst+1;
    R( i_rst, i_r_L+(lag-1)*M, i_r_D ) = 1;
end

% CAR + CARMIN = 0
i_rst = i_rst+1;
R( i_rst, i_CAR, i_r_D ) = 1;
R( i_rst, i_CARMIN, i_r_D ) = 1;

% GDP = 0
i_rst = i_rst+1;
R( i_rst, i_GDP, i_r_D ) = 1;

% PRICE = 0
i_rst = i_rst+1;
R( i_rst, i_PRICE, i_r_D ) = 1;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Equation: markup
i_rst = 0;
R( :, :, i_r_L ) = zeros( K, K );

% Q_L = 0 amd r_L + r_D = 0 at all lags
for lag = 1:P
    i_rst = i_rst+1;
    R( i_rst, i_Q_L+(lag-1)*M, i_r_L ) = 1;
    
%     i_rst = i_rst+1;
%     R( i_rst, i_r_L+(lag-1)*M, i_r_L ) = 1;
%     R( i_rst, i_r_D+(lag-1)*M, i_r_L ) = 1;
end

% CAR + CARMIN = 0
i_rst = i_rst+1;
R( i_rst, i_CAR, i_r_L ) = 1;
R( i_rst, i_CARMIN, i_r_L ) = 1;

% GDP = 0
i_rst = i_rst+1;
R( i_rst, i_GDP, i_r_L ) = 1;

% PRICE = 0
i_rst = i_rst+1;
R( i_rst, i_PRICE, i_r_L ) = 1;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Equation: quantity of loans
i_rst = 0;
R( :, :, i_Q_L ) = zeros( K, K );

% r_D = 0 at all lags
for lag = 1:P
    i_rst = i_rst+1;
    R( i_rst, i_r_D+(lag-1)*M, i_Q_L ) = 1;
end

% CAR = 0
i_rst = i_rst+1;
R( i_rst, i_CAR, i_Q_L ) = 1;

% CARMIN = 0
i_rst = i_rst+1;
R( i_rst, i_CARMIN, i_Q_L ) = 1;

% RFR = 0
i_rst = i_rst+1;
R( i_rst, i_RFR, i_Q_L ) = 1;

% NPLR = 0
i_rst = i_rst+1;
R( i_rst, i_NPL, i_Q_L ) = 1;


% Store the R matrix in the options structure
options.restrictions.linear_restrictions.R = R;
