% Name: xcorrAR.m
% Package: bmc
% From: mccallum\xcorrAnalytic.m
% Author: Roland Meeks
% Date: 28-Mar-2004
% Description: Calculate the model xcorr function using the formula given
% in Hamilton S10.2. Finds the MA representation by expanding the AR
% polynomial
% Inputs: The AR matrix of the VAR written in companion form LAMBDA.
%         The covariance matrix of the VAR OMEGA.
% Outputs: a 3-D matrix containing the covariance matrix GAMMA_s at
% xcov(:,:,s) and the correlation matrix inv(D0)*GAMMA_s*inv(D0) at corr(:,:,s).
% Notes: RM checked Mar-2006; see 'anotherXcorrTest.m'
function xcorr = xcorrAR(LAMBDA, OMEGA, num_leads_lags, iterations)

dimvar = size(LAMBDA,2);

GAMMA = zeros( dimvar, dimvar ); % initialise GAMMA as zeros #BUGBUG
PSI = eye( dimvar, dimvar ); % initialise PSI as identity

% calculate covariance mx for s = 0
for j = 1:(iterations - num_leads_lags),
    % GAMMA in Hamilton's notation is the s^th autocovariance matrix.
    % It is equal to the infinite sum of the PSI*OMEGA*PSI' matrices. In
    % our case, the expanded polynomial takes a simple form, so we can
    % easily iterate. Remember that GAMMA and PSI need to be properly
    % initialised to avoid errors, however.
    GAMMA = GAMMA + PSI*OMEGA*PSI';
    % PSI in Hamilton's notation is the MA matrix; in our case the
    % expansion of the inverse lag polynomial means that PSI_v is
    % LAMBDA^v. Since indexing in Matlab starts at 1, PSI_v = 
    % LAMBDA^(j-1).
    PSI = LAMBDA*PSI;
end; % end for

% the covariance matrix (s = 0) is stored in the num_leads_lags + 1
% position of the cross covariance matrix, because indexing starts at 1.
xcov(:,:, num_leads_lags + 1) = GAMMA;

% Calculation of the cross correlation matrices is simplified by the
% relationship GAMMA_s = LAMBDA*GAMMA_{s-1}
for s = 2:num_leads_lags + 1,
    xcov(:,:,s + num_leads_lags) = LAMBDA^(s-1)*xcov( :, :, num_leads_lags + 1 );
end; % end for

% For s < 0, we have that GAMMA_s = GAMMA'_{-s}.
for s = 1:num_leads_lags,
    % You can check this with the command, say all( xcorr_(-2) == xcorr_(2) ).
    xcov(:,:,s) = xcov(:,:,( (2*num_leads_lags + 1) - (s-1) ) )'; % for s < 0, the corr mx is the transpose 
                                                       % of the s^th GAMMA matrix
end; % end for

% the index of the zero^th cross covariance is num_leads_lags + 1, so the
% variances are the diagonal elements of this matrix, which we take the
% square root of to get the standard deviations
D0 = xcov( :, :, num_leads_lags + 1 );
D0 = diag ( sqrt( diag( D0 ) ) );

% the autocorrelation function is then
for j = 1:(2*num_leads_lags + 1),
    xcorr(:,:,j) = inv(D0)*xcov(:,:,j)*inv(D0);
end;