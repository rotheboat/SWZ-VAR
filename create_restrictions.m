% Create restriction matrices
options.restrictions = cell(1,2);

M = numel(options.names);

% F matrix
mRF = zeros( M*options.nlags + numel(options.exo_names) ... 
    + options.constant, M );

% Always include constant
mRF( end, : ) = 1;

mRF( 1, 1 ) = 1;
mRF( 1, 2 ) = 1;
mRF( 2, 2 ) = 1;
% mRF( 3, 2 ) = 1;
mRF( 2, 3 ) = 1;
mRF( 3, 3 ) = 1;

for j = 2:options.nlags
    mRF( (j-1)*M+1:j*M, : ) = ...
        mRF( 1:M, 1:M );
end

% CAR
mRF( options.nlags*M+1, 1 ) = 1;
mRF( options.nlags*M+1, 2 ) = 1;
% CAR MIN
mRF( options.nlags*M+2, 1 ) = 1;
mRF( options.nlags*M+2, 2 ) = 1;
% LLR
mRF( options.nlags*M+3, 2 ) = 1;
% 'LRGDP'
mRF( options.nlags*M+4, 3 ) = 1;
% 'LPRICES'
mRF( options.nlags*M+5, 3 ) = 1;
% 'POLRT'
mRF( options.nlags*M+6, 1:2 ) = 1;
% % 'TREND' (if present)
% mRF( options.nlags*M+7, 3 ) = 1;

options.restrictions{2} = mRF;

% A matrix 
mRA = zeros(3,3);
% Funding cost
mRA(1,1) = 1;
% Markup
mRA(1,2) = 1;
mRA(2,2) = 1;
% mRA(3,2) = 1;
% Credit demand
mRA(2,3) = 1;
mRA(3,3) = 1;
options.restrictions{1} = mRA;