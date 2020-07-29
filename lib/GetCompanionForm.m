% Returns the companion form of the autoregressive component of a VAR(p).
%
% @args:
%   mB = K x M matrix where M = number of variables in Y, K = M*p+constant
%
% @author:  Roland Meeks
% @date:    June, 2014
% @modified: 
%   July, 2020: to account for the possible presence of exogenous variables
%   in the VAR.
function mC = GetCompanionForm( mB, options )

[ ~, M ] = size( mB );

% Number of lags of the endogenous variables
p = options.nlags;

mC = [ mB( 1:M*p, : )'; eye( M*(p-1) ), zeros( M*(p-1), M ) ];

end