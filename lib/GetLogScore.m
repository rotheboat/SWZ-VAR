% Given a vector of values containing log predictive likelihood values for a number of forecast
% periods h, computes the log score.
function score = GetLogScore( v_pred_likelihood, obs )

score = sum( log( v_pred_likelihood( obs ) ) );

end