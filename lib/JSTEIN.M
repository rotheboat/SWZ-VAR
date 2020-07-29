function [jsest] = jstein(bhat)

% 	JStein() Procedure
% 
% 	A MATLAB procedure designed to compute the original
% 	James-Stein estimate of the linear model parameter Beta.
% 	The required inputs are:
% 		bhat		k-vector of estimated linear model parameters
% 	The information returned by the procedure is:
% 		jsest	k-vector of James-Stein parameter estimates
% 	The James-Stein estimator is computed as:
% 			jsest = [1- (K - 2)/ bhat'bhat]*bhat 
% 	and the number of unknown parameters must satisfy K > 2.


k = size(bhat,1);

if k < 3;	
    disp('The number of unknown parameters must exceed K = 2');
end
jsest = (1 - (k - 2)/(bhat'*bhat))*bhat;
	