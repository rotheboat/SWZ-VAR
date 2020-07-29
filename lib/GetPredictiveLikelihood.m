function plik = GetPredictiveLikelihood( density, Xmesh, Ymesh, data, index )

% Want to find the coordinate in [Xmesh Ymesh] that is as close as possible
% to the data
Xcoord = min( data(index(1)) );
Ycoord = min( data(index(2)) );

% Then want to find the density estimate corresponding to that data point.
% This is the (marginal) PL for the estimated density
plik = density( Xcoord, Ycoord );

