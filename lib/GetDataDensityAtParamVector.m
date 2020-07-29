function log_data_density_at_point = GetDataDensityAtParamVector(mY,mX,Astar,Fstar)
% Given any parameter vectors (vbstar, vgstar), evaluate the data density
% (likelihood function) given data (mY, mX)
% Author: Roland Meeks
[ Tobs, M ] = size( mY );

% Note that the determinant of A is not necessarily positive 
log_data_density_at_point = -M*Tobs*log(2*pi)/2 + Tobs*log(abs(det(Astar))) ...
    - trace( ( mY*Astar - mX*Fstar )'*( mY*Astar - mX*Fstar ) )/2;

end