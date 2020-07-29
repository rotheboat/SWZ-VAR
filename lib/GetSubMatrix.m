% Given an M-dimensional square matrix X, returns a sub-matrix
% corresponding to row/column indices given in index. Makes most sense to
% use the function if X is symmetric. Note index is a vector.
function submat = GetSubMatrix( X, index )

% find the dimension of the VAR
[M N] = size( X );

if (M ~= N)
    error( 'GetSubMatrix: function requires a square matrix' );
end;

% need an increasing sequence of indexes to correctly pick out submatrix
index = sort( index );

% create a set of coordinates
lwn = length( index );

% this creates a set of lwn^2 coordinates of the required submatrix
vcoords = [ kron(ones(1,lwn),index)', kron(index,ones(1,lwn))' ];
% this small function converts a coordinate to a position in the vectorized
% matrix y = vec(Y)
f = @( y,M ) (y(:,2)-1)*M + y(:,1);

vecX = X(:);
submat = reshape( vecX( f(vcoords,M) ), lwn, lwn );

% RM There must be a neater way to do this; you'd think that Matlab would
% let you just go matrix X[ set of coordinates ]. Meh.