function A_out = padMatrixWithZeros( A_in, index )

    % The dimensions of the matrix to be padded is M*p x M
    [ M_p, M ] = size( A_in );
    
    % The number of 'lags' in A_in
    p = M_p/M;
    
    % The number of additional 'variables' specified in 'index'
    n = length( index );
    
    % non-zero rows/columns of A_out at lag 1
    non_zero_indexes = ~ismember( 1:(M+n), index );

    X = bsxfun(@and, non_zero_indexes, non_zero_indexes' ) + 0;

    A_out = [];
    
    for i = 1:p,
        TMP = zeros( M+n, M+n );
        TMP(X~=0) = A_in( 1+(i-1)*M:M*i, : );
        A_out = [ A_out; TMP ];
    end;
end