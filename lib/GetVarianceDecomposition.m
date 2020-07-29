function mVarDec = GetVarianceDecomposition( postest )

[ draws, M, M, horizon ] = size( postest.IRF );
mVarDec = zeros( size( postest.IRF ) );

for d = 1:draws
    irf = squeeze( postest.IRF( d, :, :, : ) );
    
    total_variance = diag( irf( :, :, 1 )*irf( :, :, 1 )' );
    
    for i = 1:M
        mVarDec( d, :, i, 1 ) = diag( irf( :, i, 1 )*irf( :, i, 1 )' )./total_variance;
    end
    
    for h = 2:horizon
        total_variance = diag( irf( :, :, h )*irf( :, :, 1 )*irf( :, :, 1 )'*irf( :, :, h )' );
        
        for i = 1:M
            mVarDec( d, :, i, h ) = ...
                diag( irf( :, :, h )*irf( :, i, 1 )*irf( :, i, 1 )'*irf( :, :, h )' )./total_variance;
        end
    end
end

% END function
end

