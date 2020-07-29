% Plot the impulse-responses
for k = 1:M,
    for j = 1:M,
        subplot( M, M, j+(i-1)*M );

        plot( 1:irf_horizon, squeeze(imp_resp( j, k, : )), ...
                1:irf_horizon, zeros(1,irf_horizon), 'g' );
        title( names(j) );
        
        %plot( 1:ihor, squeeze(imp_resp(i,j,:)) );
    end;
end;