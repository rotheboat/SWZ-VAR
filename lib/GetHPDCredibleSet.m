%% Returns the lower and upper bounds of the highest posterior density set
%  for every parameter in the specified equation.
%
%  postest = a struct containing posterior draws
%  field_name = a string containing the name of the field within postest
%       containing the relevant posterior draws
%  equation_ = an integer specifying which column of the relevant matrix is
%       of interest
%  alpha = fraction giving the size of the HPD set (1-alpha), typically
%       0.05 or 0.1
%
% @returns [ mode, min, max ] of histogram of postest.field_name
function hpd_set = GetHPDCredibleSet( postest, field_name, equation_, alpha )

bin_width = 50;

[ draws, K, M ] = size( postest.B_post );

if strcmp( field_name, 'A_post' )
    n_ = M;
else
    n_ = K;
end

hpd_set = zeros( n_, 3 );

mmB = postest.(field_name);
    for i = 1:n_
    temp_ = zeros( draws, 1 );
        for j = equation_
            [num val] = hist( squeeze(mmB(:,i,j) ), bin_width ); 
            [b_, i_] = max( num );
            
            temp_(1)=num( i_ )/draws;
            h = 2;
            val_min = val( i_ );
            val_max = val( i_ );
            while temp_ < 1 - alpha
                temp_(h) = temp_(h-1) + ( num( i_+h-1 )/draws ) + ...
                    ( num( i_-h+1 )/draws );
                val_max = val( i_+h-1 );
                val_min = val( i_-h+1 );
                
                h = h + 1;
            end
            hpd_set( i, : ) = [ val( i_ ), val_min, val_max ];
        end
    end

end