function [Yraw, Zraw, cal] = GetEstimationSubsample( Yraw, Zraw, cal, syr, sqtr, eyr, eqtr )

    % get estimation indexes
    t0 = cal.getIndexFromDate( syr, sqtr );
    t1 = cal.getIndexFromDate( eyr, eqtr );
    % and reverse...
    [ y0, q0 ] = cal.getDateFromIndex( t0 );
    [ y1, q1 ] = cal.getDateFromIndex( t1 );
    dispstr = [ num2str(y0) '-' num2str(q0) ' thru ' ...
                                       num2str(y1) '-' num2str(q1) ];
    disp( 'Estimation sample:' );
    fprintf( '%s\n\n', dispstr );


    %% --------------------------DATA HANDLING---------------------------------
    % Drop any data not used in estimation. Then reset the DateHandler object.
    Yraw = Yraw( t0:t1, : );
    if 0 < size( Zraw, 2 )
        Zraw = Zraw( t0:t1, : );
    end
    cal.setEndYearMonth( y1, q1 );
    cal.setStartYearMonth( y0, q0 );

end