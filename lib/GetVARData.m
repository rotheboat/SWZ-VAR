function [Yraw, Ynames, cal] = GetVARData( data_dir, data_file, data_tab )

% BUGBUG catch errors
[num, text, ~] = xlsread( [data_dir data_file], data_tab );

% The variables that are to be modelled in the VAR are denoted names. All
% named variables must match exactly entries in Ynames for the selection to
% work. The usrdifflev variable determines whether the priors are centered
% persistent (1) or non-persistent (0) coeff. values.
Yraw = num;
Ynames = text(1,2:end); % note: first 'name' is 'qdate's; ignore
x_date_str = text(2:end,1); % Here we have formatted dates for plotting

% In any case, name the data you load 'Yraw', in order to avoid changing the
% rest of the code. Note that 'Yraw' is a matrix with T rows by M columns,
% where T is the number of time series observations (usually months or
% quarters), while M is the number of VAR dependent macro variables.
format = getInputFormat( x_date_str );

cal = DateHandler(); % instantiate DateHandler class

% edit calendar settings to match database (**not** estimation dates)
START_YEAR = format.start_year; 
START_Q = format.start_q; % DO NOT EDIT
END_YEAR = format.end_year; 
END_Q = format.end_q;     % UNLESS EXCEL CHGD
FREQ = format.frequency;

% .........................................................................
disp( 'Database information:' );
cal.setFrequency( FREQ ); 
cal.setEndYearMonth( END_YEAR, END_Q );
cal.setStartYearMonth( START_YEAR, START_Q );
cal.displayDateInfo();



end

function format = getInputFormat( date_str )
% Utility function that processes text information read in from Excel file,
% so that DateHandler can know the range of the database

    % Get rid of special characters in date_str
    date_str = erase( date_str, '-' );
    date_str = erase( date_str, ':' );
    date_str = erase( date_str, '_' );
    date_str = erase( date_str, ' ' );
    
    % Find the first/last non-empty entries in the date_str
    counter = 1;
    
    % Find the first non-empty date string
    while and( isempty( date_str{counter} ), counter < length( date_str ) ) 
        counter = counter+1;
    end
    
    % It may be that the date strings are all empty
    err_msg = 'DataFormat::There are no date strings to read';
    assert( ~eq( counter, length( date_str ) ), err_msg );
    
    ss = date_str{counter};
    
    % Find the last non-empty date string
    while and( ~isempty( date_str{counter} ), counter < length( date_str ) )
        counter = counter+1;
    end

    if contains( ss, 'Q', 'IgnoreCase', 1 )
        % If there's a Q, we have 4 sub-periods
        format.frequency = 4;
        % The position of the Q
        index = strfind( ss, 'Q' );
        % The first characters are the period
        format.start_year = str2num( ss(1:index-1) );
        format.start_q = str2num( ss( index+1:end ) );
        
        % Counter is set to the last date
        ss = date_str{counter};
        format.end_year = str2num( ss(1:index-1) );
        format.end_q = str2num( ss( index+1:end ) );

        
    elseif contains( ss, 'M', 'IgnoreCase', 1 )
        % If there's an M, we have 12 sub-periods
        format.frequency = 12;
        % The position of the M
        index = strfind( ss, 'M' );
        % The first characters are the period
        format.start_year = str2num( ss(1:index-1) );
        format.start_q = str2num( ss( index+1:end ) );
        
        % Counter is set to the last date
        ss = date_str{counter};
        format.end_year = str2num( ss(1:index-1) );
        format.end_q = str2num( ss( index+1:end ) );

    else
        % Annual or undated
        format.frequency = 1;
        format.start_year = str2num( ss );
        format.end_year = str2num( date_str{counter} );
        format.start_q = 1;
        format.end_q = 1;

    end

end