function [data_table, Ynames, cal] = GetVARData( data_dir, data_file, data_tab )
% Task to read data in from an Excel file, return the data table along with
% the variable names in the file (Ynames) and a DateHandler object for 
% later data manipulation.

% Use of xlsread() is deprecated, so now uses readtable()
data_table = readtable( [data_dir data_file], 'Sheet', data_tab );

% The column headings are a property of the data table
Ynames = data_table.Properties.VariableNames;

% We are looking for a cell array of strings giving the dates, in a column
% named "date". If that name doesn't exist, assume it's the first column
try
    x_date_str = data_table.date;
catch
    warning( 'Spreadsheet should contain a column named "date".' );
    x_date_str = data_table(:,1);
end

% Recover information on the dating used by interrogating the contents of
% x_date_str
if iscellstr( x_date_str )
    % If we were passed a cell array of strings, extract date information
    % from them. This includes undated data, if its dates are listed as
    % strings. If they are listed as numbers, we go to the else statement.
    format = getInputFormat( x_date_str );

    % edit calendar settings to match database (**not** estimation dates)
    START_YEAR = format.start_year; 
    START_Q = format.start_q;
    END_YEAR = format.end_year; 
    END_Q = format.end_q;
    FREQ = format.frequency;

else
    % If we were not passed a cell array of strings, then we were passed
    % numbers. In which case treat as 'undated'
    START_YEAR = x_date_str(1); 
    START_Q = 0; 
    END_YEAR = x_date_str(end); 
    END_Q = 0;
    FREQ = 1;
    
end

cal = DateHandler(); % instantiate DateHandler class


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
    
    % Non-empty entries in date_str return a logical true
    b_non_empty = cellfun( @(x) ~isempty(x), date_str );
    
    % It may be that the date strings are all empty. Assert they are not.
    err_msg = 'DataFormat::There are no date strings to read';
    assert( any( b_non_empty ), err_msg );

    % Find the first/last non-empty entries in the date_str. b_non_empty is
    % guaranteed to return at least one idx_non_empty since the assertion
    % was passed.
    idx_non_empty = find( b_non_empty );
    counter_first = idx_non_empty( idx_non_empty(1) );
    counter_last =  idx_non_empty( idx_non_empty(end) );

    % Let ss be the string date of the first observation
    ss = date_str{counter_first};
    
    if contains( ss, 'Q', 'IgnoreCase', 1 )
        % If there's a Q, we have 4 sub-periods
        format.frequency = 4;
        % The position of the Q
        index = strfind( ss, 'Q' );
        % The first characters are the period
        format.start_year = str2double( ss(1:index-1) );
        format.start_q = str2double( ss( index+1:end ) );
        
        % Counter is set to the last date
        ss = date_str{counter_last};
        format.end_year = str2double( ss(1:index-1) );
        format.end_q = str2double( ss( index+1:end ) );

        
    elseif contains( ss, 'M', 'IgnoreCase', 1 )
        % If there's an M, we have 12 sub-periods
        format.frequency = 12;
        % The position of the M
        index = strfind( ss, 'M' );
        % The first characters are the period
        format.start_year = str2double( ss(1:index-1) );
        format.start_q = str2double( ss( index+1:end ) );
        
        % Counter is set to the last date
        ss = date_str{counter};
        format.end_year = str2double( ss(1:index-1) );
        format.end_q = str2double( ss( index+1:end ) );

    else
        % Annual or undated
        format.frequency = 1;
        format.start_year = str2double( ss );
        format.end_year = str2double( date_str{counter} );
        format.start_q = 1;
        format.end_q = 1;

    end

end