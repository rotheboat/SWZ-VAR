% A utility class that allows us to map between a calendar date in the form of yyyy-mm and an
% observation index. The format of yyyy is the base time period under consideration, which is
% usually a year; the format of mm is observations per base time period, for example, months.
% But one could index yyyy to mean weeks and mm to mean days. The only requirement is that the
% sub-index has to occur at least once per base time period. The frequency of the data is set
% to i_freq.
%
% Ported from Ox to Matlab by Roland Meeks, March 1, 2011
classdef DateHandler < handle

    properties ( Access=protected )
        
        start_year;
        start_month;
        end_year;
        end_month;
        i_freq = 12;
        
    end;
    
    methods 
        
        %------------------------------------------------------------------------------
        %@Author: RM
        %@Description: internal utility function
        %@Returns: matrix [year,month]
        %------------------------------------------------------------------------------
        function dates = getDates( obj )
        
            month = (1:1:obj.i_freq)';
            year  = (obj.start_year:1:obj.end_year)';
            dates = [kron( year, ones(obj.i_freq,1) ) kron( ones(length(year),1), month )];
            dates = dates(obj.start_month:end, :);
            
            % drop the observations in the final year, after the last month
            dates = dates( 1:end - (obj.i_freq - obj.end_month), :);
            
        end;
    
        %------------------------------------------------------------------------------
        %@Author: Roland Meeks
        %@Description: given the index of an observation, returns the date.
        %@Arguments: integer valued index
        %@Returns: matrix [year, month]
        %------------------------------------------------------------------------------
        function [ i_year, i_month ] = getDateFromIndex( obj, index )
        
            dates = getDates( obj );
            i_year = dates( index, 1 );
            i_month = dates( index, 2 );
            
        end;
        %------------------------------------------------------------------------------
        %@Author: Roland Meeks
        %@Description: given the year and month of an observation, returns its index.
        %@Arguments: int, int
        %@Returns: int	
        %------------------------------------------------------------------------------
        function index = getIndexFromDate( obj, pyear, pmonth )
        
            dates = getDates( obj );
            % What this piece does is create a boolean matrix which is unity when
            % both the year and month are "dot equal" to the corresponding values in
            % the array 'dates'. It then uses that boolean matrix to pick an element
            % out of the index matrix created by the 'cumulate' command.
            index = ( dates(:,1) == pyear & dates(:,2) == pmonth );
            
            allindex = 1:length( index );
            
            index = allindex(index);
        end;
        %------------------------------------------------------------------------------
        %@Author: RM
        %@Description: overrides the default time limits of the entire database
        %------------------------------------------------------------------------------
        function setStartYearMonth( obj, pyear, pmonth )
        
            obj.start_year = pyear;
            obj.start_month = pmonth;
            
            if ( or( obj.i_freq < obj.start_month, obj.i_freq < obj.end_month ) )
                warning( 'DateHandler:setStartYearMonth' , 'Frequency mismatch' );
            end;
            
            if ( obj.start_year > obj.end_year )
                warning('DateHandler:setStartYearMonth', 'Start year is set after end year');
            elseif (obj.start_year == obj.end_year )
                if (obj.start_month > obj.end_month )
                    warning('DateHandler:setStartYearMonth', 'Start month is set after end month');
                end;
            end;
        
        end;
        %------------------------------------------------------------------------------
        %@Author: RM
        %@Description: overrides the default time limits of the entire database
        %------------------------------------------------------------------------------
        function setEndYearMonth( obj, pyear, pmonth )
        
            obj.end_year = pyear;
            obj.end_month = pmonth;
            
            if ( or( obj.i_freq < obj.start_month, obj.i_freq < obj.end_month ) )
                warning( 'DateHandler:setStartYearMonth' , 'Frequency mismatch' );
            end;
            
        end;
        %------------------------------------------------------------------------------
        %@Author: RM
        %@Description: overrides the default frequency (observations per year) of the 
        %              entire database
        %------------------------------------------------------------------------------
        function setFrequency( obj, pfreq )
        
            % lowest acceptable frequency is annual, one obs. per year
            if ~( pfreq < 1 )
                obj.i_freq = pfreq;
                disp( ['Data frequency set to ' num2str(pfreq) ' obs. per year'] );
            else
                error( 'Frequency must be annual or higher' );
            end;
            
        end;
        
        %------------------------------------------------------------------------------
        %@Author: Roland Meeks
        %@Description: given the year and month of an observation, returns its index.
        %@Arguments: int, int
        %@Returns: int	
        %------------------------------------------------------------------------------
        function displayDateInfo( obj )
        
            disp( ['Dates run ' num2str(obj.start_year) '-' num2str(obj.start_month) ...
                   ' to ' num2str(obj.end_year) '-' num2str(obj.end_month) ] );
            
        end;
                
    end; % end public methods
end % end classdef