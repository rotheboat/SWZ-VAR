% DeleteMeClass

classdef DeleteMeClass < handle
    
    properties (Access=protected)
        xval;
        fval;
    end;
    
    methods
    
        function set_xval( obj, a )
            obj.xval = a;
        end;
        
        function x = get_xval( obj )
            x = obj.xval;
        end;
        
        function z = divideByY( obj, y )
            x = get_xval( obj ); % use a function
            z = x/y;
        end;
        
    end;
    
end