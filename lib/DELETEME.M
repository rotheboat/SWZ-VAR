% test DeleteMeClass

clear; clc;

h = DateHandler(); 

h.setEndYearMonth( 2003, 1 )
h.setStartYearMonth( 1999, 3 )
h.setFrequency( 4 );

h.getDates()

h.getIndexFromDate( 2000, 4 )

h.displayDateInfo();

% you can pass a handle to a (anonymous) function and by default it is passed by reference
doStuffToObj = @(obj) obj.setEndYearMonth( 2004, 1 );

doStuffToObj(h);

h.displayDateInfo();