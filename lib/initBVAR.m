% Sets random number generators either to a fixed seed > 0 or to be drawn 
% afresh, each time the code is run.
function initBVAR( seed )

if ( 0 == seed )
    randn('state',sum(100*clock)); %#ok<*RAND>
    rand('twister',sum(100*clock));
else
    randn('state',seed);
    rand('twister',seed);
end
