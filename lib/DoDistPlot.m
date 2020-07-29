% Plots the marginal predictive distribution for a single variable
function DoDistPlot( postmean, Pred_mean, names, index, varargin )

% check for presence of optional argument 'type', indicating the type of
% univariate distribution that should be plotted.
if ( 5 == nargin )
    type = varargin{1};
else
    type = '';
end;

% obtain location and scale for variable of interest
pscale = sqrt( postmean.SIGMA(index, index) );
ploc = Pred_mean(index);

% figure out the limits of the x-axis
xrange = 3*sqrt(pscale);
xstep = 2*xrange/100;
xvals = [-xrange:xstep:xrange] + ploc;

% check which case is supplied
if ( strcmp( type, 't' ) ) % the predictive dens it student's t
    
    df = postmean.v_post;
    yvals = locsctpdf( xvals, df, ploc, pscale );
    plot( xvals, yvals );
    xlabel(names(index(1)));
    
elseif ( strcmp( type, 'N' ) )  % the predictive dens is normal
    
    yvals = normpdf( xvals, ploc, pscale );
    plot( xvals, yvals );
    xlabel(names(index(1)));
    
else % use the kd estimator
    
end;