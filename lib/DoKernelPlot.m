% Uses the function kde2d to compute a bivariate KD plot given a set of n
% pairs of input data.
%
% Roland Meeks 2010-11-02
function [bandwidth,density,Xmesh,Ymesh] = ...
                                    DoKernelPlot( Y_pred, w_names, index )

data = [Y_pred(:,index(1)),Y_pred(:,index(2))];
MAX = max(data,[],1); 
MIN = min(data,[],1); 
Range = MAX-MIN;
% playing with the range might be a good idea - default was 4
MAX_XY = MAX+Range/10; 
MIN_XY = MIN-Range/10;
[bandwidth,density,Xmesh,Ymesh] = kde2d(data,2^7,MIN_XY,MAX_XY);
contour( Xmesh,Ymesh,density );

xlabel(w_names(index(1)));
ylabel(w_names(index(2)));

colormap( 'copper' );

