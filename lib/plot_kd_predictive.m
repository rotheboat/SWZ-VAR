function plot_kd_predictive( Y_pred, index )

% plot predictive density
figure; M = 3;

for i=1:M,
    subplot(M,1,i);
    ksdensity( Y_pred(:,i) );
end;

% Plot the scatters
figure(2)


%contour(Xmesh,Ymesh,density)

for i = 1:M,
    for j = 1:M,
        if (j > i),
            subplot( M-1, M-1, (j-1)+(i-1)*(M-1) );
            
            data = [Y_pred(:,i),Y_pred(:,j)];
            MAX = max(data,[],1); 
            MIN = min(data,[],1); 
            Range = MAX-MIN;
            MAX_XY = MAX+Range/4; 
            MIN_XY = MIN-Range/4;
            [bandwidth,density,Xmesh,Ymesh]=kde2d(data,2^7,MIN_XY,MAX_XY);
            %surf(Xmesh,Ymesh,density,'LineStyle','none')
            contour( Xmesh,Ymesh,density );
        end;
    end;
end;

colormap autumn;

disp( 'Done.' );