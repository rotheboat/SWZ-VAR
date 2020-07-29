% This script is for producing plots of counterfactual variables
x_data = datenum( x_date_str(t0+p:t1) );
tick_interval = 16;

if plot_deviations_from_data
    cfdata_trans = cfdata - ( shiftdim( repmat( Y, 1, 1, counterfactual.num_draws ), 2 ) );
    
    for j = 1:numel( names )
        if ( ismember( names{j}, plot_integrated ) && ismember( names{j}, plot_names ) )
            for h=2:length(x_data)
                cfdata_trans(:,h,j) = cfdata(:,h,j) + cfdata(:,h-1,j);
            end;
        end;
    end;
else
    cfdata_trans = cfdata;
end;

% If the `path` option is selected, then for every draw we plot the
% resulting counterfactual path for the selected variables, with paths
% overlaid one top of each other.  Because this could be many paths, we
% select an evenly spaced 200 (if there are more than 200 draws).
if ( strcmpi( strtrim(plot_type), 'path' ) )

    figure(4); clf; 

    for draw = 1:100:counterfactual.num_draws %BUGBUG

        Y_new = squeeze( cfdata_trans( draw, :, : ) );

        counter = 1;

        for j=vars_to_plot

            subplot( plot_rows, plot_cols, counter )
            hold on;

            if ( ismember( names{j}, plot_in_diffs ) )
                D_Y_new = [ NaN(4,1); Y_new(5:T,j) - Y_new(1:T-4,j) ];
                plot( x_data, D_Y_new, 'b' );
            elseif ( ismember( names{j}, plot_integrated ) )
                % Convert variables that need it to % space
                if ( ismember( names{j}, 'DWTRIGGER' ) )
                    I_Y_new = 8.5 + cumsum( Y_new(:,j) ); % BUGBUG 1990Q3 initval hard coded
                elseif ( ismember( names{j}, 'DWPUBRAR' ) )
                    I_Y_new = 9.73 + cumsum( Y_new(:,j) ); % BUGBUG 1990Q3 initval hard coded
                else
                    I_Y_new = 100*(1 + cumsum( Y_new(:,j) ) );
                end
                plot( x_data, I_Y_new, 'b' );            
            else
                % Convert variables that need it to % space
                if ( ismember( names{j}, {'M4LHHSEC','M4LPNFC'} ) )
                    plot( x_data, 100*Y_new(:,j), 'b' );
                else
                    plot( x_data, Y_new(:,j), 'b' );
                end
            end

            axis tight;        
            counter = counter+1;
            hold off;
        end

    end

% If the `quantile` option is selected, then for every time period and
% every selected variable we compute the percentiles of the distribution of
% paths, and plot those.
elseif ( strcmpi( strtrim(plot_type), 'quantile' ) )
    
    figure(5); clf; 
    
    for j = vars_to_plot,
        % we compute the pointwise percentiles at each date
        if ( ismember( names{j}, plot_integrated ) )
            cfdata_trans(:,:,j) = (1 + cumsum( cfdata_trans(:,:,j),2 ));
        end
        
        for t = 1:T,
            % if the variable in question is to be plotted in differences
            % rather than levels, we need to compute the percentiles of the
            % differenced data, which are not the same as the differences
            % in percentiles!
            if ( ismember( names{j}, plot_in_diffs ) )
                % Because we take a 4 quarter difference, we pad the first
                % 4 periods with NaNs (which don't results in any plot)
                if ge( t, 5 )
                    qu_out = ...
                        quantile( cfdata_trans(:,t,j) - cfdata_trans(:,t-4,j), quantiles_to_plot );
                else
                    qu_out = NaN( size( quantiles_to_plot ) );
                end
            else
                % For variables in levels, we can just calculate the
                % percentiles from the data array directly, without
                % differencing 
                qu_out = ...
                   quantile( cfdata_trans( :, t, j ), quantiles_to_plot );
            end

            % Convert variables that need it to % space
            if ( ismember( names{j}, {'M4LHHSEC','M4LPNFC'} ) )
                qu_out = 100*qu_out;
            elseif ( ismember( names{j}, {'LHPALL'} ) && plot_deviations_from_data )
                qu_out = 100*qu_out;
            end
            
            % still inside the time loop, store the (three) quantiles for
            % the variables in question at time t
            path_qu16(t,j) = qu_out( 1 );
            path_qu50(t,j) = qu_out( 2 );
            path_qu84(t,j) = qu_out( 3 );
        end;
    end;

    counter = 1;
    
    for j=vars_to_plot

        subplot( plot_rows, plot_cols, counter )
        
        hold on;
        if ( ismember( names{j}, plot_in_diffs ) )
           % create an area plot that stacks the 16/84 quantiles, then
           % colors the 16 quantile white and the 84 grey, then finally
           % makes sure the figure axes are on top...
           ha = area(x_data(5:T),[path_qu16(5:end,j), ...
                path_qu84(5:end,j)-path_qu16(5:end,j)], ...
                min(path_qu16(5:end,j)),'LineStyle','none');
           plot( x_data(5:T), path_qu50(5:end,j), 'k', 'LineWidth', 2.0 ); 
        else
            ha = area(x_data,[path_qu16(:,j), path_qu84(:,j)-path_qu16(:,j)], ...
                min(path_qu16(:,j)),'LineStyle','none');
            plot( x_data, path_qu50(:,j), 'k', 'LineWidth', 2.0 );
            
        end
        
        set(ha(1),'FaceColor',[1,1,1]);
        set(ha(2),'FaceColor',[0.8,0.8,0.8]);
        set(gca,'layer','top');

        counter = counter+1;
        hold off;
        
    end    
    
else
    warning( 'Option `plot_type` should be `path` or `quantile`.' );
end

% Plot the data
counter = 1;

for j=vars_to_plot
    subplot( plot_rows, plot_cols, counter );
    hold on;

    if plot_deviations_from_data
        plot( x_data, zeros(1,length(x_data)), '-g' );
        title(long_names{j},'Interpreter','latex','FontSize',16);
    else
        if ( ismember( names{j}, plot_in_diffs ) )
            D_Y = [ NaN(4,1); Y(5:T,j) - Y(1:T-4,j) ];
            plot( x_data, D_Y, 'm', 'LineWidth', 2.0 );
            title( strcat( 'D.', names{j}));
        elseif ( ismember( names{j}, plot_integrated ) )
            if ( ismember( names{j}, 'DWTRIGGER' ) )
                I_Y = 8.5 + cumsum( Y(:,j) ); % BUGBUG 1990Q3 initval hard coded
            elseif ( ismember( names{j}, 'DWPUBRAR' ) )
                I_Y = 9.73 + cumsum( Y(:,j) ); % BUGBUG 1990Q3 initval hard coded
            else
                I_Y = 100*(1 + cumsum( Y(:,j) ) );
            end
            plot( x_data, I_Y, 'm', 'LineWidth', 2.0 );
            title( strcat( 'I.', names{j}));
        else
            % Convert variables that need it to % space
            if ( ismember( names{j}, {'M4LHHSEC','M4LPNFC'} ) )
                plot( x_data, 100*Y(:,j), 'm', 'LineWidth', 2.0 );
            else
                plot( x_data, Y(:,j), 'm', 'LineWidth', 2.0 );    
            end

            title(long_names{j},'Interpreter','latex','FontSize',16);
        end
    end
    
    hold off;
    counter = counter+1;
        
    % make a date axis
    set(gca,'XTick',x_data([1:tick_interval:T]));
    datetick( 'x', 'yy-QQ', 'keepticks' );
    %AxisHandle = get(gcf,'CurrentAxes');
    set(gca,'FontSize',12);
    axis tight;    
end

clear( 'D_Y', 'D_Y_new', 'I_Y' );

figName = 'counter-factual';
banks_format_plot;  % format and save the figure