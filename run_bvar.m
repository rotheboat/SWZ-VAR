clear
%close all;
%clf;
%clc;
disp( 'Indonesia Data' );
addpath( '.\lib\' );
initBVAR( 20081118 ); % arg ~= 0 gives rand seed
tic;

%% ---------------------------USER SETTINGS--------------------------------
%
% Define specification of the VAR model
options.constant = 1; % 1: if you desire intercepts, 0: otherwise 
options.nlags = 1; % Number of lags on dependent variables

% Dates used for modelling (may be the same as the database dates)
EST_START_YEAR = 2001; EST_START_Q = 4;
EST_END_YEAR = 2018; EST_END_Q = 3;

% Set prior type for BVAR model:
prior_name = 'Litterman';
% Set prior parameters lambda (or omit, and take defaults)
% lambda(1):
% lambda(2): smaller values squeeze towards random walk
% lambda(3): larger values squeeze lags > 1 towards zero
% lambda(4): smaller values shrink cofficients on exogenous terms to zero
options.lambda = [ 1, .5, 2, 1, 0, 0, 0 ]; % diffuse
% options.lambda = [ 0.6, 0.15, 1, 10, 0, 0, 0 ]; % loose
% options.lambda = [ 0.6, 0.1, 1, 0.1, 5, 5, 0 ]; % tight

              
% MODEL DEFINITION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_dir = 'C:\Users\Home\Box Sync\Macropru_Counterfactuals\data\';
data_file = 'indonesia.xlsx';
data_tab = 'export';
% data_dir = 'C:\Users\rmeeks\Box Sync\Macropru_Counterfactuals\data\';
% data_file = 'indonesia.xlsx';
% data_tab = 'export';

% Variable names to include in the empirical model (case sensitive)
options.names = {'DEP3M', 'LENDRATE', 'LPRIVCRED'};
% Variables names to include as exogenous covariates
options.exo_names = {'CAR_WTSUM', 'CAR_MIN_WTSUM', 'NPLR', 'LRGDP', 'LPRICES', 'POLRT'};
% Specify as a pair [minimum lag, maximum lag] (can be [0,0])
options.exo_lags = [1,1];
%'Y','P','PCOM','FF','NBR','TR'
options.long_names = {'Deposit Rate', 'Loan Rate', 'Credit'};
units = {'ppt','ppt','%','%','%','%','%'};
rescale = { 1, 1, 1, 1 };
               
% SIMULATION SETTINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                   
options.number_independent_chains = 50;
options.draws = 100;
options.burn = 100;
options.quit_multiple = 10;
options.irf_horizon = 24;
options.reject_nonstationary = false;
create_restrictions;
options.dummy_obs = false;
% presample dummy observations
presampledummyobs = 100;

% PLOT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
keep_figure = 0;     % if TRUE, newly computed IRFs are overlaid on Figure 2
second_color = 'm';     % color to plot in if you 'keep figure'
plot_names = options.names;
irf_horizon = 24;
quantiles_to_plot = [0.05 0.5 0.95];
plot_integrated = {};
filled_irf_bands = true;
plot_posterior_mode = true;

% variable to shock
shock_to = {'DEP3M', 'LENDRATE', 'LPRIVCRED'};
shock_permanent = 0;
shock_value = {0.01, 0.01, 0.01}; % normalizes period 1 response to shock_to variable. If 
                   % zero, then no normalization.

%-------------------------------NOT USED-----------------------------------
% 
options.forecasting = 0;     % 1: Compute h-step ahead predictions, 0: no prediction
options.forecast_method = 1; % 0: Direct forecasts 
                                  % 1: Iterated forecasts
options.forecast_horizon = 1; % Number of periods to forecast ahead
options.h_oos = 1;  % int > 1 do h_oos recursive out of sample forecasts; 
                     % the default is to do 1 OOS prediction.
options.h_score = 12;        % number of initial OOS periods to score (<= h_oos)
options.ksbwidth = 0.4;      % ksdensity bandwidth parameter (if plotting single series)
%---------------------------END USER SETTINGS------------------------------

disp( ' ' );
disp( 'Running Bayesian VAR' );
disp( '====================' );
disp( datestr( now ) );


%% -----------LOAD DATA & CLEANUP OPTIONS--------------------------------
[Yraw, Ynames, cal] = GetVARData( data_dir, data_file, data_tab );
[Yraw, Zraw] = ProcessVAROptions( Yraw, options, Ynames, shock_to );
[Yraw, Zraw, cal] = GetEstimationSubsample(Yraw, Zraw, cal, EST_START_YEAR, EST_START_Q, EST_END_YEAR, EST_END_Q);

% Outputs Y1 and X1 are the raw data, manipulated so that the full sample
% is available, after taking away observations to create lags and to
% account for possibly h-step ahead specifications such as y(t+h)=B y(t) +
[ Y, X, Tstar, options ] = DoRawDataManipulation( ...
                                    Yraw, Zraw, options, cal );

%% RUN BVAR--------------------------------------------------------------
initBVAR( 20121126 );
% Set the prior
[priorparams, options] = SetSimsZhaPrior( Y, X, prior_name, options );
% Obtain posterior
postest = GibbsBVAR( Y, X, priorparams, options );

% Find HPD parameter values
[~, max_index] = max( postest.loglik );
mA_hpd = squeeze( postest.A_post(max_index,:,:) );
mF_hpd = squeeze( postest.F_post(max_index,:,:) );
mB_hpd = squeeze( postest.B_post(max_index,:,:) );

% Obtain marginal likelihood
log10mdd = GetSimsZhaMDD( Y, X, priorparams, postest, options, mA_hpd, mF_hpd );
postest.chibmdd = log10mdd;


%% DISPLAY RESULTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% HPD parameters
disp( 'mA = ' )
disp( mA_hpd );
disp( 'mF = ' )
disp( mF_hpd );
disp( 'mB = ' )
disp( mB_hpd );

disp( 'System Eigenvalues: ' );
disp( sort( abs( eig( GetCompanionForm( ...
    mB_hpd, options ) ) ) ) );

fprintf( 'Log_10 Chib marginal likelihood = %.3f\n\n', log10mdd );


%% Plot IRFs
M = numel( options.names );

for plot_num = 1:numel( shock_to )
    shock_index = find( strcmp( shock_to{plot_num}, options.names ) );
    font_size = 12; % default (presentation mode default = 20)

    % avoid plotting nonsense...
    if keep_figure
        filled_irf_bands = false;
    end

    % Check whether some of the variables in the impulse responses are supposed 
    % to be 'zeroed out'
    index_zeroed_variables = [];
    plot_zero_restriction_case = 0; % default

    % Check which variables are to be plotted
    vars_to_plot = zeros( 0, numel( plot_names ) );

    if 0 < numel( plot_names )
        for i = 1:numel(plot_names)
            [b_flag, i_loc ] = ismember( plot_names(i), options.names );
            if b_flag
                vars_to_plot(i) = i_loc;
            else
                warning( ['Cannot plot response of ' plot_names{i} ': not in model' ] );
                plot_names = options.names;
                vars_to_plot = 1:M;
            end
        end
    else
        plot_names = options.names;
        vars_to_plot = 1:M;
    end

    num_vars_to_plot = length( vars_to_plot );


    % Compute quantiles of IRF distribution (all variables)
    imp = postest.IRF;
    imp_qu16 = zeros( 0, 0, irf_horizon );
    imp_qu50 = zeros( 0, 0, irf_horizon );
    imp_qu84 = zeros( 0, 0, irf_horizon );

    % Compute the impulse-response
    imp_mode = impulse( reshape( mB_hpd(1:M*options.nlags,:)', M, M, options.nlags ), ...
            inv(mA_hpd), irf_horizon );

    if ~(0==shock_value{plot_num})
        norm_constant = ...
                shock_value{plot_num}/squeeze( imp_mode(shock_index,shock_index,1) );
        imp_mode = imp_mode*norm_constant;
    end

    % Check if we need to plot a normalized shock (recall, we're plotting only
    % the response to a single shock, given by shock_index)
    if ~( 0 == shock_value{plot_num} )
        % for each draw, get the scalar size of the shock for that draw,
        % and create a normalizing scalar quantity that will set the shock
        % to shock_value
        norm_constant = shock_value{plot_num} ./ imp(:,shock_index,shock_index,1);

        % Iterate over horizon, multiplying the shock_index column in each
        % draw by norm_scalar. The matrix imp(:,:,shock_index,h) has the
        % shock_index column of the impulse matrix in transpose along each row, 
        % with each draw being a new row.
        for h = 1:irf_horizon
            imp(:,:,shock_index,h) = diag( norm_constant )*imp(:,:,shock_index,h); 
        end

    end

    % Check if some of the responses need to be integrated
    for j = 1:numel( options.names )
        if ( ismember( options.names{j}, plot_integrated ) && ismember( options.names{j}, plot_names ) )
            for h=2:irf_horizon
                imp(:,j,shock_index,h) = imp(:,j,shock_index,h) + imp(:,j,shock_index,h-1);
                if plot_posterior_mode
                    imp_mode(j,shock_index,h) = imp_mode(j,shock_index,h) + imp_mode(j,shock_index,h-1);
                end
            end
        end
    end

    % Get the specified quantiles of the posterior IRF distribution (which may
    % have been integrated)
    for j = 1:M
        for h = 1:irf_horizon
            qu_out = ...
               quantile( imp( :, j, shock_index, h ), quantiles_to_plot );
            imp_qu16(j,shock_index,h) = 100*qu_out( 1 );
            imp_qu50(j,shock_index,h) = 100*qu_out( 2 );
            imp_qu84(j,shock_index,h) = 100*qu_out( 3 );
        end
    end

    % Pointwise percentile bands ----------------------------------------------
    f = figure(plot_num);
    set(f,'Name', cat(2, 'Shock to ', shock_to{plot_num}));
    

    if (~keep_figure)
        clf(plot_num);
        line_color = 'b';
    else
        line_color = second_color;
    end

    % init
    counter = 1;

    % size of subplot; checks the number of variables to plot and adjusts
    % the size of the subplot accordingly...
    plot_rows = ceil( sqrt(num_vars_to_plot) );
    plot_cols = floor( sqrt(num_vars_to_plot) );
    if ( num_vars_to_plot > plot_rows*plot_cols )
        plot_cols = plot_rows;
    end

    if (0 < irf_horizon)    
        % do irf plot
        for j = vars_to_plot

            subplot( plot_rows, plot_cols, counter ); 
            hold on;

            if ( filled_irf_bands )
                lower_data = rescale{j}*squeeze( imp_qu16(j,shock_index,:) );
                upper_data = rescale{j}*squeeze( imp_qu84(j,shock_index,:)-imp_qu16(j,shock_index,:) );
                hb = area(1:irf_horizon,[lower_data, upper_data], ...
                                        min(lower_data),'LineStyle','none');
                plot( 1:irf_horizon, rescale{j}*squeeze( imp_qu50(j,shock_index,:) ), ...
                                        'k', 'LineWidth', 1.0 );
                if plot_posterior_mode
                    plot( 1:irf_horizon, rescale{j}*squeeze( 100*imp_mode(j,shock_index,:) ), ...
                                        'm', 'LineWidth', 1.0 );
                end

                set(hb(1),'FaceColor',[1,1,1]);
                set(hb(2),'FaceColor',[0.85,0.85,0.85]);
                set(gca,'layer','top');

                plot( 1:irf_horizon, zeros(1,irf_horizon), '-g' );                                                                       
            else
                plot( 1:irf_horizon, zeros(1,irf_horizon), '-g' );
                plot( 1:irf_horizon, rescale{j}*squeeze(imp_qu16( j, shock_index, : )), ...
                                        [':' line_color], 'LineWidth', 1.0 );
                plot( 1:irf_horizon, rescale{j}*squeeze(imp_qu50( j, shock_index, : )), ...
                                        ['-' line_color], 'LineWidth', 2.0 );
                plot( 1:irf_horizon, rescale{j}*squeeze(imp_qu84( j, shock_index, : )), ...
                                        [':' line_color], 'LineWidth', 1.0 );            
            end

            hold off;
            axis tight;        
            title( options.long_names(j),'FontSize',font_size );
            ylabel( units(j),'Interpreter','none','FontSize',14);
            if counter > (plot_rows-1)*plot_cols
                xlabel('Quarters since shock','FontSize',14);
            end
            set(gca,'FontSize',14);

            counter = counter+1;
        end

    end
end

