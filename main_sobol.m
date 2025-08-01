%   This script uses the UQLab toolbox to perform a variance-based global
%   sensitivity analysis (GSA) on the nonlocal Fisher-KPP model. 

%   DEPENDENCIES:
%   - MATLAB
%   - UQLab (www.uqlab.com)
%   - The solver files: 'uq_ff_wrapper.m' (for the Front-Fixing method) and
%   'uq_ft_rk_wrapper.m' (for the FTRK method)

clearvars; close all; clc;
rng(101, 'twister');
   uqlab;


% 1. Define Fixed Parameters
FixedParams.M = 100;  
FixedParams.T = 5;   

% 2. Computational Model (using the wrapper)
%ModelOpts.mFile = 'uq_ft_rk_wrapper'; 
ModelOpts.mFile = 'uq_ff_wrapper'; 
ModelOpts.Parameters = FixedParams;    
ModelOpts.isVectorized = false;       

myModel = uq_createModel(ModelOpts);


% 3. Probabilistic Input Model
N_distr = 3; % 1 - Uniform, 2 - Gaussian, 3 - Lognormal


h0_distr =      {'Uniform',[0.5 3];     'Gaussian', [1.75 0.41667];     'Lognormal', [0 0.5]};
mu_distr =      {'Uniform',[0.1 2];     'Gaussian', [1.05 0.31667];     'Lognormal', [-0.5 0.5]};
alpha2_distr =  {'Uniform',[1 4];       'Gaussian', [2.5 0.5];          'Lognormal', [0.8 0.25]};
r_distr =       {'Uniform',[0.5 2];     'Gaussian', [1.25 0.25];        'Lognormal', [0 0.25]};
k_distr =       {'Uniform',[0.25 4];    'Gaussian', [2.125 0.625];      'Lognormal', [0.5 0.35]};
sigma_distr =   {'Uniform',[0.1 2];     'Gaussian', [1.05 0.31667];     'Lognormal', [-0.5 0.5]};


InputOpts.Marginals(1).Name = 'h0';    
InputOpts.Marginals(1).Type = h0_distr{N_distr,1};
InputOpts.Marginals(1).Parameters = h0_distr{N_distr,2};


InputOpts.Marginals(2).Name = 'mu';    
InputOpts.Marginals(2).Type = mu_distr{N_distr,1};
InputOpts.Marginals(2).Parameters = mu_distr{N_distr,2};


InputOpts.Marginals(3).Name = 'alpha2'; 
InputOpts.Marginals(3).Type = alpha2_distr{N_distr,1};
InputOpts.Marginals(3).Parameters = alpha2_distr{N_distr,2};


InputOpts.Marginals(4).Name = 'r';      
InputOpts.Marginals(4).Type = r_distr{N_distr,1};
InputOpts.Marginals(4).Parameters = r_distr{N_distr,2};


InputOpts.Marginals(5).Name = 'k';      
InputOpts.Marginals(5).Type = k_distr{N_distr,1};
InputOpts.Marginals(5).Parameters = k_distr{N_distr,2};


InputOpts.Marginals(6).Name = 'sigma';   
InputOpts.Marginals(6).Type = sigma_distr{N_distr,1};
InputOpts.Marginals(6).Parameters = sigma_distr{N_distr,2};

myInput = uq_createInput(InputOpts);

fprintf('Input created\n');

output_dir = 'SA_Figures';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end



fprintf('MC strated...\n')
tic
SobolOpts.Type = 'Sensitivity';
SobolOpts.Method = 'Sobol';
SobolOpts.Sobol.Order = 1;
SobolOpts.Sobol.SampleSize = 1000;
mySobolAnalysisMC = uq_createAnalysis(SobolOpts);
mySobolResultsMC = mySobolAnalysisMC.Results;
fprintf('\t done!\n')
toc


fprintf('PCE strated...\n')
tic
PCEOpts.Type = 'Metamodel';
PCEOpts.MetaType = 'PCE';
PCEOpts.FullModel = myModel;
PCEOpts.Degree = 5;
PCEOpts.ExpDesign.NSamples = 200;
myPCE = uq_createModel(PCEOpts);

mySobolAnalysisPCE = uq_createAnalysis(SobolOpts);
mySobolResultsPCE = mySobolAnalysisPCE.Results;
fprintf('\t done!\n')
toc

fprintf('LRA strated...\n')
tic
LRAOpts.Type = 'Metamodel';
LRAOpts.MetaType = 'LRA';
LRAOpts.FullModel = myModel;
LRAOpts.Rank = 1:20;
LRAOpts.Degree = 1:20;
LRAOpts.ExpDesign.NSamples = 200;
myLRA = uq_createModel(LRAOpts);
mySobolAnalysisLRA = uq_createAnalysis(SobolOpts);
mySobolResultsLRA = mySobolAnalysisLRA.Results;
fprintf('\t done!\n')
toc
uq_print(mySobolAnalysisMC)
uq_print(mySobolAnalysisPCE)
uq_print(mySobolAnalysisLRA)

SobolTotal = [...
    mySobolResultsMC.Total ...
    mySobolResultsPCE.Total ...
    mySobolResultsLRA.Total];
SobolFirstOrder = [...
    mySobolResultsMC.FirstOrder ...
    mySobolResultsPCE.FirstOrder ...
    mySobolResultsLRA.FirstOrder];

% --- Plotting Section ---
all_distribution_types = {'Uniform', 'Gaussian', 'Lognormal'};
current_distribution_name = all_distribution_types{N_distr};
num_vars = length(InputOpts.Marginals);
fig_total = uq_figure('Name', sprintf('Total Sobol Indices (%s Distributions)', current_distribution_name));
uq_bar(1:num_vars, SobolTotal, 0.8); 
ylim([0 1]);
xlim([0 num_vars + 1]);
xlabel('Variable name');
ylabel('Total Sobol'' indices ($S_{Ti}$)');
title(sprintf('Total Sobol'' Indices - %s Inputs', current_distribution_name));
set(gca, 'XTick', 1:num_vars, 'XTickLabel', mySobolResultsPCE.VariableNames);
uq_legend({sprintf('MC (N=%d)', mySobolResultsMC.Cost), ...
           sprintf('PCE (N_{ED}=%d)', myPCE.ExpDesign.NSamples), ...
           sprintf('LRA (N_{ED}=%d)', myLRA.ExpDesign.NSamples)}, ...
           'Location', 'northeast');
grid on;

filename_total_eps = fullfile(output_dir, sprintf('Sobol_Total_Indices_%s.eps', current_distribution_name));
print(fig_total, filename_total_eps, '-depsc', '-r300'); % Save as EPS, 300 dpi
fprintf('Saved Total Sobol figure to: %s\n', filename_total_eps);

fig_first = uq_figure('Name', sprintf('First-order Sobol Indices (%s Distributions)', current_distribution_name));
uq_bar(1:num_vars, SobolFirstOrder, 0.8); % barWidth directly in uq_bar
ylim([0 1]);
xlim([0 num_vars + 1]);
xlabel('Variable name');
ylabel('First-order Sobol'' indices ($S_i$)');
title(sprintf('First-order Sobol'' Indices - %s Inputs', current_distribution_name));
set(gca, 'XTick', 1:num_vars, 'XTickLabel', mySobolResultsPCE.VariableNames);
uq_legend({sprintf('MC (N=%d)', mySobolResultsMC.Cost), ...
           sprintf('PCE (N_{ED}=%d)', myPCE.ExpDesign.NSamples), ...
           sprintf('LRA (N_{ED}=%d)', myLRA.ExpDesign.NSamples)}, ...
           'Location', 'northeast');
grid on;

filename_first_eps = fullfile(output_dir, sprintf('Sobol_First_Order_Indices_%s.eps', current_distribution_name));
print(fig_first, filename_first_eps, '-depsc', '-r300'); % Save as EPS, 300 dpi
fprintf('Saved First-order Sobol figure to: %s\n', filename_first_eps);



% --- Combined First-Order and Total Sobol Indices Plot ---
num_vars = length(InputOpts.Marginals);
parameter_names = mySobolResultsPCE.VariableNames; % Get from one of the results

fig_combined = uq_figure('Name', sprintf('Sobol Indices (%s Distributions)', current_distribution_name));
hold on;


plot_data = zeros(num_vars, 6);
for i = 1:num_vars
    plot_data(i, 1) = SobolFirstOrder(i, 1); % MC S_i
    plot_data(i, 2) = SobolTotal(i, 1);      % MC S_Ti
    plot_data(i, 3) = SobolFirstOrder(i, 2); % PCE S_i
    plot_data(i, 4) = SobolTotal(i, 2);      % PCE S_Ti
    plot_data(i, 5) = SobolFirstOrder(i, 3); % LRA S_i
    plot_data(i, 6) = SobolTotal(i, 3);      % LRA S_Ti
end

b = bar(1:num_vars, plot_data, 'grouped'); 
hold off;

color_Si_MC   = [0 0.4470 0.7410]; % Blue
color_STi_MC  = [0.3010 0.7450 0.9330]; % Light Blue
color_Si_PCE  = [0.8500 0.3250 0.0980]; % Orange
color_STi_PCE = [0.9290 0.6940 0.1250]; % Light Orange/Yellow
color_Si_LRA  = [0.4660 0.6740 0.1880]; % Green
color_STi_LRA = [0.7000 0.8500 0.4000]; % Light Green

set(b(1), 'FaceColor', color_Si_MC, 'DisplayName', sprintf('S_i MC (N=%d)', mySobolResultsMC.Cost));
set(b(2), 'FaceColor', color_STi_MC, 'DisplayName', sprintf('S_{Ti} MC'));
set(b(3), 'FaceColor', color_Si_PCE, 'DisplayName', sprintf('S_i PCE (N_{ED}=%d)', myPCE.ExpDesign.NSamples));
set(b(4), 'FaceColor', color_STi_PCE, 'DisplayName', sprintf('S_{Ti} PCE'));
set(b(5), 'FaceColor', color_Si_LRA, 'DisplayName', sprintf('S_i LRA (N_{ED}=%d)', myLRA.ExpDesign.NSamples));
set(b(6), 'FaceColor', color_STi_LRA, 'DisplayName', sprintf('S_{Ti} LRA'));

% Set axes limits
ylim([0 1.05]); 
xlim([0.5 num_vars + 0.5]);

xlabel('Input Variable');
ylabel('Sobol'' Indices');
%title(sprintf('First-Order and Total Sobol'' Indices - %s Inputs', current_distribution_name));
set(gca, 'XTick', 1:num_vars, 'XTickLabel', parameter_names);

legend('show', 'Location', 'northeast', 'FontSize', 8); % Adjust FontSize if needed
grid on;
box on;

filename_combined_eps = fullfile(output_dir, sprintf('Sobol_Combined_Indices_%s.eps', current_distribution_name));
print(fig_combined, filename_combined_eps, '-depsc', '-r300'); % Save as EPS, 300 dpi
fprintf('Saved Combined Sobol figure to: %s\n', filename_combined_eps);