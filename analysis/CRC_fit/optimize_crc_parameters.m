function [summary_table, parameter_table] = optimize_crc_parameters(case_indices, num_starts, max_iterations)
%OPTIMIZE_CRC_PARAMETERS Refit selected CRC cases without replacing the baseline.
% No arguments fits all eight cases. Example: optimize_crc_parameters(1,3,120).
if nargin < 1; case_indices = 1:8; end
if nargin < 2; num_starts = 3; end
if nargin < 3; max_iterations = 120; end
module_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(fileparts(module_dir),'src')));
cases = crc_cases();
cases = cases(case_indices);
settings = struct('num_starts',num_starts,'max_iterations',max_iterations, ...
    'objective','ordinary-volume-SSE');
case_ids = string({cases.id});
optimized_parameters = cell(1,numel(cases));
summary_rows = cell(numel(cases),6);
parameter_rows = cell(0,6);
for i = 1:numel(cases)
    fit = fit_crc_case(cases(i), NumStarts=num_starts, ...
        MaxIterations=max_iterations);
    optimized_parameters{i} = fit.optimized_params;
    summary_rows(i,:) = {fit.case_id,fit.optimized.metrics.rmse, ...
        fit.optimized.metrics.r2,fit.solver.resnorm,fit.solver.exitflag,fit.solver.method};
    for j = 1:numel(fit.parameter_names)
        parameter_rows(end+1,:) = {fit.case_id,fit.parameter_names{j}, ...
            fit.optimizer_initial_values(j),fit.solver.x(j), ...
            fit.lower_bound(j),fit.upper_bound(j)}; %#ok<AGROW>
    end
    fprintf('CRC optimization finished: %s\n',fit.case_id);
end
summary_table = cell2table(summary_rows,'VariableNames', ...
    {'case_id','rmse_cm3','r2','objective_value','exitflag','solver'});
parameter_table = cell2table(parameter_rows,'VariableNames', ...
    {'case_id','parameter','initial_value','fitted_value','lower_bound','upper_bound'});
output_dir = fullfile(module_dir,'results');
if ~exist(output_dir,'dir'); mkdir(output_dir); end
save(fullfile(output_dir,'refit_results.mat'),'case_ids','optimized_parameters', ...
    'settings','summary_table','parameter_table');
writetable(summary_table,fullfile(output_dir,'refit_summary.csv'));
writetable(parameter_table,fullfile(output_dir,'refit_parameters.csv'));
end
