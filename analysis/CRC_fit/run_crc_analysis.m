function [results, summary_table] = run_crc_analysis()
%RUN_CRC_ANALYSIS Recalculate all eight CRC cases with their final parameters.
module_dir = fileparts(mfilename('fullpath'));
analysis_dir = fileparts(module_dir);
addpath(genpath(fullfile(analysis_dir,'src')));
cases = crc_cases();
fixed = load(fullfile(module_dir,'input_data','final_parameters.mat'));
output_dir = fullfile(module_dir,'results');
if ~exist(output_dir,'dir'); mkdir(output_dir); end

results = struct([]);
summary_rows = cell(numel(cases),3);
parameter_rows = cell(0,3);
for i = 1:numel(cases)
    assert(string(cases(i).id)==fixed.case_ids(i),'CRC case order differs from the saved parameters.');
    params = fixed.final_parameters{i};
    evaluation = evaluate_crc_case(cases(i),params);
    results(i).case_id = cases(i).id;
    results(i).parameters = params;
    results(i).evaluation = evaluation;
    summary_rows(i,:) = {cases(i).id,evaluation.metrics.rmse,evaluation.metrics.r2};
    for j = 1:numel(fixed.parameter_names{i})
        name = fixed.parameter_names{i}{j};
        parameter_rows(end+1,:) = {cases(i).id,name,get_parameter(params,name)}; %#ok<AGROW>
    end
end
summary_table = cell2table(summary_rows,'VariableNames',{'case_id','rmse_cm3','r2'});
parameter_table = cell2table(parameter_rows,'VariableNames',{'case_id','parameter','value'});
save(fullfile(output_dir,'crc_results.mat'),'results','summary_table','parameter_table');
writetable(summary_table,fullfile(output_dir,'crc_fit_summary.csv'));
writetable(parameter_table,fullfile(output_dir,'crc_final_parameters.csv'));
fprintf('CRC: reproduced %d cases using fixed final parameters.\n',numel(cases));
end
