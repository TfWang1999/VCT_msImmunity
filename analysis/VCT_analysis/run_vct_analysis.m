function [report, patient_table, auc_table] = run_vct_analysis(cohort_root)
%RUN_VCT_ANALYSIS Analyze VCT cohorts using RECIST-inspired diameter thresholds.
% Baseline-relative endpoint classes: PR <= -30%, SD between -30% and +20%,
% PD >= +20%. PR denotes response; CR is not assigned separately.

module_dir = fileparts(mfilename('fullpath'));
analysis_dir = fileparts(module_dir);
addpath(genpath(fullfile(analysis_dir,'src')));
if nargin < 1 || strlength(string(cohort_root)) == 0
    cohort_root = fullfile(module_dir,'results','cohort');
end
cohort_root = char(cohort_root);
results_dir = fullfile(module_dir,'results');
if ~exist(results_dir,'dir'); mkdir(results_dir); end
design = load(fullfile(module_dir,'input_data','VP_parameter_matrix.mat'),'samples');
expected_ids = (1:size(design.samples,1))';

folder_names = {'VP_Xelox','VP_Xelox_IR','VP_Xelox_antipd1','VP_Xelox_IR_antipd1'};
therapy_files = {'Therapy_xelox.mat','Therapy_xelox_IR.mat', ...
    'Therapy_xelox_antiPD1.mat','Therapy_xelox_antiPD1_IR.mat'};
regimen_names = {'Xelox','Xelox + RT','Xelox + anti-PD-1','Xelox + RT + anti-PD-1'};
score_names = {'msImmunity','tsImmunity','molsImmunity'};

report.regimen_names = regimen_names;
report.score_names = score_names;
report.response_definition = 'equivalent-diameter change <= -30%';
report.progression_definition = 'equivalent-diameter change >= +20%';
patient_rows = cell(0,11);
auc_rows = cell(0,6);

for regimen_index = 1:4
    regimen_dir = fullfile(cohort_root,folder_names{regimen_index});
    files = dir(fullfile(regimen_dir,'VP*.mat'));
    patient_numbers = arrayfun(@(item) sscanf(item.name,'VP%d.mat'),files(:));
    [patient_numbers,order] = sort(patient_numbers);
    if ~isequal(patient_numbers,expected_ids)
        error('run_vct_analysis:IncompleteCohort', ...
            'Incomplete virtual cohort for %s. Run run_virtual_cohort first.',folder_names{regimen_index});
    end
    files = files(order);
    therapy = load_single_variable(fullfile(module_dir,'input_data',therapy_files{regimen_index}));
    n = numel(files);
    scores = zeros(n,3);
    volume_change = zeros(n,1);
    diameter_change = zeros(n,1);
    denominators = zeros(n,2);
    for patient_index = 1:n
        archive_item = load_single_variable(fullfile(files(patient_index).folder,files(patient_index).name));
        initial_volume = archive_item.Vtumor_list(1);
        final_volume = archive_item.Vtumor_list(end);
        volume_change(patient_index) = final_volume ./ initial_volume - 1;
        diameter_change(patient_index) = (final_volume ./ initial_volume).^(1/3) - 1;
        [scores(patient_index,:),~,diagnostics] = Quantify(archive_item.params,therapy);
        denominators(patient_index,:)=[diagnostics.ms_denominator,diagnostics.ts_denominator];
        assert(all(isfinite(scores(patient_index,:))) && ...
            all(isfinite(archive_item.Vtumor_list)) && all(archive_item.Vtumor_list>0), ...
            'run_vct_analysis:InvalidResult','Scores and tumor volumes must be finite; volumes must be positive.');
    end
    responder = diameter_change <= -0.30;
    progressive = diameter_change >= 0.20;
    response_class = repmat({'SD'},n,1);
    response_class(responder) = {'PR'};
    response_class(progressive) = {'PD'};
    report.regimens(regimen_index).name = regimen_names{regimen_index};
    report.regimens(regimen_index).patient_id = patient_numbers;
    report.regimens(regimen_index).scores = scores;
    report.regimens(regimen_index).volume_change = volume_change;
    report.regimens(regimen_index).diameter_change = diameter_change;
    report.regimens(regimen_index).responder = responder;
    report.regimens(regimen_index).response_class = response_class;
    report.regimens(regimen_index).response_rate = mean(responder);
    report.regimens(regimen_index).score_denominators = denominators;
    report.regimens(regimen_index).spearman_ms = corr(scores(:,1),diameter_change,'Type','Spearman');
    for score_index = 1:3
        [~,~,~,auc] = roc_curve(responder,scores(:,score_index));
        report.auc(score_index,regimen_index) = auc;
        auc_rows(end+1,:) = {score_names{score_index},regimen_names{regimen_index}, ... %#ok<AGROW>
            auc,n,sum(responder),mean(responder)};
    end
    for patient_index = 1:n
        patient_rows(end+1,:) = {regimen_names{regimen_index},patient_numbers(patient_index), ... %#ok<AGROW>
            scores(patient_index,1),scores(patient_index,2),scores(patient_index,3), ...
            100*volume_change(patient_index),100*diameter_change(patient_index), ...
            responder(patient_index),response_class{patient_index}, ...
            denominators(patient_index,1),denominators(patient_index,2)};
    end
    fprintf('%s: response rate %.1f%%; AUC %.4f %.4f %.4f\n', ...
        regimen_names{regimen_index},100*mean(responder),report.auc(:,regimen_index));
    fprintf('  Spearman rho (msImmunity vs diameter change): %.6f\n',report.regimens(regimen_index).spearman_ms);
end

patient_table = cell2table(patient_rows,'VariableNames', ...
    {'regimen','patient_id','msImmunity','tsImmunity','molsImmunity', ...
    'volume_change_percent','diameter_change_percent','responder','response_class', ...
    'ms_denominator','ts_denominator'});
auc_table = cell2table(auc_rows,'VariableNames', ...
    {'score','regimen','auc','cohort_size','responder_count','response_rate'});
writetable(patient_table,fullfile(results_dir,'vct_patient_results.csv'));
writetable(auc_table,fullfile(results_dir,'vct_auc.csv'));
save(fullfile(results_dir,'vct_analysis_results.mat'),'report','patient_table','auc_table');
end
