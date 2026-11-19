function abundance_table = calculate_abundance_robustness(input_file)
% Recompute msImmunity AUC at the five abundance constants in Fig. S8A.
abundance_values = [0.3,0.4,0.5,0.6,0.7];
baseline_index = find(abundance_values==0.5);
module_dir = fileparts(mfilename('fullpath'));
analysis_dir = fileparts(module_dir);
addpath(genpath(fullfile(analysis_dir,'src')));
vct_dir = fullfile(analysis_dir,'VCT_analysis');
input_dir = fullfile(vct_dir,'input_data');
if nargin < 1
    input_file = fullfile(vct_dir,'results','vct_patient_results.csv');
end
design = load(fullfile(input_dir,'VP_parameter_matrix.mat'),'samples');
patients = readtable(input_file,'TextType','string');
n = size(design.samples,1);
assert(height(patients)==4*n,'Expected four treatment records per virtual patient.');
regimen_names = ["Xelox","Xelox + RT","Xelox + anti-PD-1","Xelox + RT + anti-PD-1"];
therapy_files = {'Therapy_xelox.mat','Therapy_xelox_IR.mat', ...
    'Therapy_xelox_antiPD1.mat','Therapy_xelox_antiPD1_IR.mat'};
n_values = numel(abundance_values);
auc_ms = zeros(n_values,numel(regimen_names));
for k = 1:numel(regimen_names)
    current = sortrows(patients(patients.regimen==regimen_names(k),:),'patient_id');
    assert(isequal(current.patient_id,(1:n)'),'Patient IDs must match the saved design.');
    therapy = load_single_variable(fullfile(input_dir,therapy_files{k}));
    scores = zeros(n,n_values);
    for i = 1:n
        params = apply_virtual_patient_parameters(design.samples(i,:));
        all_scores = Quantify(params,therapy,abundance_values);
        scores(i,:) = all_scores(:,1)';
    end
    reference = current.msImmunity;
    assert(all(abs(scores(:,baseline_index)-reference)./max(abs(reference),realmin)<1e-10), ...
        'The 0.5 calculation must reproduce the saved msImmunity scores.');
    for j = 1:n_values
        [~,~,~,auc_ms(j,k)] = roc_curve(logical(current.responder),scores(:,j));
    end
    fprintf('Cell abundance: %s, %d patients evaluated.\n',regimen_names(k),n);
end
abundance_table = table(repelem(regimen_names(:),n_values), ...
    repmat(abundance_values(:),numel(regimen_names),1),auc_ms(:), ...
    'VariableNames',{'regimen','abundance_constant','auc_ms'});
output_dir = fullfile(module_dir,'results');
if ~exist(output_dir,'dir'); mkdir(output_dir); end
writetable(abundance_table,fullfile(output_dir,'abundance_robustness_summary.csv'));
end
