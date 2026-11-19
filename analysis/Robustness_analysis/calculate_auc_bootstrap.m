function auc_table = calculate_auc_bootstrap(input_file)
%CALCULATE_AUC_BOOTSTRAP Paired patient bootstrap; numerical output only.
% A shared patient draw preserves all four regimens and all three scores.
bootstrap_replicates = 5000;
bootstrap_seed = 20261107;
module_dir = fileparts(mfilename('fullpath'));
analysis_dir = fileparts(module_dir);
addpath(genpath(fullfile(analysis_dir,'src')));
if nargin < 1
    input_file = fullfile(analysis_dir,'VCT_analysis','results','vct_patient_results.csv');
end
output_dir = fullfile(module_dir,'results');
if ~exist(output_dir,'dir'); mkdir(output_dir); end
regimen_names = ["Xelox","Xelox + RT","Xelox + anti-PD-1","Xelox + RT + anti-PD-1"];
score_names = ["msImmunity","tsImmunity","molsImmunity"];
patients = readtable(input_file,'TextType','string');
patient_ids = sort(unique(patients.patient_id));
n = numel(patient_ids);
assert(height(patients)==4*n,'Expected exactly four records per patient.');
scores = zeros(n,3,4);
labels = false(n,4);
auc_estimate = zeros(3,4);
for k=1:4
    current = sortrows(patients(patients.regimen==regimen_names(k),:),'patient_id');
    assert(isequal(current.patient_id,patient_ids),'Missing or duplicate patient/regimen record.');
    assert(all(ismember(current.responder,[0,1])),'Response labels must be binary.');
    assert(all(current.responder==(current.diameter_change_percent<=-30)), ...
        'Saved response labels differ from the fixed equivalent-diameter definition.');
    scores(:,:,k) = current{:,cellstr(score_names)};
    labels(:,k) = logical(current.responder);
    assert(all(isfinite(scores(:,:,k)),'all'),'All three scores must be finite.');
    assert(any(labels(:,k)) && any(~labels(:,k)),'Both outcome classes are required.');
    for m=1:3
        [~,~,~,auc_estimate(m,k)] = roc_curve(labels(:,k),scores(:,m,k));
    end
end
% Unstratified patient-level resampling preserves the joint response pattern
% across all four regimens. Separate responder-stratified draws per regimen
% would break that pairing.
rng(bootstrap_seed,'twister');
bootstrap_indices = uint32(randi(n,n,bootstrap_replicates));
bootstrap_auc = nan(bootstrap_replicates,3,4);
for b=1:bootstrap_replicates
    indices = double(bootstrap_indices(:,b));
    for k=1:4
        sampled_labels = labels(indices,k);
        if ~any(sampled_labels) || all(sampled_labels)
            continue % Undefined AUC; retain NaN and report the count below.
        end
        for m=1:3
            [~,~,~,bootstrap_auc(b,m,k)] = roc_curve(sampled_labels,scores(indices,m,k));
        end
    end
end
auc_rows = cell(12,9);
for k=1:4
    valid = all(isfinite(bootstrap_auc(:,:,k)),2);
    assert(any(valid),'No valid bootstrap replicates for a regimen.');
    auc_ci = prctile(bootstrap_auc(valid,:,k),[2.5,97.5],1);
    for m=1:3
        auc_rows((k-1)*3+m,:) = {regimen_names(k),score_names(m),auc_estimate(m,k), ...
            auc_ci(1,m),auc_ci(2,m),n,sum(labels(:,k)),sum(valid),sum(~valid)};
    end
end
auc_table = cell2table(auc_rows,'VariableNames',{'regimen','score','auc','ci_low','ci_high', ...
    'patient_count','responder_count','valid_bootstrap_count','undefined_bootstrap_count'});
writetable(auc_table,fullfile(output_dir,'auc_bootstrap_summary.csv'));
fprintf('AUC bootstrap: %d paired patients, %d replicates.\n',n,bootstrap_replicates);
end
