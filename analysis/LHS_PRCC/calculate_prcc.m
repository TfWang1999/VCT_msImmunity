function result_table = calculate_prcc(rerun_simulations, workers)
%CALCULATE_PRCC Numerical PRCC for the saved 2,000 x 76 P01 LHS design.
% Default: recompute PRCC and intervals from the saved terminal volumes.
% Pass true to rerun the same samples in all five scenarios (no new draws).
if nargin < 1; rerun_simulations = false; end
if nargin < 2; workers = 0; end
bootstrap_replicates = 500;

module_dir = fileparts(mfilename('fullpath'));
analysis_dir = fileparts(module_dir);
addpath(genpath(fullfile(analysis_dir,'src')));
design = load(fullfile(module_dir,'input_data','p01_prcc_design.mat'));
output_dir = fullfile(module_dir,'results');
if ~exist(output_dir,'dir'); mkdir(output_dir); end
X = design.parameter_samples;
[n,d] = size(X);
scenario_names = design.scenario_names;
if rerun_simulations
    terminal_volume = zeros(n,numel(scenario_names));
    for k = 1:numel(scenario_names)
        experiment = design.experiments{k};
        values = zeros(n,1);
        parfor (i = 1:n, workers)
            params = design.reference_params;
            for j = 1:d
                params = set_parameter(params,design.parameter_names(j),X(i,j));
            end
            params.TSM.b = params.TSM.a/10;
            [~,volume] = run_simulation(experiment,params);
            values(i) = volume(end);
        end
        terminal_volume(:,k) = values;
        fprintf('PRCC: %s, %d saved samples reproduced.\n',scenario_names(k),n);
    end
else
    saved = load(fullfile(module_dir,'input_data','p01_terminal_volume.mat'),'terminal_volume');
    terminal_volume = saved.terminal_volume;
end
assert(isequal(size(terminal_volume),[n,numel(scenario_names)]));
prcc_values = compute_prcc(X,terminal_volume);

% Each resampled row keeps all five outcomes together. Re-rank after drawing.
rng(design.random_seed+1,'twister');
indices = randi(n,n,bootstrap_replicates);
bootstrap_prcc = zeros(d,numel(scenario_names),bootstrap_replicates);
parfor (b = 1:bootstrap_replicates, workers)
    ids = indices(:,b);
    bootstrap_prcc(:,:,b) = compute_prcc(X(ids,:),terminal_volume(ids,:));
end
interval = prctile(bootstrap_prcc,[2.5,97.5],3);
ci_lower = interval(:,:,1);
ci_upper = interval(:,:,2);
rank_within_scenario = zeros(size(prcc_values));
for k = 1:numel(scenario_names)
    rank_within_scenario(:,k) = tiedrank(-abs(prcc_values(:,k)));
end
ns = numel(scenario_names);
result_table = table(repmat((1:d)',ns,1),repmat(design.parameter_names,ns,1), ...
    repmat(design.parameter_symbols,ns,1),repelem(scenario_names',d), ...
    repmat(design.reference_values,ns,1),repmat(design.lower_bound,ns,1), ...
    repmat(design.upper_bound,ns,1),prcc_values(:),abs(prcc_values(:)), ...
    ci_lower(:),ci_upper(:),rank_within_scenario(:),'VariableNames', ...
    {'parameter_index','parameter','parameter_symbol','scenario','P01_reference_value', ...
     'lower_bound','upper_bound','PRCC','abs_PRCC','CI95_lower','CI95_upper','rank_within_scenario'});
save(fullfile(output_dir,'p01_prcc_results.mat'),'prcc_values','ci_lower','ci_upper', ...
    'result_table');
writetable(result_table,fullfile(output_dir,'p01_prcc_sensitivity.csv'));
fprintf('PRCC: saved all %d parameters x %d scenarios and 95%% intervals.\n',d,ns);
end

function rho = compute_prcc(X,Y)
% Each inverse rank-correlation matrix contains 76 inputs and ONE outcome.
% A coefficient therefore adjusts for the other 75 input ranks.
n = size(X,1); d = size(X,2);
rx = zscore(tiedrank(X),0,1);
ry = zscore(tiedrank(Y),0,1);
rxx = (rx'*rx)/(n-1);
rxy = (rx'*ry)/(n-1);
rho = zeros(d,size(Y,2));
for k = 1:size(Y,2)
    R = [rxx,rxy(:,k);rxy(:,k)',1];
    assert(rcond(R)>1e-12,'Rank correlation matrix is singular.');
    precision = R\eye(d+1);
    rho(:,k) = -precision(1:d,end)./sqrt(diag(precision(1:d,1:d))*precision(end,end));
end
assert(all(isfinite(rho) & abs(rho)<=1+1e-12,'all'),'Invalid PRCC.');
end
