function output_root = run_virtual_cohort(workers)
%RUN_VIRTUAL_COHORT Replay the saved cohort; never draw a new population.
% workers=0 (default) runs serially. A positive value enables parfor workers.
if nargin < 1; workers = 0; end
module_dir = fileparts(mfilename('fullpath'));
analysis_dir = fileparts(module_dir);
addpath(genpath(fullfile(analysis_dir,'src')));
design = load(fullfile(module_dir,'input_data','VP_parameter_matrix.mat'));
output_root = fullfile(module_dir,'results','cohort');
folder_names = {'VP_Xelox','VP_Xelox_IR','VP_Xelox_antipd1','VP_Xelox_IR_antipd1'};
therapy_files = {'Therapy_xelox.mat','Therapy_xelox_IR.mat', ...
    'Therapy_xelox_antiPD1.mat','Therapy_xelox_antiPD1_IR.mat'};
samples = design.samples;
initial_fraction = design.initial_tumor_fraction;
for k = 1:4
    regimen_dir = fullfile(output_root,folder_names{k});
    if ~exist(regimen_dir,'dir'); mkdir(regimen_dir); end
    therapy = load_single_variable(fullfile(module_dir,'input_data',therapy_files{k}));
    parfor (i = 1:size(samples,1), workers)
        params = apply_virtual_patient_parameters(samples(i,:));
        experiment = build_virtual_patient_experiment(therapy,params,initial_fraction(i));
        [time,volume,~,~,params] = run_simulation(experiment,params);
        % Only the trajectory and parameter snapshot needed by the analysis.
        x = struct('t_sol',time,'Vtumor_list',volume,'params',params);
        save(fullfile(regimen_dir,sprintf('VP%d.mat',i)),'-fromstruct',struct('x',x));
    end
    fprintf('%s: %d saved patients reproduced.\n',folder_names{k},size(samples,1));
end
end
