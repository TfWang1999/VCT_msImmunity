function [samples, unit_lhs, metadata] = generate_lhs_parameter_matrix(sample_size, random_seed)
%GENERATE_LHS_PARAMETER_MATRIX Draw stratified samples from the Beta marginals.
% Immune scenario: lower T-cell influx/expansion and pL/ph, lower p2,
% and higher IL-12/IL-2 activation thresholds within existing bounds.
% rhoL/rhoh/rho8 remain derived from their underlying molecular parameters.
% The d0/lambda0 ranges are informed by the CRC fits.
% run_virtual_cohort reads the saved patient design directly.

if nargin < 2; random_seed = 12345; end
validateattributes(sample_size, {'numeric'}, {'scalar','integer','positive'});
rng(random_seed, 'twister');

% Active bounds are exported in input_data/vct_parameter_design.csv.
% q_A is the effective tumor-site concentration increment per anti-PD-1 dose (nM).
spec = [ ... % lower, upper, beta alpha, beta beta
 .25 .42 4 3; 5e10 1e11 6 4; 1e-9 1e-8 2 6; 5e-9 5e-8 6 4;
 10 30 4 6; 5e8 5e9 2 4; 5e9 5e10 2 4; 5 200 2 4;
 5e8 5e9 6 4; .15 2 2 4; 5 200 2 4; 5e8 5e9 6 4;
 .15 2 2 4; .2 1.5 2 2; .1 10 2 4; 5e7 5e8 4 6;
 5e9 5e10 2 4; .1 .8 2 3; .01 .08 2 2; 5 30 2 2;
 2e4 3e4 4 2; 2e3 1e4 4 2; 5 30 2 2; 1e3 3e3 4 2;
 200 800 2 4; 1 100 2 3; 1 100 2 3; 1 100 2 3;
 10 50 2 2; 1 100 2 3; 1 100 2 3; 10 50 2 2;
 1 5 2 2; 2e4 3e4 4 2; 2e3 1e4 4 2; 200 400 2 2;
 1e3 4e3 2 2; 200 400 4 2; .5 2 2 2; 1e3 4e3 4 2;
 400 2e3 4 2; 200 1e3 4 2; 1 100 4 2; 1 100 4 2;
 50 200 4 2; 1 100 4 2; 1 100 2 2; .2 2 2 2;
 1 10 1 2; 1e-9 1e-7 1 4; 1e-9 1e-7 1 6; 1e-9 1e-7 1 6;
 1e-7 1e-6 4 2; 1 100 2 2; .025 .045 2 2; .1 .35 2 3;
 .1 .4 2 4; .625 3.125 2 4; .1 1 2 3; .05 .15 2 3;
 .01 .25 2 4; .01 .1 2 3; .1 1 2 3; 1 30 2 3;
 5 60 2 3; 1e3 1e6 1 6; .08 1 4 2; .01 .3 4 2;
 .008 .2 1 4; .001 .1 1 4; .003 .3 1 4; .02 .8 1 4;
 1e-8 1e-6 1 6; 1e-8 1e-6 1 6; 1e-10 1e-8 1 6; 1e-7 1e-6 4 2];

parameter_names = {
 'TSM.r_C','TSM.C_max','TSM.eta_h','TSM.eta_8','TSM.lambda_DC','TSM.D_0','TSM.K_DC', ...
 'TSM.lambda_T8','TSM.T_80','TSM.T8h','TSM.lambda_Th','TSM.T_40','TSM.Thh','TSM.lambda_Tr', ...
 'TSM.lambda_M1','TSM.M0','TSM.K_MC','TSM.sigma1','TSM.sigma2','par_cd4.f_1','par_cd4.s_I', ...
 'par_cd4.K1','par_cd4.f_2','par_cd4.s_T','par_cd4.K2','par_Th.a_IFNg_0','par_Th.a_IL2_0', ...
 'par_Th.a_TNFa_0','par_Th.K_Tb','par_Tr.a_TGFb_0','par_Tr.a_IL10_0','par_Tr.K_FOXP', ...
 'par_Mac.f_1','par_Mac.s_I','par_Mac.K1','par_Mac.f_2','par_Mac.s_T','par_Mac.K2', ...
 'par_Mac.f_3','par_Mac.s_tgfbr','par_Mac.s_10','par_Mac.K3','par_M1.b_IL12_0', ...
 'par_M1.b_TNFa_0','par_M1.Km1','par_M2.b_IL10_0','par_M2.b_TGFb_0','par_M2.Km2', ...
 'CY.p_1','CY.p_8','CY.p_h','CY.p_L','CY.K_P','TSM.epC','TSM.a','par_DNA.d_0', ...
 'par_DNA.lambda_0','par_DNA.lambda_P','par_DNA.lambda_F','par_DNA.rho_F','par_DNA.k_b', ...
 'par_DNA.alpha_A','par_DNA.K_A_protein','par_DNA.K_D','PK.q_A','CY.p_2','CY.K_IL12', ...
 'CY.K_IL2','CY.K_IL10','CY.K_IFNg','CY.K_TNFa','CY.K_TGFb','CY.lambda_DC_IL12', ...
 'CY.lambda_C_TGFb','CY.lambda_C_IL10','CY.lambda_T8_IFNg'};

if size(spec,1) ~= 76 || numel(parameter_names) ~= 76
    error('generate_lhs_parameter_matrix:Specification', 'Parameter specification must contain 76 rows.');
end

reference = init_global_params();
reference_values = zeros(1,76);
sampling_scale = repmat({'linear'},1,76);
for parameter_index = 1:76
    fields = strsplit(parameter_names{parameter_index},'.');
    reference_values(parameter_index) = reference.(fields{1}).(fields{2});
end
analysis_dir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
crc_file = fullfile(analysis_dir,'CRC_fit','results','crc_analysis_results.mat');
% lambda_P/lambda_F/rho_F retain original relative heterogeneity after the
% common CRC coupling changes 50->1.5625, 1->0.5 and 0.2->0.1, respectively.

unit_lhs = zeros(sample_size, 76);
samples = zeros(sample_size, 76);
for parameter_index = 1:76
    stratified = ((0:sample_size-1)' + rand(sample_size,1)) ./ sample_size;
    unit_lhs(:,parameter_index) = stratified(randperm(sample_size));
    beta_quantile = betainv(unit_lhs(:,parameter_index), spec(parameter_index,3), spec(parameter_index,4));
    bounds = spec(parameter_index,1:2);
    if strcmp(sampling_scale{parameter_index},'log10')
        samples(:,parameter_index) = 10.^(log10(bounds(1)) + diff(log10(bounds))*beta_quantile);
    else
        samples(:,parameter_index) = bounds(1) + diff(bounds)*beta_quantile;
    end
end

metadata.parameter_names = parameter_names;
metadata.lower_bound = spec(:,1)';
metadata.upper_bound = spec(:,2)';
metadata.beta_alpha = spec(:,3)';
metadata.beta_beta = spec(:,4)';
metadata.random_seed = random_seed;
metadata.reference_values = reference_values;
metadata.sampling_scale = sampling_scale;
metadata.crc_result_file = crc_file;
metadata.variable_parameter_count = sum(spec(:,2)>spec(:,1));
metadata.method = 'Stratified probability sampling with bounded Beta marginals';
metadata.prior_basis = 'Exploratory immune scenario: lambda_T8/lambda_Th/T8h/Thh Beta(2,4) and ph/pL Beta(1,6), within unchanged bounds; rhoL/rhoh/rho8 remain derived. Retain p2 Beta(1,6) on [1e3,1e6], K_IL12 Beta(4,2) on [0.08,1], K_IL2 Beta(4,2) on [0.01,0.30], and other population marginals with reduced p1 and CRC-informed treatment ranges. All marginals are linear. Lower checkpoint coefficients can oppose the effect of lower T-cell influx/expansion. These are simulation priors, not estimated population distributions; no AUC objective or patient selection is used.';
metadata.r_C_prior_basis = 'Beta(4,3) on [0.25,0.42] day^-1.';
end
