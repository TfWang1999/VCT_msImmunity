function params = apply_crc_manual_parameters(case_data)
%APPLY_CRC_MANUAL_PARAMETERS Apply CRC starting values and shared fixed inputs.

params = init_global_params();
params.TSM.r_C = case_data.manual.r_C;
params.TSM.C_max = case_data.manual.C_max;
params.PK.q_A = case_data.manual.q_A;
params.par_DNA.d_0 = case_data.manual.d_0;
params.par_DNA.lambda_0 = case_data.manual.lambda_0;
params.TSM.a = case_data.manual.a;

% Shared chemotherapy coupling coefficients.
params.par_DNA.lambda_F = 0.5;
params.par_DNA.lambda_P = 1.5625;
params.par_DNA.rho_F = 0.1;
params.par_DNA.D_damage = (params.par_DNA.lambda_F * params.par_DNA.F_ave + ...
    params.par_DNA.lambda_P * params.par_DNA.P_ave) * ...
    (1 + params.par_DNA.rho_F * params.par_DNA.F_ave) / params.par_DNA.k_b;
params.par_DNA.A_protein = params.par_DNA.alpha_A * params.par_DNA.D_damage / ...
    (params.par_DNA.D_damage + params.par_DNA.K_D) / params.par_DNA.mu_A;

% Patient-specific checkpoint parameters; remaining inputs are common.
params.CY.p_1 = case_data.manual.p_1;
params.CY.p_2 = case_data.manual.p_2;
params.CY.p_8 = 1e-7;
params.CY.p_h = 1e-7;
params.CY.p_L = 3e-9;
params.CY.K_P = 1e-7;
params.TSM.eta_h = 1e-8;
params.TSM.eta_8 = 4e-8;
end

