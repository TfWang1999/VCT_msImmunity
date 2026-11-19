function params = apply_mouse_manual_parameters(case_data)
%APPLY_MOUSE_MANUAL_PARAMETERS Reproduce submitted mouse-model parameter sets.

params = init_global_params();
params.TSM.r_C = case_data.manual.r_C;
params.TSM.C_max = case_data.manual.C_max;
params.TSM.T_40 = case_data.manual.T_40;
params.TSM.T_80 = case_data.manual.T_80;
params.CY.p_1 = case_data.manual.p_1;
params.CY.p_2 = case_data.manual.p_2;
params.TSM.epC = case_data.manual.epC;
params.PK.q_A = case_data.manual.q_A; % Effective tumor-site increment per anti-PD-1 dose, nM
params.PK.A_half_days=223/24; % Boche et al. 2022, Fig. 2B; mouse reference assumption
params.TSM.eta_h = 1e-8;
params.TSM.eta_8 = 4e-8;
end

