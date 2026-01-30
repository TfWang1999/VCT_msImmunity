function params = init_global_params()
%% FUNCTION NAME: init_global_params
%
% DESCRIPTION: 
%     Initializes global parameters for the multiscale tumor-immune model,
%     including cytokine dynamics, signaling pathways, and tissue-scale constants.
%
% OUTPUTS:
%     params - A structure containing nested structs for different biological scales:
%             .CY: Cytokine microenvironment parameters
%             .par_cd4: Intracellular signaling for CD4+ T cells
%             .par_Th /.par_Tr: T-helper and regulatory T cell parameters
%             .par_Mac /.par_M1 /.par_M2: Macrophage polarization parameters
%             .TSM: Tissue-scale dynamics and population constants
%             .par_DNA: DNA damage and repair parameters for chemotherapy
%

%---------------------------------------------------------

%% CY: Cytokine Microenvironment Parameters
CY = struct();
CY.F_ave = 0.367;
CY.P_ave = 0.0236;
CY.d_IL12 = 3.3; % unit: /day
CY.K_IL12 = 0.8;
CY.d_IL2 = 23.8; 
CY.K_IL2 = 0.237;
CY.d_IL10 = 5.94;
CY.K_IL10 = 0.0875;
CY.d_IFNg = 34.70;
CY.K_IFNg = 0.018;
CY.d_TNFa = 55.00;
CY.K_TNFa = 0.03;
CY.d_TGFb = 198.00;
CY.K_TGFb = 0.206;
CY.lambda_DC_IL12 = 9e-7; 
CY.lambda_C_TGFb = 1.1e-7; 
CY.lambda_C_IL10 = 1.3e-10; 
CY.lambda_T8_IFNg = 2.5e-7;

%% CD4: CD4+ T Cell model Parameters (Molecular Scale)
par_cd4 = struct();
par_cd4.p_I = 20;   % Binding rate of [IFNg] 
par_cd4.p_T = 10;   % Binding rate of [TGFb]
par_cd4.K1 = 3000;
par_cd4.K2 = 300;
par_cd4.s_I = 27000;  % Maximum binding of
par_cd4.s_T = 1500;   % Maximum binding of
par_cd4.d_ifnr = 6;   % Internalization/degradation rate of
par_cd4.d_tgfr = 6;   % Internalization/degradation rate of
par_cd4.f_1 = 20;     % Generation rate of
par_cd4.f_2 = 20;     % Generation rate of [FOXP3]
par_cd4.k_F = 30;     % Inhibition coefficient of [FOXP3] on
par_cd4.k_T = 30;     % Inhibition coefficient of on [FOXP3]
par_cd4.d_F = 0.15;   % Degradation rate of [FOXP3]
par_cd4.d_T = 0.15;   % Degradation rate of

%% Th and Tr: T Cell Subset Functional Parameters
% Th (Helper T cells)
par_Th = struct();
% baseline production rate of cytokines by Th
par_Th.a_IFNg_0 = 100; 
par_Th.a_IL2_0 = 20;
par_Th.a_TNFa_0 = 20;
% The half-saturation constant
par_Th.K_Tb = 30;
% proportion of Th cells that secrete cytokines
par_Th.IFNg_proportion = 0.2; 
par_Th.IL2_proportion = 0.2;
par_Th.TNFa_proportion = 0.2;

% Tr (Regulatory T cells)
par_Tr = struct();
% baseline production rate of cytokines by Tr
par_Tr.a_TGFb_0 = 10;
par_Tr.a_IL10_0 = 10;
% The half-saturation constant
par_Tr.K_FOXP = 30;
% proportion of Tr cells that secrete cytokines
par_Tr.TGFb_proportion = 0.1;
par_Tr.IL10_proportion = 0.1;

%% Mac: TAM model Parameters (Molecular Scale)
par_Mac = struct();
% Binding rates
par_Mac.p_g = 20;   % [IFNg] binding rate /nM h
par_Mac.p_a = 63;   % [TNFa] binding rate
par_Mac.p_t = 20;   % [TGFb] binding rate
par_Mac.p_10 = 40;  % [IL10] binding rate
% Constants
par_Mac.K1 = 5000;
par_Mac.K2 = 300;
par_Mac.K3 = 500;
% Maximum binding capacities
par_Mac.s_I = 27000;    
par_Mac.s_T = 1500;     
par_Mac.s_10 = 800;     
par_Mac.s_tgfbr = 1500;
% Internalization degradation rates
par_Mac.d_ifnr = 6/60; 
par_Mac.d_tgfr = 6/60; 
par_Mac.d_il10r = 6/60;
par_Mac.d_tnfr = 6/60; 
% Generation rates
par_Mac.f_1 = 4;  % pSTAT1   
par_Mac.f_2 = 300; % NFkB
par_Mac.f_3 = 1;   % pSTAT3   
% The half-saturation constant
par_Mac.k_S = 10;    
par_Mac.d_s1 = 0.15;  
par_Mac.d_nf = 5;  
par_Mac.d_s3 = 0.15;  

%% M1 and M2: Polarization Parameters
% M1 Macrophages
par_M1 = struct();
% Half-saturation constant
par_M1.Km1 = 100; 
% baseline production rate of cytokines by M1
par_M1.b_IL12_0 = 20;
par_M1.b_TNFa_0 = 20;
par_M1.b_IFNg_0 = 20;
% proportion of M1 that secrete cytokines
par_M1.IL12_proportion = 0.2;
par_M1.TNFa_proportion = 0.2;
par_M1.IFNg_proportion = 0.2;

% M2 Macrophages
par_M2 = struct();
par_M2.Km2 = 1; % Half-saturation parameter
% baseline production rate of cytokines by M2
par_M2.b_IL10_0 = 20;
par_M2.b_TGFb_0 = 20;
% proportion of M1 that secrete cytokines
par_M2.IL10_proportion = 0.2;
par_M2.TGFb_proportion = 0.2;

%% TSM: Tissue Scale Dynamics 
TSM = struct();
% Equation 1: Tumor growth and killing
TSM.r_C = 0.49;    % Tumor growth rate (/day)
TSM.C_max = 3e10;  % Carrying capacity (max tumor cell count)
TSM.d_C = 0.17;    % Natural tumor death rate
TSM.eta_h = 1e-8;  % Killing efficiency of Th cells
TSM.eta_8 = 4e-8;  % Killing efficiency of CD8+ T cells

% Equation 2: Dendritic Cells (DCs)
TSM.lambda_DC = 20;           % DC activation rate
TSM.D_0 = 1.94e7;             % Basal number of naive DCs
TSM.K_DC = TSM.C_max / 2;     % Half-saturation parameter for activation
TSM.d_D = 0.2;                % DC death rate

% Equation 3: CD8+ T Cells (T8)
TSM.lambda_T8 = 24.9; % Activation rate of CD8+ T cells
TSM.T_80 = 9.39e8;    % Basal number of CD8+ T cells
TSM.T8h = 0.25;       % Proliferation rate of T8 cells
TSM.ep8 = 1;          % Relative PD-L1 expression: T8 vs Th
TSM.epC = 50;         % Relative PD-L1 expression: C vs Th
TSM.d_8 = 0.18;       % Death rate of T8 cells

% Equation 4: Helper T Cells (Th)
TSM.lambda_Th = 27.96; % Activation rate of Th cells
TSM.T_40 = 1.88e9;     % Basal number of CD4+ T cells
TSM.Thh = 0.25;        % Proliferation rate of Th cells
TSM.d_h = 0.197;       % Death rate of Th cells

% Equation 5: Regulatory T Cells (Tr)
TSM.lambda_Tr = 0.35; % Recruitment rate of Tregs
TSM.d_r = 0.2;        % Death rate of Tregs

% Equation 6 & 7: Macrophages (M1 & M2)
TSM.lambda_M1 = 1.65;
TSM.M0 = 1e8;              % Basal number of macrophages
TSM.K_MC = TSM.C_max / 2;  % Recruitment half-saturation
TSM.sigma1 = 0.2;          % Phenotype switch rate: M1 --> M2
TSM.sigma2 = 0.02;         % Phenotype switch rate: M2 --> M1
TSM.d_M1 = 0.15;           % Death rate of M1
TSM.d_M2 = 0.15;           % Death rate of M2

%% Therapy Parameters
% RT parameters
TSM.a = 0.03;   
TSM.b = TSM.a / 10;
TSM.R_D = 0.6;    TSM.R_8 = 0.6;
TSM.R_h = 0.6;    TSM.R_r = 0.7;
TSM.R_1 = 0.3;    TSM.R_2 = 0.3;
TSM.Dose = 0; % Dosage placeholder

% Chemotherapy (DNA Damage Dynamics)
par_DNA.F_ave = 0.367;
par_DNA.P_ave = 0.0236;
par_DNA.lambda_P = 50;
par_DNA.lambda_F = 1;
par_DNA.k_b = 0.05;
par_DNA.alpha_A = 0.05;
par_DNA.K_D = 10;
par_DNA.mu_A = 0.035;
par_DNA.rho_F = 0.2;
par_DNA.K_A_protein = 0.5;
par_DNA.D_damage = (par_DNA.lambda_F*par_DNA.F_ave + par_DNA.lambda_P*par_DNA.P_ave) * (1+par_DNA.rho_F*par_DNA.F_ave) / par_DNA.k_b;
par_DNA.A_protein = par_DNA.alpha_A * par_DNA.D_damage/(par_DNA.D_damage+par_DNA.K_D) / par_DNA.mu_A;
par_DNA.d_0 = 0.25;
par_DNA.lambda_0 = 0.3;

TSM.k_D = 0.6;
TSM.k_8 = 0.6;     TSM.k_h = 0.6;
TSM.k_r = 0.6;     TSM.k_1 = 0.6;
TSM.k_2 = 0.6;   

% Immunotherapy (PD-1 / PD-L1 Axis)
CY.p_1 = 50;   % PD-1 PD-L1 binding constant
CY.p_2 = 1e6;  % Anti-PD-1 PD-L1 binding constant
CY.p_8 = 1e-7; % PD-1 expression rate of T8
CY.p_h = 1e-7; % PD-1 expression rate of Th
CY.p_L = 3e-9; % PD-L1 expression rate
CY.K_P = 1e-7; % Inhibitory constant
CY.A = 50;     % Anti-PD-1 concentration (30~100 mg/L)

%% Integration into params output
params.par_DNA = par_DNA;
params.CY = CY;
params.par_cd4 = par_cd4;
params.par_Th = par_Th;
params.par_Tr = par_Tr;
params.par_Mac = par_Mac;
params.par_M1 = par_M1;
params.par_M2 = par_M2;
params.TSM = TSM;

end