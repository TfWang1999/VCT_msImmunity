function S = generate_lhs_parameter_matrix(sample_size)
%% FUNCTION NAME: generate_lhs_parameter_matrix
%
% DESCRIPTION: 
%     Generates a Latin Hypercube Sampling (LHS) matrix for 76 key parameters 
%     spanning molecular, microenvironmental, and tissue scales. 
%     Each row represents a unique "Virtual Patient" parameter set, where 
%     variability is modeled using Beta distributions to capture 
%     population-level heterogeneity.[3, 6, 5]
%
% INPUTS:
%     sample_size - (int) Number of virtual patients to generate (e.g., 1000)
%
% OUTPUTS:
%     S - (Matrix) Matrix of size [sample_size x 76] containing sampled parameters
%
% DEPENDENCIES:
%     This function requires 'Beta_Matrix.m' to perform the underlying 
%     Beta distribution transformations .
%
%---------------------------------------------------------

    % Total number of parameters identified as sources of heterogeneity
    total_params = 76;
    S = zeros(sample_size, total_params);
    
    %% --- Section 1: Tissue Scale Parameters ---
    
    % Cancer Cell Dynamics (Equation 1: dC/dt)
    S(:,1) = Beta_Matrix(0.4, 0.6, 3, 3, sample_size);     % TSM.r_C: Proliferation rate of tumor cells
    S(:,2) = Beta_Matrix(5e10, 1e11, 6, 4, sample_size);   % TSM.C_max: Carrying capacity of tumor cells
    S(:,3) = Beta_Matrix(1e-9, 1e-8, 2, 6, sample_size);   % TSM.eta_h: Killing rate of tumor cell by Th
    S(:,4) = Beta_Matrix(5e-9, 5e-8, 6, 4, sample_size);   % TSM.eta_8: Killing rate of tumor cell by T8
    
    % Dendritic Cell Dynamics (Equation 2: dD/dt)
    S(:,5) = Beta_Matrix(10, 30, 4, 6, sample_size);       % TSM.lambda_DC: Activation rate of DCs
    S(:,6) = Beta_Matrix(5e8, 5e9, 2, 4, sample_size);     % TSM.D_0: Initial immature DCs
    S(:,7) = Beta_Matrix(5e9, 5e10, 2, 4, sample_size);    % TSM.K_DC: Half-saturation constant for DC activation
    
    % CD8+ T Cell Dynamics (Equation 3: dT8/dt)
    S(:,8) = Beta_Matrix(5, 200, 2, 2, sample_size);       % TSM.lambda_T8: Activation rate of CD8 T cells
    S(:,9) = Beta_Matrix(5e8, 5e9, 6, 4, sample_size);     % TSM.T80: Initial naive CD8 T cells
    S(:,10) = Beta_Matrix(0.15, 2, 2, 2, sample_size);     % TSM.T8h: Proliferation rate of CD8 T cells
    
    % Helper T Cell Dynamics (Equation 4: dTh/dt)
    S(:,11) = Beta_Matrix(5, 200, 2, 2, sample_size);      % TSM.lambda_Th: Activation rate of Th cells
    S(:,12) = Beta_Matrix(5e8, 5e9, 6, 4, sample_size);    % TSM.T_40: Initial naive CD4 T cells
    S(:,13) = Beta_Matrix(0.15, 2, 2, 2, sample_size);     % TSM.Thh: Proliferation rate of Th cells
    
    % Regulatory T Cell (Treg) Dynamics
    S(:,14) = Beta_Matrix(0.2, 1.5, 2, 2, sample_size);    % lambda_Tr: Activation rate of Tregs
    
    % Tumor-Associated Macrophage (TAM) Dynamics (Equation 5/6)
    S(:,15) = Beta_Matrix(0.1, 10, 2, 4, sample_size);     % TSM.lambda_M: Recruitment rate of macrophages
    S(:,16) = Beta_Matrix(5e7, 5e8, 4, 6, sample_size);    % TSM.M0: Initial monocyte/macrophage baseline
    S(:,17) = Beta_Matrix(5e9, 5e10, 2, 4, sample_size);   % TSM.K_MC: Half-saturation constant for TAM recruitment
    S(:,18) = Beta_Matrix(0.1, 0.8, 2, 3, sample_size);    % TSM.sigma1: Phenotype change rate (M1 --> M2)
    S(:,19) = Beta_Matrix(0.01, 0.08, 2, 2, sample_size);  % TSM.sigma2: Phenotype change rate (M2 --> M1)
    
    %% --- Section 2: Intracellular CD4+ T Cell Module ---
    
    % IFN-gamma binding and receptor dynamics
    S(:,20) = Beta_Matrix(5, 30, 2, 2, sample_size);       % par_cd4.f_1: T-bet production rate
    S(:,21) = Beta_Matrix(2e4, 3e4, 4, 2, sample_size);    % par_cd4.s_I: IFNgR density on CD4+ surface
    S(:,22) = Beta_Matrix(2e3, 1e4, 4, 2, sample_size);    % par_cd4.K1: Half-saturation of IFNg-R
    
    % TGF-beta binding and receptor dynamics
    S(:,23) = Beta_Matrix(5, 30, 2, 2, sample_size);       % par_cd4.f_2: FOXP3 production rate
    S(:,24) = Beta_Matrix(1e3, 3e3, 4, 2, sample_size);    % par_cd4.s_T: TGFbR density on CD4+ surface
    S(:,25) = Beta_Matrix(200, 800, 2, 4, sample_size);    % par_cd4.K2: Half-saturation of TGFb-R
    
    % Cytokine baseline production rates
    S(:,26) = Beta_Matrix(1, 100, 2, 3, sample_size);      % par_Th.a_IFNg_0: Baseline IFNg by Th
    S(:,27) = Beta_Matrix(1, 100, 2, 3, sample_size);      % par_Th.a_IL2_0: Baseline IL2 by Th
    S(:,28) = Beta_Matrix(1, 100, 2, 3, sample_size);      % par_Th.a_TNFa_0: Baseline TNFa by Th
    S(:,29) = Beta_Matrix(10, 50, 2, 2, sample_size);      % par_Th.K_Tb: Half-saturation of T-bet
    S(:,30) = Beta_Matrix(1, 100, 2, 3, sample_size);      % par_Tr.a_TGFb_0: Baseline TGFb by Tr
    S(:,31) = Beta_Matrix(1, 100, 2, 3, sample_size);      % par_Tr.a_IL10_0: Baseline IL10 by Tr
    S(:,32) = Beta_Matrix(10, 50, 2, 2, sample_size);      % par_Tr.K_FOXP: Half-saturation of FOXP3
    
    %% --- Section 3: Intracellular Macrophage (TAM) Module ---
    
    S(:,33) = Beta_Matrix(1, 5, 2, 2, sample_size);        % par_Mac.f_1: pSTAT1 production rate
    S(:,34) = Beta_Matrix(2e4, 3e4, 4, 2, sample_size);    % par_Mac.s_I: IFNgR density on TAM surface
    S(:,35) = Beta_Matrix(2e3, 1e4, 4, 2, sample_size);    % par_Mac.K1: Half-saturation of IFNg-R
    S(:,36) = Beta_Matrix(200, 400, 2, 2, sample_size);    % par_Mac.f_2: NFkB production rate
    S(:,37) = Beta_Matrix(1000, 4000, 2, 2, sample_size);  % par_Mac.s_T: TNFaR density on TAM surface
    S(:,38) = Beta_Matrix(200, 400, 4, 2, sample_size);    % par_Mac.K2: Half-saturation of TNFa-R
    S(:,39) = Beta_Matrix(0.5, 2, 2, 2, sample_size);      % par_Mac.f_3: pSTAT3 production rate
    S(:,40) = Beta_Matrix(1000, 4000, 4, 2, sample_size);  % par_Mac.s_tgfbr: TGFbR density on TAM surface
    S(:,41) = Beta_Matrix(400, 2000, 4, 2, sample_size);   % par_Mac.s_10: IL10R density on TAM surface
    S(:,42) = Beta_Matrix(200, 1000, 4, 2, sample_size);   % par_Mac.K3: Half-saturation for IL10-R and TGFb-R
    
    % Polarization phenotypes
    S(:,43) = Beta_Matrix(1, 100, 4, 2, sample_size);      % par_M1.b_IL12_0: Baseline IL12 by M1
    S(:,44) = Beta_Matrix(1, 100, 4, 2, sample_size);      % par_M1.b_TNFa_0: Baseline TNFa by M1
    S(:,45) = Beta_Matrix(50, 200, 4, 2, sample_size);     % par_M1.Km1: Half-saturation Michaelis constant for m1
    S(:,46) = Beta_Matrix(1, 100, 4, 2, sample_size);      % par_M2.b_IL10_0: Baseline IL10 by M2
    S(:,47) = Beta_Matrix(1, 100, 2, 2, sample_size);      % par_M2.b_TGFb_0: Baseline TGFb by M2
    S(:,48) = Beta_Matrix(0.2, 2, 2, 2, sample_size);      % par_M2.Km2: Half-saturation Michaelis constant for m2
    
    %% --- Section 4: Anti-PD-1 Treatment & Therapy Module ---
    
    % PD-1/PD-L1 signaling parameters
    S(:,49) = Beta_Matrix(1, 100, 1, 2, sample_size);      % CY.p_1: PD-1/PD-L1 binding rate constant
    S(:,50) = Beta_Matrix(1e-9, 1e-7, 1, 4, sample_size);  % CY.p_8: PD-1 expression rate on T8
    S(:,51) = Beta_Matrix(1e-9, 1e-7, 1, 4, sample_size);  % CY.p_h: PD-1 expression rate on Th
    S(:,52) = Beta_Matrix(1e-9, 1e-7, 1, 4, sample_size);  % CY.p_L: PD-L1 expression rate
    S(:,53) = Beta_Matrix(1e-7, 1e-6, 4, 2, sample_size);  % CY.K_P: Inhibitory constant for PD1-PDL1 complex
    S(:,54) = Beta_Matrix(1, 100, 2, 2, sample_size);      % TSM.ep_C: PD-L1 amplification in tumor vs T cells
    
    % Radiotherapy and Chemotherapy (DNA Damage Model)
    S(:,55) = Beta_Matrix(0.025, 0.045, 2, 2, sample_size); % TSM.a: LQ model constant for radiotherapy
    S(:,56) = Beta_Matrix(0.1, 0.25, 2, 2, sample_size);   % par_DNA.d_0: Induced death rate by DNA damage
    S(:,57) = Beta_Matrix(0.15, 0.3, 2, 2, sample_size);   % par_DNA.lambda_0: Growth inhibition by DNA damage
    S(:,58) = Beta_Matrix(20, 100, 2, 4, sample_size);     % par_DNA.lambda_P: Damage coefficient (Oxaliplatin)
    S(:,59) = Beta_Matrix(0.2, 2, 2, 3, sample_size);      % par_DNA.lambda_F: Damage coefficient (Capecitabine)
    S(:,60) = Beta_Matrix(0.1, 0.3, 2, 3, sample_size);    % par_DNA.rho_F: 5-FU repair inhibition coefficient
    S(:,61) = Beta_Matrix(0.01, 0.25, 2, 4, sample_size);  % par_DNA.k_b: Recovery rate of repair proteins
    S(:,62) = Beta_Matrix(0.01, 0.1, 2, 3, sample_size);   % par_DNA.alpha_A: Pro-apoptotic protein activation rate
    S(:,63) = Beta_Matrix(0.1, 1, 2, 3, sample_size);      % par_DNA.K_A_protein: Half-saturation of pro-apoptotic proteins
    S(:,64) = Beta_Matrix(1, 30, 2, 3, sample_size);       % par_DNA.K_D: Half-saturation of DNA damage sites
    S(:,65) = Beta_Matrix(10, 100, 2, 3, sample_size);     % CY.A: Anti-PD1 dose concentration (mg/L)
    S(:,66) = Beta_Matrix(1e3, 1e6, 2, 6, sample_size);    % CY.p_2: Binding constant for anti-PD1/PD1
    
    %% --- Section 5: Extracellular Microenvironmental Scale ---
    
    % Cytokine half-saturation Michaelis constants
    S(:,67) = Beta_Matrix(0.08, 1, 1, 2, sample_size);     % CY.K_IL12
    S(:,68) = Beta_Matrix(0.01, 0.3, 1, 4, sample_size);    % CY.K_IL2
    S(:,69) = Beta_Matrix(0.008, 0.2, 1, 4, sample_size);   % CY.K_IL10
    S(:,70) = Beta_Matrix(0.001, 0.1, 1, 4, sample_size);   % CY.K_IFNg
    S(:,71) = Beta_Matrix(0.003, 0.3, 1, 4, sample_size);   % CY.K_TNFa
    S(:,72) = Beta_Matrix(0.02, 0.8, 1, 4, sample_size);    % CY.K_TGFb
    
    % Cytokine production rates by specific subsets
    S(:,73) = Beta_Matrix(1e-8, 1e-6, 1, 6, sample_size);  % CY.lambda_DC_IL12: Production of IL12 by DCs
    S(:,74) = Beta_Matrix(1e-8, 1e-6, 1, 6, sample_size);  % CY.lambda_C_TGFb: Production of TGFb by tumor cells
    S(:,75) = Beta_Matrix(1e-10, 1e-8, 1, 6, sample_size); % CY.lambda_C_IL10: Production of IL10 by tumor cells
    S(:,76) = Beta_Matrix(1e-7, 1e-6, 4, 2, sample_size);  % CY.lambda_T8_IFNg: Production of IFNg by CD8 T cells
end

function y=Beta_Matrix(a,b,c,d,e)

y=a+(b-a)*betarnd(c,d,[1,e]);

end