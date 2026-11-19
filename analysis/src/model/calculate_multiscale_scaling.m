function [Dim_TSM,new_params] = calculate_multiscale_scaling(V_tumor, params, therapy)
%% FUNCTION NAME: calculate_multiscale_scaling
%
% DESCRIPTION: 
%     1. Calculates molecular-scale coupling parameters dependent on tumor volume.
%     2. Generates the dimensionless parameter set required for ODE simulation.
%     This function acts as the bridge between raw physical parameters and the solver.
%
% INPUTS:
%     V_tumor - Current macroscopic tumor volume (cm^3)
%     params  - Global parameter structure (from init_global_params)
%
% OUTPUTS:
%     Dim_TSM    - Struct containing scaled, dimensionless parameters for the tissue-scale ODE system
%     new_params - Updated physical parameter structure including derived coupling terms
%
%---------------------------------------------------------

    % Extract sub-structures
    par_DNA = params.par_DNA;
    CY = params.CY;
    par_cd4 = params.par_cd4;
    par_Th = params.par_Th;
    par_Tr = params.par_Tr;
    par_Mac = params.par_Mac;
    par_M1 = params.par_M1;
    par_M2 = params.par_M2;
    TSM = params.TSM;

    %% Fixed Physical Constants and Conversion Factors
    M_IFNg = 17;        % kDa
    M_IL2 = 15.5;       % kDa
    M_TNFa = 17;        % kDa
    M_TGFb = 12.5;      % kDa
    M_IL10 = 18.5;      % kDa
    M_IL12 = (35+40)/2; % kDa
    NA = 6.02214076e23; % Avogadro's number
    
    % Time and Unit conversion (day -> sec, kDa -> g/mol, g -> ng)
    Tc = 24*60*60*1000*1e9; 

    %% Helper T Cell (Th) Signaling Coupling
    par_Th.Tbet_ave = par_cd4.f_1 * par_cd4.s_I / (2 * (par_cd4.s_I + par_cd4.K1) * par_cd4.d_T);
    par_Th.a_IFNg_ave = par_Th.a_IFNg_0 * par_Th.Tbet_ave / (par_Th.Tbet_ave + par_Th.K_Tb);
    par_Th.a_IL2_ave = par_Th.a_IL2_0 * par_Th.Tbet_ave / (par_Th.Tbet_ave + par_Th.K_Tb);
    par_Th.a_TNFa_ave = par_Th.a_TNFa_0 * par_Th.Tbet_ave / (par_Th.Tbet_ave + par_Th.K_Tb);
    
    % Volume-dependent production rates for Th cytokines
    par_Th.lambda_IFNg = par_Th.a_IFNg_ave * M_IFNg * Tc * par_Th.IFNg_proportion / NA / V_tumor;
    par_Th.lambda_IL2 = par_Th.a_IL2_ave * M_IL2 * Tc * par_Th.IL2_proportion / NA / V_tumor;
    par_Th.lambda_TNFa = par_Th.a_TNFa_ave * M_TNFa * Tc * par_Th.TNFa_proportion / NA / V_tumor;

    %% Regulatory T Cell (Treg) Signaling Coupling
    par_Tr.FOXP_ave = par_cd4.f_2 * par_cd4.s_T / (2 * (par_cd4.s_T + par_cd4.K2) * par_cd4.d_F);
    par_Tr.a_TGFb_ave = par_Tr.a_TGFb_0 * par_Tr.FOXP_ave / (par_Tr.FOXP_ave + par_Tr.K_FOXP);
    par_Tr.a_IL10_ave = par_Tr.a_IL10_0 * par_Tr.FOXP_ave / (par_Tr.FOXP_ave + par_Tr.K_FOXP);
    
    % Volume-dependent production rates for Treg cytokines
    par_Tr.lambda_TGFb = par_Tr.a_TGFb_ave * M_TGFb * Tc * par_Tr.TGFb_proportion / NA / V_tumor;
    par_Tr.lambda_IL10 = par_Tr.a_IL10_ave * M_IL10 * Tc * par_Tr.IL10_proportion / NA / V_tumor;

    %% M1 Macrophage Polarization Coupling
    par_M1.m1_ave = (par_Mac.f_1 * par_Mac.f_2 * par_Mac.s_I * par_Mac.s_T) /...
        (4 * (par_Mac.s_I + par_Mac.K1) * (par_Mac.s_T + par_Mac.K2) * par_Mac.d_nf * par_Mac.d_s1);
    par_M1.b_IL12_ave = par_M1.b_IL12_0 * (par_M1.m1_ave / (par_M1.m1_ave + par_M1.Km1));
    par_M1.b_TNFa_ave = par_M1.b_TNFa_0 * (par_M1.m1_ave / (par_M1.m1_ave + par_M1.Km1));
    
    % Volume-dependent production rates for M1 cytokines
    par_M1.lambda_IL12 = par_M1.b_IL12_ave * M_IL12 * Tc * par_M1.IL12_proportion / NA / V_tumor;
    par_M1.lambda_TNFa = par_M1.b_TNFa_ave * M_TNFa * Tc * par_M1.TNFa_proportion / NA / V_tumor;

    %% M2 Macrophage Polarization Coupling
    par_M2.m2_ave = (par_Mac.f_3 * (par_Mac.s_tgfbr + par_Mac.s_10)) /...
        (2 * (par_Mac.s_10 + par_Mac.s_tgfbr + par_Mac.K3) * par_Mac.d_s3);
    par_M2.b_IL10_ave = par_M2.b_IL10_0 * (par_M2.m2_ave / (par_M2.m2_ave + par_M2.Km2));
    par_M2.b_TGFb_ave = par_M2.b_TGFb_0 * (par_M2.m2_ave / (par_M2.m2_ave + par_M2.Km2));
    
    % Volume-dependent production rates for M2 cytokines
    par_M2.lambda_IL10 = par_M2.b_IL10_ave * M_IL10 * Tc * par_M2.IL10_proportion / NA / V_tumor;
    par_M2.lambda_TGFb = par_M2.b_TGFb_ave * M_TGFb * Tc * par_M2.TGFb_proportion / NA / V_tumor;

    %% Tissue-Level Derived Coupling Constants
    TSM.alpha1 = par_M1.lambda_IL12 / CY.d_IL12 / CY.K_IL12;
    TSM.alphaD = CY.lambda_DC_IL12 / CY.d_IL12 / CY.K_IL12;
    TSM.alpha2 = par_M2.lambda_IL10 / CY.d_IL10 / CY.K_IL10 + par_M2.lambda_TGFb / CY.d_TGFb / CY.K_TGFb;
    TSM.alphaC = CY.lambda_C_IL10 / CY.d_IL10 / CY.K_IL10 + CY.lambda_C_TGFb / CY.d_TGFb / CY.K_TGFb;
    TSM.alphaR = par_Tr.lambda_IL10 / CY.d_IL10 / CY.K_IL10 + par_Tr.lambda_TGFb / CY.d_TGFb / CY.K_TGFb;
    TSM.K_h = CY.d_IL2 * CY.K_IL2 / par_Th.lambda_IL2;

    TSM.beta1 = par_M1.lambda_IL12 / CY.d_IL12 / CY.K_IL12;
    TSM.betaD = CY.lambda_DC_IL12 / CY.d_IL12 / CY.K_IL12;
    TSM.beta2 = par_M2.lambda_IL10 / CY.d_IL10 / CY.K_IL10 + par_M2.lambda_TGFb / CY.d_TGFb / CY.K_TGFb;
    TSM.betaC = CY.lambda_C_IL10 / CY.d_IL10 / CY.K_IL10 + CY.lambda_C_TGFb / CY.d_TGFb / CY.K_TGFb;
    TSM.betaR = par_Tr.lambda_IL10 / CY.d_IL10 / CY.K_IL10 + par_Tr.lambda_TGFb / CY.d_TGFb / CY.K_TGFb;

    TSM.gammaC = CY.lambda_C_TGFb / CY.d_TGFb / CY.K_TGFb;
    TSM.gammaR = par_Tr.lambda_TGFb / CY.d_TGFb / CY.K_TGFb;
    TSM.gamma2 = par_M2.lambda_TGFb / CY.d_TGFb / CY.K_TGFb;
    TSM.gammah = par_Th.lambda_IFNg / CY.d_IFNg / CY.K_IFNg;
    TSM.gamma8 = CY.lambda_T8_IFNg / CY.d_IFNg / CY.K_IFNg;

    TSM.delta1 = par_M1.lambda_TNFa / CY.d_TNFa / CY.K_TNFa;
    TSM.delta2 = par_M2.lambda_TGFb / CY.d_TGFb / CY.K_TGFb + par_M2.lambda_IL10 / CY.d_IL10 / CY.K_IL10;
    TSM.deltah = par_Th.lambda_IFNg / CY.d_IFNg / CY.K_IFNg + par_Th.lambda_TNFa / CY.d_TNFa / CY.K_TNFa;
    TSM.delta8 = CY.lambda_T8_IFNg / CY.d_IFNg / CY.K_IFNg;
    TSM.deltaR = par_Tr.lambda_TGFb / CY.d_TGFb / CY.K_TGFb + par_Tr.lambda_IL10 / CY.d_IL10 / CY.K_IL10;
    TSM.deltaC = CY.lambda_C_IL10 / CY.d_IL10 / CY.K_IL10 + CY.lambda_C_TGFb / CY.d_TGFb / CY.K_TGFb;

    %% Therapy Modulation and System PD-1/PD-L1 Logic
    % Chemotherapy model
    par_DNA.D_damage = (par_DNA.lambda_F * par_DNA.F_ave + par_DNA.lambda_P * par_DNA.P_ave) *...
        (1 + par_DNA.rho_F * par_DNA.F_ave) / par_DNA.k_b;
    par_DNA.A_protein = par_DNA.alpha_A * par_DNA.D_damage / (par_DNA.D_damage + par_DNA.K_D) / par_DNA.mu_A;
    
    TSM.Chem_death = par_DNA.d_0 * par_DNA.A_protein / (par_DNA.A_protein + par_DNA.K_A_protein);
    TSM.growth = par_DNA.lambda_0 * par_DNA.D_damage / (par_DNA.D_damage + par_DNA.K_D);
    % PD-1/PD-L1 model
    TSM.rho8 = CY.p_1 * CY.p_L * CY.p_8 / CY.K_P;
    TSM.rhoL = CY.p_1 * CY.p_L;
    TSM.rhoh = CY.p_1 * CY.p_L * CY.p_h / CY.K_P;
    TSM.p_2 = CY.p_2;

    %% Dimension-less Parameter Transformation (Final Set for Simulation)
    Dim_TSM = struct();
    
    % Tumor growth dynamics (C)
    Dim_TSM.rC = TSM.r_C;
    Dim_TSM.growth =  TSM.growth;
    Dim_TSM.dC = TSM.d_C;
    Dim_TSM.etah = TSM.eta_h * TSM.T_40;
    Dim_TSM.eta8 = TSM.eta_8 * TSM.T_80;
    
    % Dendritic Cell dynamics (D)
    Dim_TSM.lambdaDC = TSM.lambda_DC;
    Dim_TSM.KDc = TSM.K_DC / TSM.C_max;
    Dim_TSM.dD = TSM.d_D;
    
    % CD8+ T Cell dynamics (T8)
    Dim_TSM.lambdaT8 = TSM.lambda_T8;
    Dim_TSM.lambdaT8h = TSM.T8h;
    Dim_TSM.alpha1 = TSM.alpha1 * TSM.M0;
    Dim_TSM.alphaD = TSM.alphaD * TSM.D_0;
    Dim_TSM.alpha2 = TSM.alpha2 * TSM.M0;
    Dim_TSM.alphaC = TSM.alphaC * TSM.C_max;
    Dim_TSM.alphaR = TSM.alphaR * TSM.T_40;
    Dim_TSM.Kh = TSM.K_h / TSM.T_40;
    Dim_TSM.d8 = TSM.d_8;
    
    % Helper T Cell dynamics (Th)
    Dim_TSM.lambdaTh = TSM.lambda_Th;
    Dim_TSM.lambdaThh = TSM.Thh;
    Dim_TSM.beta1 = TSM.beta1 * TSM.M0;
    Dim_TSM.betaD = TSM.betaD * TSM.D_0;
    Dim_TSM.beta2 = TSM.beta2 * TSM.M0;
    Dim_TSM.betaC = TSM.betaC * TSM.C_max;
    Dim_TSM.betaR = TSM.betaR * TSM.T_40;
    Dim_TSM.dh = TSM.d_h;
    
    % Regulatory T Cell dynamics (Tr)
    Dim_TSM.lambdaTr = TSM.lambda_Tr;
    Dim_TSM.gammaC = TSM.gammaC * TSM.C_max;
    Dim_TSM.gammaR = TSM.gammaR * TSM.T_40;
    Dim_TSM.gamma2 = TSM.gamma2 * TSM.M0;
    Dim_TSM.gammah = TSM.gammah * TSM.T_40;
    Dim_TSM.gamma8 = TSM.gamma8 * TSM.T_80;
    Dim_TSM.dr = TSM.d_r;
    
    % Macrophage dynamics (M1 & M2)
    Dim_TSM.lambdaM1 = TSM.lambda_M1;
    Dim_TSM.Kmc = TSM.K_MC / TSM.C_max;
    Dim_TSM.sigma1 = TSM.sigma1;
    Dim_TSM.sigma2 = TSM.sigma2;
    Dim_TSM.delta1 = TSM.delta1 * TSM.M0;
    Dim_TSM.delta2 = TSM.delta2 * TSM.M0;
    Dim_TSM.deltah = TSM.deltah * TSM.T_40;
    Dim_TSM.delta8 = TSM.delta8 * TSM.T_80;
    Dim_TSM.deltaR = TSM.deltaR * TSM.T_40;
    Dim_TSM.deltaC = TSM.deltaC * TSM.C_max;
    Dim_TSM.d1 = TSM.d_M1;
    Dim_TSM.d2 = TSM.d_M2;

    % Therapy specific constants
    Dim_TSM.a = TSM.a; Dim_TSM.b = TSM.b;
    Dim_TSM.R_D = TSM.R_D; Dim_TSM.R_8 = TSM.R_8;
    Dim_TSM.R_h = TSM.R_h; Dim_TSM.R_r = TSM.R_r;
    Dim_TSM.R_1 = TSM.R_1; Dim_TSM.R_2 = TSM.R_2;
    Dim_TSM.Dose = TSM.Dose;
    Dim_TSM.Chem_death = TSM.Chem_death;
    Dim_TSM.k_D = TSM.k_D; Dim_TSM.k_8 = TSM.k_8;
    Dim_TSM.k_h = TSM.k_h; Dim_TSM.k_r = TSM.k_r;
    Dim_TSM.k_1 = TSM.k_1; Dim_TSM.k_2 = TSM.k_2;
    
    % Immunotherapy scaling
    Dim_TSM.p_2 = TSM.p_2;
    Dim_TSM.A = TSM.A;
    Dim_TSM.rho8 = TSM.rho8 * TSM.T_40 * TSM.T_80;
    Dim_TSM.ep8 = TSM.ep8 * TSM.T_80 / TSM.T_40;
    Dim_TSM.epC = TSM.epC * TSM.C_max / TSM.T_40;
    Dim_TSM.rhoL = TSM.rhoL * TSM.T_40;
    Dim_TSM.rhoh = TSM.rhoh * TSM.T_40 * TSM.T_40;

    %% Re-package updated physical params
    new_params = params;
    new_params.par_DNA = par_DNA;
    new_params.CY = CY;
    new_params.par_cd4 = par_cd4;
    new_params.par_Th = par_Th;
    new_params.par_Tr = par_Tr;
    new_params.par_Mac = par_Mac;
    new_params.par_M1 = par_M1;
    new_params.par_M2 = par_M2;
    new_params.TSM = TSM;
    % Build one shared concentration/effect profile for simulation and scores.
    % Therapy columns: RT Gy, oral treatment day, platinum infusion, antibody dose.
    % Each injection entry is a dose event, not a persistent treatment flag.
    if nargin >= 3
        pk = params.PK;
        grid = linspace(0,1,round(1/pk.step_days)+1)';
        kF = log(2)/pk.F_half_days;
        kP = log(2)/pk.P_half_days;
        kA = log(2)/pk.A_half_days;
        % Solve Tmax = log(ka/kF)/(ka-kF), excluding the removable root at ka=kF.
        ka = fzero(@(x) log(x/kF)/(x-kF)-pk.F_formation_days, [kF*1.0001, kF*20]);
        qF = pk.q_F;
        qP = pk.q_P;
        qA = pk.q_A;
        profile = zeros(numel(grid),7,size(therapy,1));
        G0=0; F0=0; P0=0; A0=0;
        for day_index=1:size(therapy,1)
            F = F0*exp(-kF*grid) + ka*G0/(ka-kF)*(exp(-kF*grid)-exp(-ka*grid));
            Gend = G0*exp(-ka);
            if therapy(day_index,2)>0
                for oral_time=[0,0.5]
                    u=max(grid-oral_time,0);
                    F=F+therapy(day_index,2)*qF*ka/(ka-kF)*(exp(-kF*u)-exp(-ka*u));
                    Gend=Gend+therapy(day_index,2)*qF*exp(-ka*(1-oral_time));
                end
            end
            u=grid;
            P=P0*exp(-kP*u)+therapy(day_index,3)*qP/(pk.P_infusion_days*kP) ...
                .* (1-exp(-kP*min(u,pk.P_infusion_days))) .*exp(-kP*max(u-pk.P_infusion_days,0));
            A=(A0+therapy(day_index,4)*qA)*exp(-kA*u);
            damage=(par_DNA.lambda_F*F+par_DNA.lambda_P*P).*(1+par_DNA.rho_F*F)/par_DNA.k_b;
            growth_fraction=damage./(damage+par_DNA.K_D);
            protein=par_DNA.alpha_A*growth_fraction/par_DNA.mu_A;
            death_fraction=protein./(protein+par_DNA.K_A_protein);
            profile(:,:,day_index)=[F,P,A,par_DNA.d_0*death_fraction, ...
                par_DNA.lambda_0*growth_fraction,death_fraction,growth_fraction];
            G0=Gend; F0=F(end); P0=P(end); A0=A(end);
        end
        new_params.drug_profile = profile;
        new_params.drug_grid = grid;
    end

end