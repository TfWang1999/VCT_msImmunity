function dYdt = Dimless_tumor(~, Y, params)
%% FUNCTION NAME: Dimless_tumor
%
% DESCRIPTION: 
%     Calculates the derivatives for the dimensionless multiscale tumor-immune 
%     dynamics model, integrating radiotherapy, chemotherapy, and immunotherapy.
%
% INPUTS:
%     ~      - Time (unused in this autonomous system)
%     Y      - Vector of state variables
%     params - Struct containing calibrated parameters extracted from init_global_params
%
% OUTPUTS:
%     dYdt   - Vector of rates of change for each cell population
%

%---------------------------------------------------------

    % Extract State Variables
    C = Y(1);  % Cancer cells
    D = Y(2);  % Dendritic Cells (DCs)
    T8 = Y(3); % CD8+ T cells
    Th = Y(4); % Helper T cells
    Tr = Y(5); % Regulatory T cells
    M1 = Y(6); % M1-like Macrophages
    M2 = Y(7); % M2-like Macrophages
    
    % Extract Parameters from Struct (Standard Tissue Parameters)
    rC = params.rC - params.growth;
    dC = params.dC;
    etah = params.etah;
    eta8 = params.eta8;
    lambdaDC = params.lambdaDC;
    KDc = params.KDc;
    dD = params.dD;
    lambdaT8 = params.lambdaT8;
    lambdaT8h = params.lambdaT8h;          
    ep8 = params.ep8;
    epC = params.epC;                     
    alpha1 = params.alpha1;                
    alphaD = params.alphaD;
    alpha2 = params.alpha2;                
    alphaC = params.alphaC;
    alphaR = params.alphaR;
    Kh = params.Kh;
    d8 = params.d8;                      
    lambdaTh = params.lambdaTh;
    lambdaThh = params.lambdaThh;          
    beta1 = params.beta1;                  
    betaD = params.betaD;
    beta2 = params.beta2;                  
    betaC = params.betaC;
    betaR = params.betaR;
    dh = params.dh;                        
    lambdaTr = params.lambdaTr;
    gammaC = params.gammaC;
    gammaR = params.gammaR;                
    gamma2 = params.gamma2;               
    gammah = params.gammah;
    gamma8 = params.gamma8;
    dr = params.dr;                       
    lambdaM1 = params.lambdaM1;
    Kmc = params.Kmc;
    sigma1 = params.sigma1;
    sigma2 = params.sigma2;
    delta1 = params.delta1;                
    delta2 = params.delta2;                
    deltah = params.deltah;                
    delta8 = params.delta8;
    deltaR = params.deltaR;                
    deltaC = params.deltaC;
    d1 = params.d1;
    d2 = params.d2;

    % Therapy Modulation Parameters
    R_1 = params.R_1; R_2 = params.R_2;
    R_8 = params.R_8; R_h = params.R_h;
    R_D = params.R_D; R_r = params.R_r;
    a = params.a;     b = params.b;
    Dose = params.Dose;
    
    % Chemotherapy effect
    Chem_death = params.Chem_death;
    k_D = params.k_D; k_8 = params.k_8;
    k_h = params.k_h; k_r = params.k_r;
    k_1 = params.k_1; k_2 = params.k_2;
    
    % Immunotherapy effect (PD-1 blockade)
    A = params.A;         % Anti-PD-1 concentration
    rhoL = params.rhoL;   p_2 = params.p_2;
    rhoh = params.rhoh;   rho8 = params.rho8;

    %% --- Calculation of Derivatives ---

    % dC/dt: Cancer Cells (Logistic growth - immune killing - therapy induced death)
    dCdT = rC * C * (1 - C) - (etah * Th * C + eta8 * T8 * C) - dC * C...
        - (1 - exp(-a * Dose - b * Dose^2)) * C - Chem_death * C;

    % dD/dt: Dendritic Cells (Activation - death - therapy side effects)
    dDdT = lambdaDC * C / (KDc + C) - dD * D...
        - R_D * (1 - exp(-a * Dose - b * Dose^2)) * D - k_D * Chem_death * C;

    % dT8/dt: CD8+ T Cells (Activation modulated by M1/D and inhibited by PD-1/PD-L1)
    dT8dT = (lambdaT8 * (alpha1 * M1 + alphaD * D) / (1 + alpha1 * M1 + alphaD * D)...
        / (1 + alpha2 * M2 + alphaC * C + alphaR * Tr)...
        + lambdaT8h * (Th / (Th + Kh)) * T8) *...
        (1 / (1 + (rho8 * T8 * (Th + ep8 * T8 + epC * C) / (1 + rhoL * (Th + ep8 * T8 + epC * C) + p_2 * A))))...
        - d8 * T8 - R_8 * (1 - exp(-a * Dose - b * Dose^2)) * T8 - k_8 * Chem_death * T8;

    % dTh/dt: Helper T Cells (Activation modulated by microenvironment and PD-1 blockade)
    dThdT = (lambdaTh * (beta1 * M1 + betaD * D) / (1 + beta1 * M1 + betaD * D)...
        / (1 + beta2 * M2 + betaC * C + betaR * Tr)...
        + lambdaThh * (Th / (Th + Kh)) * Th) *...
        (1 / (1 + (rhoh * T8 * (Th + ep8 * T8 + epC * C) / (1 + rhoL * (Th + ep8 * T8 + epC * C) + p_2 * A))))...
        - dh * Th - R_h * (1 - exp(-a * Dose - b * Dose^2)) * Th - k_h * Chem_death * Th;

    % dTr/dt: Regulatory T Cells (Recruitment from M2/C/Tr and inhibited by Th/T8)
    dTrdT = lambdaTr * (gamma2 * M2 + gammaC * C + gammaR * Tr) / (1 + gamma2 * M2 + gammaC * C + gammaR * Tr)...
        * (1 / (1 + gammah * Th + gamma8 * T8)) - dr * Tr...
        - R_r * (1 - exp(-a * Dose - b * Dose^2)) * Tr - k_r * Chem_death * Tr;

    % dM1/dt: M1 Macrophages (Recruitment - polarization switch - death)
    dM1dT = lambdaM1 * C / (Kmc + C)...
        + sigma2 * M2 * (deltah * Th + delta8 * T8 + delta1 * M1) / (1 + deltah * Th + delta8 * T8 + delta1 * M1)...
        - sigma1 * M1 * (deltaC * C + deltaR * Tr + delta2 * M2) / (1 + deltaC * C + deltaR * Tr + delta2 * M2)...
        - d1 * M1 - R_1 * (1 - exp(-a * Dose - b * Dose^2)) * M1 - k_1 * Chem_death * M1;

    % dM2/dt: M2 Macrophages (Recruitment - polarization switch - death)
    dM2dT = sigma1 * M1 * (deltaC * C + deltaR * Tr + delta2 * M2) / (1 + deltaC * C + deltaR * Tr + delta2 * M2)...
        - sigma2 * M2 * (deltah * Th + delta8 * T8 + delta1 * M1) / (1 + deltah * Th + delta8 * T8 + delta1 * M1)...
        - d2 * M2 - R_2 * (1 - exp(-a * Dose - b * Dose^2)) * M2 - k_2 * Chem_death * M2;

    % Return derivative vector
    dYdt = [dCdT; dDdT; dT8dT; dThdT; dTrdT; dM1dT; dM2dT];
end