function [Score_matrix,Par_matrix,diagnostics] = Quantify(params,Therapy_martix,abundance_approx,allow_undefined_scores)

% Quantify - Calculates immune score metrics (msImmunity, tsImmunity, molsImmunity) based on a multi-scale mathematical model
%
%
% Input:
% params: A structure containing global parameters, DNA damage parameters, and cytokine parameters
% Therapy_martix: A treatment protocol matrix recording the dose and duration of radiotherapy (RT) and chemotherapy (Chem)
% abundance_approx: Optional cell-abundance approximation(s), default 0.5.
% Only the approximate steady-state cell abundances use this value. The
% integration weights and tissue-only molecular approximation remain fixed.
% A vector returns one score row per value, sharing the same PK calculation.
% Scores are evaluated directly from the formulas without a positivity gate.
% allow_undefined_scores: Optional logical, default false. Sensitivity scans
% may request NaN for nonpositive or nonfinite score denominators.
%
% Output:
% Score_matrix: A three-level score matrix [msImmunity, tsImmunity, molsImmunity]
% Par_matrix: A dimensionless parameter matrix after scale adjustment

if nargin < 3; abundance_approx = 0.5; end
if nargin < 4; allow_undefined_scores = false; end
validateattributes(allow_undefined_scores,{'logical'},{'scalar'},mfilename,'allow_undefined_scores');
validateattributes(abundance_approx,{'numeric'}, ...
    {'real','finite','vector','nonempty','>',0,'<',1},mfilename,'abundance_approx');

if isempty(Therapy_martix)
    therapy=zeros(1,4);
else
    therapy=Therapy_martix(:,2:5);
end
[Dim_TSM,params] = calculate_multiscale_scaling(1,params,therapy);
   % Dim_TSM = Par_all.Dim_TSM;
    params.TSM.r_C;
    rC = Dim_TSM.rC;
    dC =  Dim_TSM.dC;
    etah =  Dim_TSM.etah;
    eta8 =  Dim_TSM.eta8;
    lambdaDC =  Dim_TSM.lambdaDC;
    KDc =  Dim_TSM.KDc;
    dD =  Dim_TSM.dD;
    lambdaT8 =  Dim_TSM.lambdaT8;
    lambdaT8h =  Dim_TSM.lambdaT8h;          
    ep8 =  Dim_TSM.ep8;
    epC =  Dim_TSM.epC;                     
    alpha1 =  Dim_TSM.alpha1;                
    alphaD =  Dim_TSM.alphaD;
    alpha2 =  Dim_TSM.alpha2;                
    alphaC =  Dim_TSM.alphaC;
    alphaR =  Dim_TSM.alphaR;
    Kh =  Dim_TSM.Kh;
    d8 =  Dim_TSM.d8;                      
    lambdaTh =  Dim_TSM.lambdaTh;
    lambdaThh =  Dim_TSM.lambdaThh;          
    beta1 =  Dim_TSM.beta1;                  
    betaD =  Dim_TSM.betaD;
    beta2 =  Dim_TSM.beta2;                  
    betaC =  Dim_TSM.betaC;
    betaR =  Dim_TSM.betaR;
    dh =  Dim_TSM.dh;                        
    lambdaTr =  Dim_TSM.lambdaTr;
    gammaC =  Dim_TSM.gammaC;
    gammaR =  Dim_TSM.gammaR;                
    gamma2 =  Dim_TSM.gamma2;               
    gammah =  Dim_TSM.gammah;
    gamma8 =  Dim_TSM.gamma8;
    dr =  Dim_TSM.dr;                       
    lambdaM1 =  Dim_TSM.lambdaM1;
    Kmc =  Dim_TSM.Kmc;
    sigma1 =  Dim_TSM.sigma1;
    sigma2 =  Dim_TSM.sigma2;
    delta1 =  Dim_TSM.delta1;                
    delta2 =  Dim_TSM.delta2;                
    deltah =  Dim_TSM.deltah;                
    delta8 =  Dim_TSM.delta8;
    deltaR =  Dim_TSM.deltaR;                
    deltaC =  Dim_TSM.deltaC;
    d1 =  Dim_TSM.d1;
    d2 =  Dim_TSM.d2;
    % Therapy
  
    a =  Dim_TSM.a;
    b =  Dim_TSM.b;
   
    % chemotherapy
    Chem_death =  Dim_TSM.Chem_death;
    % k_C =  Dim_TSM.k_C;
    k_D =  Dim_TSM.k_D;
    k_8 =  Dim_TSM.k_8;
    k_h =  Dim_TSM.k_h;
    k_r =  Dim_TSM.k_r;
    k_1 =  Dim_TSM.k_1;
    k_2 =  Dim_TSM.k_2;
    % Immunotherapy
    A =      Dim_TSM.A; % Anti-PD-1
    rhoL =   Dim_TSM.rhoL;
    p_2 =    Dim_TSM.p_2;
    rhoh =   Dim_TSM.rhoh;
    rho8 =   Dim_TSM.rho8;

   
% Approximate steady-state values. Reuse the unchanged PK profile for a scan.
abundance_values = abundance_approx(:);
Score_matrix = zeros(numel(abundance_values),3);
for abundance_index = 1:numel(abundance_values)
    abundance_approx = abundance_values(abundance_index);
    C_m = (rC-dC)/rC;
    T8_m = (lambdaT8*abundance_approx*abundance_approx * (1 + lambdaT8h*abundance_approx) )   / d8  *rhoL/rho8;
    Th_m = (lambdaTh*abundance_approx*abundance_approx * (1 + lambdaThh*abundance_approx) )  / dh  *rhoL/rhoh;
    M1_m = lambdaM1  * C_m/(Kmc + C_m) * (sigma2*abundance_approx + d2)/(sigma1*abundance_approx*d2+sigma2*abundance_approx*d1+d1*d2);
    M2_m = lambdaM1  * C_m/(Kmc + C_m) * sigma1*abundance_approx/(sigma1*abundance_approx*d2+sigma2*abundance_approx*d1+d1*d2);
    Tr_m = lambdaTr*abundance_approx*abundance_approx / dr;
    DC_m = lambdaDC*C_m/(KDc + C_m) / dD;

% Average nonlinear effects over the same treatment horizon as the ODE.
    par_DNA=params.par_DNA;
    TSM=params.TSM;
    profile=params.drug_profile;
    n_days=size(profile,3);
    weights=ones(size(profile,1),1); weights([1,end])=0.5;
    weights=weights*params.drug_grid(2)/n_days;
    death_Chem=sum(sum(profile(:,4,:).*weights,1),3);
    growth_ave=sum(sum(profile(:,5,:).*weights,1),3);
    death_RT=mean(1-exp(-a*therapy(:,1)-b*therapy(:,1).^2));
    ligand=Th_m+ep8*T8_m+epC*C_m;
    blockade_denominator=1+rhoL*ligand+p_2*profile(:,3,:);
    sup_8=sum(sum((1./(1+rho8*T8_m*ligand./blockade_denominator)).*weights,1),3);
    sup_h=sum(sum((1./(1+rhoh*Th_m*ligand./blockade_denominator)).*weights,1),3);
    molecular_chem=sum(sum(profile(:,6,:).*profile(:,7,:).*weights,1),3);
    % The tissue-only comparator retains the fixed 0.5 molecular approximation.
    % Apply it on chemotherapy treatment days, without adding molecular PK terms.
    chemo_fraction=mean(therapy(:,2)>0 | therapy(:,3)>0);
    net_growth=rC-growth_ave-dC;
    tissue_net_growth=TSM.r_C-par_DNA.lambda_0*0.5*chemo_fraction-dC;
    diagnostics = struct('ms_denominator',net_growth,'ts_denominator',tissue_net_growth);

    % define of T8_score and Th_score
    T8_score =eta8* (lambdaT8 * M_pro(alpha1*M1_m + alphaD * DC_m) * M_sup(alpha2*M2_m+alphaR*Tr_m+alphaC*C_m) / d8 ...
    * (1 + lambdaT8h * Th_m/(Th_m+Kh)/ d8 ) )*  sup_8;

    Th_score =etah* (lambdaTh * M_pro(beta1*M1_m + betaD * DC_m) * M_sup(beta2*M2_m+betaR*Tr_m+betaC*C_m)/ dh...
    * (1 + lambdaThh * Th_m/(Th_m+Kh)/ dh ) )* sup_h;
    
    
    pro_1 = (gamma8*T8_m + gammah*Th_m);
    pro_2 = (delta1 *M1_m+ deltah *Th_m + delta8*T8_m);

    sup_1 =  (gamma2*M2_m + gammaR*Tr_m + gammaC*C_m);
    sup_2 = (deltaR*Tr_m +deltaC*C_m +delta2*M2_m);
    
    F_m = 1 + M_pro(pro_1) + M_pro(pro_2) + M_sup(sup_1) +M_sup(sup_2);
    
    
    msImmunity =  F_m*( (T8_score + Th_score)/net_growth)...
        + (death_Chem/net_growth...
        + death_RT/net_growth) ;




     %% calaulation of tsImmunity
    T8_score_tsm = eta8* (lambdaT8 / d8) * (1 + lambdaT8h / d8 );
    Th_score_tsm = etah* (lambdaTh / dh) * (1 + lambdaThh / dh );
    tsImmunity=(T8_score_tsm+Th_score_tsm+par_DNA.d_0*0.5*chemo_fraction+death_RT) ...
        /tissue_net_growth;

    %% Molecular-only score: average the product of the two drug-response fractions.
    molsImmunity =( 1+(gamma8+gammah)/(1+gamma8+gammah)+(delta1+deltah+delta8)/(1+delta1+deltah+delta8)...
        + 1/(1+gamma2+gammaR+gammaC) +1/(1+deltaR+deltaC+delta2) ) *...
        molecular_chem...
        *(rhoL/rhoh+rhoL/rho8);


  
  

    if allow_undefined_scores
        invalid_denominator = ~isfinite([net_growth,tissue_net_growth]) | [net_growth,tissue_net_growth]<=0;
        if invalid_denominator(1); msImmunity = NaN; end
        if invalid_denominator(2); tsImmunity = NaN; end
    end
    Score_matrix(abundance_index,:) = [msImmunity,tsImmunity,molsImmunity];
end


    Par_matrix = cell2mat(struct2cell(Dim_TSM))';

    % Immune_score = rC/(rC -dC);


end

function M_pro = M_pro(x)
    M_pro = x/(1+x);
end
function M_sup = M_sup(x)
    M_sup = 1/(1+x);
end
