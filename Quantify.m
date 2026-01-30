function [Score_matrix,Par_matrix] = Quantify(params,Therapy_martix)


[Dim_TSM,params] = calculate_multiscale_scaling(1,params);
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

   
% Approximate steady-state values
    C_m = (rC-dC)/rC;
    T8_m = (lambdaT8*0.5*0.5 * (1 + lambdaT8h*0.5) )   / d8  *rhoL/rho8;
    Th_m = (lambdaTh*0.5*0.5 * (1 + lambdaThh*0.5) )  / dh  *rhoL/rhoh; 
    M1_m = lambdaM1  * C_m/(Kmc + C_m) * (sigma2*0.5 + d2)/(sigma1*0.5*d2+sigma2*0.5*d1+d1*d2);
    M2_m = lambdaM1  * C_m/(Kmc + C_m) * sigma1*0.5/(sigma1*0.5*d2+sigma2*0.5*d1+d1*d2);
    Tr_m = lambdaTr*0.5*0.5 / dr;
    DC_m = lambdaDC*C_m/(KDc + C_m) / dD;

% Time-Weighted Therapeutic Parameters
    % growth_ave = 0;
    % death_RT = 0;
    % death_Chem = 0;

    par_DNA = params.par_DNA;
    CY = params.CY;
    TSM = params.TSM;

   
    % CY=Par_all.CY; par_DNA = Par_all.par_DNA;
    par_DNA.D_damage = (par_DNA.lambda_F*CY.F_ave + par_DNA.lambda_P*CY.P_ave) * (1+par_DNA.rho_F*CY.F_ave) / par_DNA.k_b;
    par_DNA.A_protein = par_DNA.alpha_A * par_DNA.D_damage/(par_DNA.D_damage+par_DNA.K_D) / par_DNA.mu_A;
    Chem_death_1 = par_DNA.d_0 * par_DNA.A_protein / (par_DNA.A_protein + par_DNA.K_A_protein);
    growth_1 =par_DNA.lambda_0 * par_DNA.D_damage/(par_DNA.D_damage+par_DNA.K_D);

    CY.F_ave = 0;
    par_DNA.D_damage = (par_DNA.lambda_F*CY.F_ave + par_DNA.lambda_P*CY.P_ave) * (1+par_DNA.rho_F*CY.F_ave) / par_DNA.k_b;
    par_DNA.A_protein = par_DNA.alpha_A * par_DNA.D_damage/(par_DNA.D_damage+par_DNA.K_D) / par_DNA.mu_A;
    Chem_death_2 = par_DNA.d_0 * par_DNA.A_protein / (par_DNA.A_protein + par_DNA.K_A_protein);
    growth_2 =par_DNA.lambda_0 * par_DNA.D_damage/(par_DNA.D_damage+par_DNA.K_D);
    
    death_RT =0;
    if ~isempty(Therapy_martix)
        if sum(Therapy_martix(:,2)) > 0
            dose = Therapy_martix(:,2);
            death_RT = mean( (1 - exp(-a .* dose - b .* dose.^2)));
        end
    death_Chem =0;
        if sum(Therapy_martix(:,4)) > 0
            Chem_death_ave = Chem_death_1 * 2/3 + Chem_death_2 * 1/3;
            growth_ave = growth_1 * 2/3 + growth_2 * 1/3;
            death_Chem = Chem_death_ave;
        end
    end
    
    
  
    %% calaulation of msImmunity
    % PD1-PD-L1 suppression
    sup_8 = M_sup( (rho8*T8_m*( Th_m+ep8*T8_m+epC*C_m )) / (1+rhoL*(Th_m+ep8*T8_m+epC*C_m)+p_2*A) );
    sup_h = M_sup( rhoh*Th_m*( Th_m+ep8*T8_m+epC*C_m )/(1+rhoL*(Th_m+ep8*T8_m+epC*C_m)+p_2*A ) );
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
    
    
    msImmunity =  F_m*( (T8_score + Th_score)/(rC - growth_ave -dC))...
        + (death_Chem/(rC - growth_ave-dC)...
        + death_RT/(rC - growth_ave-dC)) ;




     %% calaulation of tsImmunity
    T8_score_tsm = eta8* (lambdaT8 / d8) * (1 + lambdaT8h / d8 );
    Th_score_tsm = etah* (lambdaTh / dh) * (1 + lambdaThh / dh );
    tsImmunity =  ( (T8_score_tsm + Th_score_tsm + par_DNA.d_0 *1/2)/(TSM.r_C-par_DNA.lambda_0*1/2-dC) ) + death_RT/(TSM.r_C-par_DNA.lambda_0*1/2-dC) ;
    
    %% calaulation of molsImmunity
    par_DNA.D_damage = (par_DNA.lambda_F*CY.F_ave + par_DNA.lambda_P*CY.P_ave) * (1+par_DNA.rho_F*CY.F_ave) / par_DNA.k_b;
    par_DNA.A_protein = par_DNA.alpha_A * par_DNA.D_damage/(par_DNA.D_damage+par_DNA.K_D) / par_DNA.mu_A;
    DNA_damage =  par_DNA.A_protein / (par_DNA.A_protein + par_DNA.K_A_protein);
    A_pro = par_DNA.D_damage/(par_DNA.D_damage+par_DNA.K_D);

    %% calaulation of molsImmunity
    molsImmunity =( 1+(gamma8+gammah)/(1+gamma8+gammah)+(delta1+deltah+delta8)/(1+delta1+deltah+delta8)...
        + 1/(1+gamma2+gammaR+gammaC) +1/(1+deltaR+deltaC+delta2) ) *...
        par_DNA.A_protein / (par_DNA.A_protein + par_DNA.K_A_protein) * par_DNA.D_damage/(par_DNA.D_damage+par_DNA.K_D)...
        *(rhoL/rhoh+rhoL/rho8);


  
  

    Score_matrix = [msImmunity,tsImmunity,molsImmunity];


    Par_matrix = cell2mat(struct2cell(Dim_TSM))';

    % Immune_score = rC/(rC -dC);


end

function M_pro = M_pro(x)
    M_pro = x/(1+x);
end
function M_sup = M_sup(x)
    M_sup = 1/(1+x);
end