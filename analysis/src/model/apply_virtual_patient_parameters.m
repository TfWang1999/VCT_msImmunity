function params = apply_virtual_patient_parameters(sample_row)
%APPLY_VIRTUAL_PATIENT_PARAMETERS Map all 76 sampled values to active fields.

if numel(sample_row) ~= 76
    error('apply_virtual_patient_parameters:Size', 'Expected exactly 76 sampled values.');
end
S = sample_row(:)';
params = init_global_params();
TSM = params.TSM; CY = params.CY; par_DNA = params.par_DNA;
par_cd4 = params.par_cd4; par_Th = params.par_Th; par_Tr = params.par_Tr;
par_Mac = params.par_Mac; par_M1 = params.par_M1; par_M2 = params.par_M2;

TSM.r_C=S(1); TSM.C_max=S(2); TSM.eta_h=S(3); TSM.eta_8=S(4);
TSM.lambda_DC=S(5); TSM.D_0=S(6); TSM.K_DC=S(7);
TSM.lambda_T8=S(8); TSM.T_80=S(9); TSM.T8h=S(10);
TSM.lambda_Th=S(11); TSM.T_40=S(12); TSM.Thh=S(13); TSM.lambda_Tr=S(14);
TSM.lambda_M1=S(15); TSM.M0=S(16); TSM.K_MC=S(17); TSM.sigma1=S(18); TSM.sigma2=S(19);
par_cd4.f_1=S(20); par_cd4.s_I=S(21); par_cd4.K1=S(22);
par_cd4.f_2=S(23); par_cd4.s_T=S(24); par_cd4.K2=S(25);
par_Th.a_IFNg_0=S(26); par_Th.a_IL2_0=S(27); par_Th.a_TNFa_0=S(28); par_Th.K_Tb=S(29);
par_Tr.a_TGFb_0=S(30); par_Tr.a_IL10_0=S(31); par_Tr.K_FOXP=S(32);
par_Mac.f_1=S(33); par_Mac.s_I=S(34); par_Mac.K1=S(35);
par_Mac.f_2=S(36); par_Mac.s_T=S(37); par_Mac.K2=S(38);
par_Mac.f_3=S(39); par_Mac.s_tgfbr=S(40); par_Mac.s_10=S(41); par_Mac.K3=S(42);
par_M1.b_IL12_0=S(43); par_M1.b_TNFa_0=S(44); par_M1.Km1=S(45);
par_M2.b_IL10_0=S(46); par_M2.b_TGFb_0=S(47); par_M2.Km2=S(48);
CY.p_1=S(49); CY.p_8=S(50); CY.p_h=S(51); CY.p_L=S(52); CY.K_P=S(53);
TSM.epC=S(54); TSM.a=S(55); TSM.b=S(55)/10;
par_DNA.d_0=S(56); par_DNA.lambda_0=S(57); par_DNA.lambda_P=S(58);
par_DNA.lambda_F=S(59); par_DNA.rho_F=S(60); par_DNA.k_b=S(61);
par_DNA.alpha_A=S(62); par_DNA.K_A_protein=S(63); par_DNA.K_D=S(64);
params.PK.q_A=S(65); CY.p_2=S(66); CY.K_IL12=S(67); CY.K_IL2=S(68);
CY.K_IL10=S(69); CY.K_IFNg=S(70); CY.K_TNFa=S(71); CY.K_TGFb=S(72);
CY.lambda_DC_IL12=S(73); CY.lambda_C_TGFb=S(74);
CY.lambda_C_IL10=S(75); CY.lambda_T8_IFNg=S(76);

params.par_DNA=par_DNA; params.CY=CY; params.par_cd4=par_cd4;
params.par_Th=par_Th; params.par_Tr=par_Tr; params.par_Mac=par_Mac;
params.par_M1=par_M1; params.par_M2=par_M2; params.TSM=TSM;
end

