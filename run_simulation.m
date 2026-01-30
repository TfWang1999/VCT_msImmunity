function  [t_sol, Vtumor_list,Y_sol,Y_sol_number,parameters] = run_simulation(experimental_data, params)
%% FUNCTION NAME: run_multiscale_simulation
%
% DESCRIPTION: 
%     Executes the main multiscale simulation loop. It integrates population 
%     dynamics with dynamic volume feedback, where macroscopic tumor size 
%     re-scales molecular concentrations at each daily time step.
%
% INPUTS:
%     experimental_data - Struct containing Y0, t, and Therapy schedule
%     params            - Global parameter structure
%
% OUTPUTS:
%     t_sol        - Time points of simulation results
%     Vtumor_list  - Calculated tumor volume at each time point (cm^3)
%     Y_sol        - Dimensionless state variable trajectories
%     Y_sol_number - Absolute cell count trajectories
%     parameters   - Final state of parameters used in the last step
%
%---------------------------------------------------------

    CY = params.CY; 
    TSM = params.TSM; 
    
    % M_vector: Characteristic volume constants for cell types
    M_vector = [2572.44078451444, 175.015709799539, 4849.04826081585, 179.594380030216]; 
    
    Y0 = experimental_data.Y0;
    t_total = experimental_data.t;
    n_steps = length(t_total) - 1;
    
    % Convert initial cell counts to macroscopic tumor volume (cm^3)
    % Scaling constant 1e-12/0.37 accounts for cell volume and intracellular fraction
    Vtumor_0 = sum(bsxfun(@times, [Y0(1),Y0(3)+Y0(4)+Y0(5),Y0(6)+Y0(7),Y0(2)],M_vector)) * 1e-12 / 0.37;
    
    % Handle control group (no treatment) defaults
    if ~isfield(experimental_data, 'Therapy') || isempty(experimental_data.Therapy) 
        params.TSM.A = 0;
        params.par_DNA.F_ave = 0;
        params.par_DNA.P_ave = 0;
    end
    
    scale_vector = [TSM.C_max, TSM.D_0, TSM.T_80, TSM.T_40, TSM.T_40, TSM.M0, TSM.M0];
    
    % Pre-allocate storage for trajectories
    Y_sol = zeros(n_steps+1, 7);
    Y_sol(1,:) = bsxfun(@rdivide, Y0, scale_vector);
    t_sol = zeros(n_steps+1, 1);
    t_sol(1) = t_total(1);
    Vtumor_list = zeros(n_steps+1, 1);
    Vtumor_list(1) = Vtumor_0;
    Y_sol_number = zeros(n_steps+1, 7);
    Y_sol_number(1,:) = Y0;
    
    opts_15s = odeset('RelTol', 1e-4, 'AbsTol', 1e-7, 'MaxStep', 0.1);
    opts_15s = odeset('RelTol', 1e-4, 'AbsTol', 1e-7, 'MaxStep', 0.1);
    set = 0;% 0: ode45, 1: switch to ode15s for stiff systems

    %% Main Simulation Loop (Step-wise Integration)
    for i = 1:n_steps
        t_span = [t_total(i), t_total(i+1)];
        
        % Update therapy dose and concentrations based on schedule
        if isfield(experimental_data, 'Therapy') && ~isempty(experimental_data.Therapy) 
            % RT dose /Gy
            params.TSM.Dose = experimental_data.Therapy(i,1);
            
            % Set Anti-PD-1 state
            if experimental_data.Therapy(i,4) == 1
                params.TSM.A = params.CY.A;
            else
                params.TSM.A = 0;
            end
            
            % Set Chemotherapy component F (e.g., Capecitabine)
            if experimental_data.Therapy(i,2) == 1
                params.par_DNA.F_ave = CY.F_ave;
            else
                params.par_DNA.F_ave = 0;
            end
            
            % Set Chemotherapy component P (e.g., Oxaliplatin)
            if experimental_data.Therapy(i,3) == 1
                params.par_DNA.P_ave = CY.P_ave;
            else
                params.par_DNA.P_ave = 0;
            end
        else
           params.par_DNA.F_ave = 0;
           params.par_DNA.P_ave = 0;
           params.TSM.Dose = 0;
        end

        % CRITICAL STEP: Dynamic re-scaling of parameters based on current volume
        [Dim_TSM, ~] = calculate_multiscale_scaling(Vtumor_list(i), params);
        
        % Adaptive numerical integration to handle potential system stiffness
        if set == 0
        [t_temp, Y_temp, set] = adaptive_ode_solver(@(t,Y) Dimless_tumor(t, Y, Dim_TSM), t_span, Y_sol(i,:), 5);
        elseif solver_mode == 1
        [t_temp, Y_temp, set] = ode15s(@(t,Y) Dimless_tumor(t, Y, Dim_TSM), t_span, Y_sol(i,:), opts_15s);
        end

        % Store dimensionless results
        Y_sol(i+1,:) = Y_temp(end,:);
        t_sol(i+1) = t_temp(end);
       
        % Calculate absolute cell numbers and update macroscopic tumor volume
        current_Y = bsxfun(@times, Y_temp(end,:), scale_vector);
        Y_sol_number(i+1,:) = current_Y;
        
        % Volumetric update for the next time step
        Vtumor_new = sum(bsxfun(@times, [current_Y(1),current_Y(3)+current_Y(4)+current_Y(5),...
            current_Y(6)+current_Y(7),current_Y(2)],M_vector)) * 1e-12 / 0.37;
        
        Vtumor_list(i+1) = Vtumor_new;
    end
    
    params.Dim_TSM = Dim_TSM;
    parameters = params;
end

function  [t, Y, set] = adaptive_ode_solver(ode_fun, t_span, Y0, time_limit)
%% INTERNAL HELPER: adaptive_ode_solver
% Attempts to solve using ode45 but switches to ode15s if computation exceeds limit.
    opts_45 = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
    tic;
    [t, Y] = ode45(ode_fun, t_span, Y0, opts_45);
    elapsed_time = toc;
    set = 0;
    
    if elapsed_time > time_limit
        fprintf('Computational overhead detected (%.2fs). Switching to stiff solver (ode15s)...\n', elapsed_time);
        set = 1;
        opts_15s = odeset('RelTol', 1e-4, 'AbsTol', 1e-7, 'MaxStep', 0.1);
        [t, Y] = ode15s(ode_fun, t_span, Y0, opts_15s);
    end
end