function exp_data = initialize_experiment_conditions(Therapy_matrix, params)
%% FUNCTION NAME: initialize_experiment_conditions
%
% DESCRIPTION: 
%     Sets initial state variables and simulation timeframes for in silico trials.
%     Initial tumor cell counts are randomized to represent baseline patient 
%     heterogeneity, while immune populations start at a standard baseline.
%
% INPUTS:
%     Therapy_matrix - Matrix containing scheduling
%     params         - Global parameter structure (from init_global_params)
%
% OUTPUTS:
%     exp_data - Structure containing:
%              .Y0: Initial population vector
%              .t: Time vector for simulation steps
%              .Therapy: Extracted therapy schedule
%              .Volume: Pre-allocated volume matrix
%
%---------------------------------------------------------

    % Set fixed seed for reproducibility of random patient initialization
    rng(123);

    % Pre-allocate empty matrix for volume tracking
    Volume_matrix = [];

    % Randomize initial cancer cell count (30%-60% of carrying capacity C_max)
    % This represents the primary source of initial tumor burden heterogeneity
    C0 = params.TSM.C_max * (0.3 + 0.3 * rand); 

    % Baseline counts for immune cell populations (set to 1e6 per cell type)
    D0  = 1e6; % Dendritic Cells
    Th0 = 1e6; % Helper T cells
    T80 = 1e6; % CD8+ T cells
    Tr0 = 1e6; % Regulatory T cells
    M10 = 1e6; % M1-like Macrophages
    M20 = 1e6; % M2-like Macrophages

    % Construct the experimental data structure
    exp_data = struct(...
        'Y0', [C0, D0, Th0, T80, Tr0, M10, M20],...
        't', 0:1:max(Therapy_matrix(:,1)), ...
        'Therapy',Therapy_matrix(:,2:5), ...
        'Volume', Volume_matrix); 
end