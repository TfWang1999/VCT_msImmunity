function result = fit_crc_case(case_data, options)
%FIT_CRC_CASE Fit patient-specific CRC parameters by ordinary least squares.
% Ordinary volume least squares: minimizing SSE maximizes R2 for each patient.

arguments
    case_data (1,1) struct
    options.NumStarts (1,1) double {mustBeInteger,mustBePositive} = 3
    options.MaxIterations (1,1) double {mustBeInteger,mustBePositive} = 120
    options.Display (1,:) char = 'off'
    options.FitProfile (1,1) string = ""
    options.InitialValues (1,:) double = []
    options.InitialParameters (1,1) struct = struct()
    options.LowerBounds (1,:) double = []
    options.UpperBounds (1,:) double = []
end

baseline_params = apply_crc_manual_parameters(case_data);
if strlength(options.FitProfile) == 0
    options.FitProfile = case_data.fit_profile;
end
switch lower(options.FitProfile)
    case "fixed_structure_pd1"
        % Keep all non-fitted inputs identical to the supplied common model.
        if ~isempty(fieldnames(options.InitialParameters))
            baseline_params = options.InitialParameters;
        end
        parameter_names = {'TSM.r_C','par_DNA.d_0','par_DNA.lambda_0','CY.p_1'};
        lower_bound = [0.25,0.10,0.10,1];
        upper_bound = [0.60,0.35,0.40,10];
        if any(case_data.experimental_data.Therapy(:,4)>0)
            parameter_names{end+1} = 'CY.p_2';
            lower_bound(end+1) = 1e3;
            upper_bound(end+1) = 1e6;
        end
    case {"fixed_a4","fixed_a_pd1_5"}
        parameter_names = {'TSM.r_C','TSM.C_max','par_DNA.d_0','par_DNA.lambda_0'};
        manual_cmax = baseline_params.TSM.C_max;
        lower_bound = [0.25, 1.02 * case_data.experimental_data.Y0(1), 0.15, 0.20];
        upper_bound = [0.60, 2.0 * manual_cmax, 0.35, 0.40];
        if options.FitProfile == "fixed_a_pd1_5"
            parameter_names{end+1} = 'CY.p_2';
            lower_bound(end+1) = 1e3;
            upper_bound(end+1) = 1e9;
        end
    otherwise
        error('fit_crc_case:FitProfile', ...
            'Unknown CRC fitting profile: %s.', options.FitProfile);
end
x0 = cellfun(@(name) get_parameter(baseline_params, name), parameter_names);
if options.FitProfile == case_data.fit_profile
    if isempty(options.LowerBounds); options.LowerBounds = case_data.fit_lower_bounds; end
    if isempty(options.UpperBounds); options.UpperBounds = case_data.fit_upper_bounds; end
end
if ~isempty(options.LowerBounds)
    assert(numel(options.LowerBounds) == numel(x0), ...
        'fit_crc_case:LowerBounds', 'Provide one lower bound per fitted parameter.');
    lower_bound = options.LowerBounds;
end
if ~isempty(options.UpperBounds)
    assert(numel(options.UpperBounds) == numel(x0), ...
        'fit_crc_case:UpperBounds', 'Provide one upper bound per fitted parameter.');
    upper_bound = options.UpperBounds;
end
start_values = x0;
if ~isempty(fieldnames(options.InitialParameters))
    start_values = cellfun(@(name) get_parameter(options.InitialParameters,name), parameter_names);
end
if ~isempty(options.InitialValues)
    assert(numel(options.InitialValues) == numel(x0), ...
        'fit_crc_case:InitialValues', 'Provide one initial value per fitted parameter.');
    start_values = options.InitialValues;
end
% Baseline values and saved optima can lie on a bound; start strictly inside.
margin = 1e-6 * (upper_bound-lower_bound);
start_values = max(lower_bound+margin, min(upper_bound-margin,start_values));

data = case_data.observations(:, 2);
residual_count = numel(data);
pd1_index = find(strcmp(parameter_names,'CY.p_2'));

    function residual = objective(z)
        x = z .* x0;
        x(pd1_index) = 10 .^ z(pd1_index);
        try
            candidate = baseline_params;
            for j = 1:numel(parameter_names)
                candidate = set_parameter(candidate, parameter_names{j}, x(j));
            end
            evaluation = evaluate_crc_case(case_data, candidate);
            residual = evaluation.predicted(:) - data(:);
            if any(~isfinite(residual)) || any(evaluation.volume < 0)
                residual = repmat(1e6, residual_count, 1);
            end
        catch
            residual = repmat(1e6, residual_count, 1);
        end
    end

% Relative coordinates for rates/cell counts; logarithmic coordinates for p2.
start_coordinates = start_values ./ x0;
lower_coordinates = lower_bound ./ x0;
upper_coordinates = upper_bound ./ x0;
start_coordinates(pd1_index) = log10(start_values(pd1_index));
lower_coordinates(pd1_index) = log10(lower_bound(pd1_index));
upper_coordinates(pd1_index) = log10(upper_bound(pd1_index));
solver = bounded_least_squares(@objective, start_coordinates, ...
    lower_coordinates, upper_coordinates, ...
    NumStarts=options.NumStarts, MaxIterations=options.MaxIterations, ...
    Display=options.Display);
physical_x = solver.x .* x0;
physical_x(pd1_index) = 10 .^ solver.x(pd1_index);
physical_starts = solver.start_solutions .* x0;
physical_starts(:,pd1_index) = 10 .^ solver.start_solutions(:,pd1_index);
solver.x = physical_x;
solver.start_solutions = physical_starts;

optimized_params = baseline_params;
for j = 1:numel(parameter_names)
    optimized_params = set_parameter(optimized_params, parameter_names{j}, solver.x(j));
end

result.case_id = case_data.id;
result.parameter_names = parameter_names;
result.initial_values = x0;
result.optimizer_initial_values = start_values;
result.lower_bound = lower_bound;
result.upper_bound = upper_bound;
result.baseline_params = baseline_params;
result.optimized_params = optimized_params;
result.baseline = evaluate_crc_case(case_data, baseline_params);
result.optimized = evaluate_crc_case(case_data, optimized_params);
result.solver = solver;
result.objective = 'ordinary-volume-SSE';
result.fit_profile = options.FitProfile;
result.experimental_data = case_data.experimental_data;
end
