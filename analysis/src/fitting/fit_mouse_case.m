function result = fit_mouse_case(case_data, options)
%FIT_MOUSE_CASE Jointly fit all arms by ordinary tumor-volume least squares.

arguments
    case_data (1,1) struct
    options.NumStarts (1,1) double {mustBeInteger,mustBePositive} = 3
    options.MaxIterations (1,1) double {mustBeInteger,mustBePositive} = 120
    options.Display (1,:) char = 'off'
    options.InitialParameters (1,1) struct = struct()
end

baseline_params = apply_mouse_manual_parameters(case_data);
% q_A is the fixed per-injection concentration increment; p_2 is fitted.
if strcmp(case_data.id, 'CT26')
    % Every CT26 arm explicitly fixes epC to either 100 or 1, so the global
    % epC value is not identifiable from this four-arm experiment.
    parameter_names = {'TSM.r_C','TSM.C_max','CY.p_1','CY.p_2'};
else
    parameter_names = {'TSM.r_C','TSM.C_max','CY.p_1','CY.p_2','TSM.epC'};
end
x0 = cellfun(@(name) get_parameter(baseline_params, name), parameter_names);
if strcmp(case_data.id, 'CT26')
    lower_bound = [0.20, max(1.02 * case_data.Y0(1), 0.5*x0(2)), 0.1*x0(3), 0.01*x0(4)];
    upper_bound = [0.80, 2.0*x0(2), 10*x0(3), 100*x0(4)];
else
    lower_bound = [0.20, max(1.02 * case_data.Y0(1), 0.5*x0(2)), 0.1*x0(3), 0.01*x0(4), 0.1*x0(5)];
    upper_bound = [0.80, 2.0*x0(2), 10*x0(3), 100*x0(4), 10*x0(5)];
end
baseline_evaluation = evaluate_mouse_case(case_data, baseline_params);
residual_count = numel(baseline_evaluation.observed);
initial_params = options.InitialParameters;
if isempty(fieldnames(initial_params))
    initial_params = baseline_params;
end
start_values = cellfun(@(name) get_parameter(initial_params,name),parameter_names);
margin = 1e-6 * (upper_bound - lower_bound);
start_values = max(lower_bound + margin,min(upper_bound - margin,start_values));
pd1_index = find(strcmp(parameter_names,'CY.p_2'));

    function z = to_coordinates(x)
        z = x ./ x0;
        z(pd1_index) = log10(x(pd1_index));
    end

    function x = to_parameters(z)
        x = z .* x0;
        x(:,pd1_index) = 10 .^ z(:,pd1_index);
    end

    function candidate = fitted_parameters(x)
        candidate = initial_params;
        for j = 1:numel(parameter_names)
            candidate = set_parameter(candidate,parameter_names{j},x(j));
        end
    end

    function residual = objective(z)
        try
            evaluation = evaluate_mouse_case(case_data,fitted_parameters(to_parameters(z)));
            residual = evaluation.predicted(:) - evaluation.observed(:);
            if any(~isfinite(residual)) || any(evaluation.predicted < 0)
                residual = repmat(1e6, residual_count, 1);
            end
        catch
            residual = repmat(1e6, residual_count, 1);
        end
    end

solver = bounded_least_squares(@objective,to_coordinates(start_values), ...
    to_coordinates(lower_bound),to_coordinates(upper_bound), ...
    NumStarts=options.NumStarts, MaxIterations=options.MaxIterations, ...
    Display=options.Display);

solver.x = to_parameters(solver.x);
solver.start_solutions = to_parameters(solver.start_solutions);
optimized_params = fitted_parameters(solver.x);

result.case_id = case_data.id;
result.fixed_q_A_nM = optimized_params.PK.q_A;
result.parameter_names = parameter_names;
result.initial_values = x0;
result.lower_bound = lower_bound;
result.upper_bound = upper_bound;
result.baseline_params = baseline_params;
result.optimized_params = optimized_params;
result.baseline = baseline_evaluation;
result.optimized = evaluate_mouse_case(case_data, optimized_params);
result.solver = solver;
result.objective = 'ordinary-volume-SSE';
end
