function result = bounded_least_squares(residual_function, x0, lower_bound, upper_bound, options)
%BOUNDED_LEAST_SQUARES Toolbox-aware bounded nonlinear least squares.
% Uses lsqnonlin when available; otherwise uses a bounded transformation and
% fminsearch, which is included in base MATLAB.

arguments
    residual_function (1,1) function_handle
    x0 (1,:) double
    lower_bound (1,:) double
    upper_bound (1,:) double
    options.NumStarts (1,1) double {mustBeInteger,mustBePositive} = 3
    options.MaxIterations (1,1) double {mustBeInteger,mustBePositive} = 150
    options.MaxFunctionEvaluations (1,1) double {mustBeInteger,mustBePositive} = 3000
    options.Display (1,:) char = 'off'
    options.RandomSeed (1,1) double = 20260817
    options.PreferLsqnonlin (1,1) logical = true
end

x0 = x0(:)';
lower_bound = lower_bound(:)';
upper_bound = upper_bound(:)';
if any(~isfinite([x0 lower_bound upper_bound])) || ...
        any(lower_bound >= upper_bound) || any(x0 <= lower_bound) || any(x0 >= upper_bound)
    error('bounded_least_squares:InvalidBounds', ...
        'Initial values must be finite and strictly inside finite bounds.');
end

has_lsqnonlin = options.PreferLsqnonlin && exist('lsqnonlin', 'file') == 2;
if has_lsqnonlin
    rng(options.RandomSeed, 'twister');
    solver_options = optimoptions('lsqnonlin', ...
        'Display', options.Display, ...
        'MaxIterations', options.MaxIterations, ...
        'MaxFunctionEvaluations', options.MaxFunctionEvaluations, ...
        'FiniteDifferenceStepSize',1e-4, ...
        'TypicalX',x0, ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8);
    transform = @(z) lower_bound + (upper_bound - lower_bound) ./ (1 + exp(-z));
    z0 = log((x0 - lower_bound) ./ (upper_bound - x0));
    best_resnorm = Inf;
    start_solutions = zeros(options.NumStarts, numel(x0));
    start_resnorm = zeros(options.NumStarts, 1);
    start_exitflag = zeros(options.NumStarts, 1);
    for start_index = 1:options.NumStarts
        if start_index == 1
            start = x0;
        else
            start = transform(z0 + 0.75 .* randn(size(z0)));
        end
        [candidate_x, candidate_resnorm, candidate_residual, candidate_exitflag, candidate_output] = ...
            lsqnonlin(residual_function, start, lower_bound, upper_bound, solver_options);
        start_solutions(start_index,:) = candidate_x;
        start_resnorm(start_index) = candidate_resnorm;
        start_exitflag(start_index) = candidate_exitflag;
        if candidate_resnorm < best_resnorm
            best_resnorm = candidate_resnorm;
            x = candidate_x;
            resnorm = candidate_resnorm;
            residual = candidate_residual;
            exitflag = candidate_exitflag;
            output = candidate_output;
        end
    end
    result.method = 'lsqnonlin-multistart';
else
    rng(options.RandomSeed, 'twister');
    transform = @(z) lower_bound + (upper_bound - lower_bound) ./ (1 + exp(-z));
    to_unbounded = @(x) log((x - lower_bound) ./ (upper_bound - x));
    z0 = to_unbounded(x0);
    solver_options = optimset( ...
        'Display', options.Display, ...
        'MaxIter', options.MaxIterations, ...
        'MaxFunEvals', options.MaxFunctionEvaluations, ...
        'TolFun', 1e-8, ...
        'TolX', 1e-8);

    best_value = Inf;
    best = struct();
    start_solutions = zeros(options.NumStarts, numel(x0));
    start_resnorm = zeros(options.NumStarts, 1);
    start_exitflag = zeros(options.NumStarts, 1);
    for start_index = 1:options.NumStarts
        if start_index == 1
            start = z0;
        else
            start = z0 + 0.75 .* randn(size(z0));
        end
        objective = @(z) sum(residual_function(transform(z)).^2);
        [z, value, exitflag, output] = fminsearch(objective, start, solver_options);
        start_solutions(start_index,:) = transform(z);
        start_resnorm(start_index) = value;
        start_exitflag(start_index) = exitflag;
        if value < best_value
            best_value = value;
            best.z = z;
            best.exitflag = exitflag;
            best.output = output;
        end
    end
    x = transform(best.z);
    residual = residual_function(x);
    resnorm = sum(residual.^2);
    exitflag = best.exitflag;
    output = best.output;
    result.method = 'fminsearch-bounded-transform';
end

result.x = x(:)';
result.residual = residual(:);
result.resnorm = resnorm;
result.exitflag = exitflag;
result.output = output;
result.start_solutions = start_solutions;
result.start_resnorm = start_resnorm;
result.start_exitflag = start_exitflag;
end
