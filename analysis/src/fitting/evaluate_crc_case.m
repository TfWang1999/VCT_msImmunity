function evaluation = evaluate_crc_case(case_data, params)
%EVALUATE_CRC_CASE Simulate one CRC patient and evaluate measurement fit.

[time, volume, states, state_counts] = run_simulation(case_data.experimental_data, params);
indices = round(case_data.observations(:, 1));
if any(indices < 1) || any(indices > numel(volume))
    error('evaluate_crc_case:MeasurementIndex', ...
        'Measurement index is outside the simulated trajectory for %s.', case_data.id);
end

observed = case_data.observations(:, 2);
predicted = volume(indices);
evaluation.time = time;
evaluation.volume = volume;
evaluation.states = states;
evaluation.state_counts = state_counts;
evaluation.observed = observed;
evaluation.predicted = predicted;
evaluation.metrics = fit_metrics(observed, predicted);
end

