function evaluation = evaluate_mouse_case(case_data, params)
%EVALUATE_MOUSE_CASE Simulate every experimental arm for one mouse dataset.

n_steps = numel(case_data.time) - 1;
evaluation.groups = repmat(struct(), 1, numel(case_data.groups));
all_observed = [];
all_predicted = [];

for group_index = 1:numel(case_data.groups)
    group = case_data.groups(group_index);
    group_params = params;
    if isfinite(group.epC_override)
        group_params.TSM.epC = group.epC_override;
    end

    therapy = zeros(n_steps, 4);
    if isfinite(group.treatment_start_row)
        dose_rows=round(case_data.anti_dose_times-case_data.time(1))+1;
        dose_rows=dose_rows(dose_rows>=1 & dose_rows<=n_steps);
        therapy(dose_rows,4)=1;
    end
    experimental_data = struct('Y0', case_data.Y0, 't', case_data.time, 'Therapy', therapy);
    [time, volume] = run_simulation(experimental_data, group_params);
    volume_mm3 = volume .* 1e3;

    valid = group.time >= min(time) & group.time <= max(time) & isfinite(group.mean);
    predicted = interp1(time, volume_mm3, group.time(valid), 'linear');
    observed = group.mean(valid);

    evaluation.groups(group_index).name = group.name;
    evaluation.groups(group_index).time = time;
    evaluation.groups(group_index).volume_mm3 = volume_mm3;
    evaluation.groups(group_index).data_time = group.time(valid);
    evaluation.groups(group_index).observed = observed;
    evaluation.groups(group_index).predicted = predicted;

    all_observed = [all_observed; observed(:)]; %#ok<AGROW>
    all_predicted = [all_predicted; predicted(:)]; %#ok<AGROW>
end

evaluation.observed = all_observed;
evaluation.predicted = all_predicted;
evaluation.metrics = fit_metrics(all_observed, all_predicted);
end

