function cases = crc_cases()
%CRC_CASES Construct the eight CRC fitting cases from all valid observations.
% Rows whose mean/lower/upper tumor volumes are all zero are missing-value
% placeholders rather than biological observations and are excluded.

analysis_dir = fileparts(fileparts(fileparts(mfilename('fullpath'))));

source = load(fullfile(analysis_dir, 'CRC_fit', 'input_data', 'patient_all.mat'));
treatment = readtable(fullfile(analysis_dir, 'CRC_fit', ...
    'input_data', 'crc_treatment_schedule.csv'), 'TextType', 'string');

% Calibrated initial cancer-cell counts; fixed during each patient fit.
initial_cancer_cells = [17416920743.78532 31450000000 50000000000 30000000000 ...
    7327625826.5462618 7617928875.6086159 24700000000 6656970062.5397663];
% Clinical pathology values in combine.tex, tab:CRC_para; not fit targets.
ki67_percent = [70 90 NaN NaN 60 80 60 70];
fit_profiles = repmat("fixed_structure_pd1",1,8);
manual_values = [
    0.3671109302256212 100000000000.0 33.082024526724666 0.35000000000000003 0.4 0.03
    0.29158604406441996 100000000000.0 33.082024526724666 0.1 0.1 0.03
    0.261590273107821 100000000000.0 33.082024526724666 0.10715263548004622 0.13666920567122362 0.03
    0.2500246942789317 100000000000.0 33.082024526724666 0.10264455750105642 0.30142985340830813 0.03
    0.2500025599363937 100000000000.0 33.082024526724666 0.10037028526632884 0.318074893025754 0.03
    0.2500024372942667 100000000000.0 33.082024526724666 0.19997239365163003 0.20518288110012825 0.03
    0.2500000015871197 100000000000.0 33.082024526724666 0.1001380268941187 0.360733989868351 0.03
    0.25000000007108053 100000000000.0 33.082024526724666 0.12472529787688241 0.39311110466507493 0.03
    ];

manual_p1 = [1 1 1.0358998368980907 2.7336656877272816 1.7938317959049028 3.8654464125663486 8.096545835338363 9.450640823388223];
manual_p2 = [1000000.0 1000000.0 144928.2331144978 1000000.0 1000000.0 1000000.0 1000000.0 1000000.0];

% Daily treatment events exported from Supplementary Data.xlsx, Table 5.
% The historical MAT drug flags describe sustained effects, not injections.
cases = repmat(struct(), 1, 8);
for i = 1:8
    patient_matrix = source.(sprintf('P%d', i));
    valid_volume = all(isfinite(patient_matrix(:, 6:8)), 2) & ...
        any(patient_matrix(:, 6:8) ~= 0, 2);
    rows = find(valid_volume);
    raw_observations = patient_matrix(rows, [1 6 7 8]);
    interval_lower = min(raw_observations(:, 3), raw_observations(:, 4));
    interval_upper = max(raw_observations(:, 3), raw_observations(:, 4));
    cases(i).id = sprintf('P%02d', i);
    cases(i).ki67_percent = ki67_percent(i);
    cases(i).fit_profile = fit_profiles(i);
    cases(i).fit_lower_bounds = [];
    cases(i).fit_upper_bounds = [];
    schedule = treatment(treatment.case_id == string(cases(i).id), :);
    assert(isequal(schedule.study_day, patient_matrix(:,1)), ...
        'crc_cases:TreatmentDays', 'Treatment days do not match %s.', cases(i).id);
    cases(i).source_matrix = patient_matrix;
    cases(i).measurement_rows = rows;
    cases(i).excluded_zero_rows = find(~valid_volume & ...
        all(patient_matrix(:, 6:8) == 0, 2));
    cases(i).observations = [raw_observations(:, 1:2), interval_lower, interval_upper];
    cases(i).experimental_data = struct( ...
        'Y0', [initial_cancer_cells(i), repmat(1e6, 1, 6)], ...
        't', 0:1:max(patient_matrix(:, 1)), ...
        'Therapy', schedule{:,3:6}, ...
        'Volume', cases(i).observations);

    cases(i).manual.r_C = manual_values(i, 1);
    cases(i).manual.C_max = manual_values(i, 2);
    cases(i).manual.q_A = manual_values(i, 3);
    cases(i).manual.d_0 = manual_values(i, 4);
    cases(i).manual.lambda_0 = manual_values(i, 5);
    cases(i).manual.a = manual_values(i, 6);
    cases(i).manual.p_1 = manual_p1(i);
    cases(i).manual.p_2 = manual_p2(i);
end
cases(2).experimental_data.Y0([3 4]) = 5e5;
end
