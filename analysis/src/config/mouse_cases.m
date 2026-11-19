function cases = mouse_cases()
%MOUSE_CASES Mouse validation datasets and submitted manual parameter sets.
% Tumor volumes are in mm^3; model output is converted from cm^3 when fitted.
% Effective tumor-site concentration increments per anti-PD-1 dose (nM), converted from the
% 200 nM reference exposure using each schedule and a 223-hour half-life.
q_A = [44.75927982091127 44.75927982091127 104.4383195821263 29.83951988060751];

control_1 = [
     9 20.87130 22.77000 19.65200 26.156810 28.98030
    11 45.491390 67.85680 87.194340 28.533310 73.181590
    13 105.7660 78.911710 154.82880 63.379220 93.722110
    15 243.74710 169.42450 250.18030 181.23350 206.8960
    17 471.91490 364.9570 328.45930 436.46040 406.5600
    19 716.47490 512.83510 599.83340 609.11390 772.95070
    21 1069.0380 873.6200 897.24170 940.56620 835.0290];
anti_1 = [
     9 18.081500 22.770000 18.629850 22.195640 23.652000
    11 22.528000 58.176800 53.925090 68.295440 29.600000
    13 30.142840 78.911710 102.68950 116.43310 53.120000
    15 69.764370 120.19330 129.49080 153.56380 72.074130
    17 241.60920 283.32580 240.83940 313.94060 177.56800
    19 349.73690 512.83510 519.46930 391.78480 258.97640
    21 731.57850 938.57340 831.00350 594.23560 482.83250];
cases = make_case('LLC', 2.5e6, 9:23, [1 0.1 0.001 0.001 0.1 0.1 0.1], ...
    [0.52 4e8 20 75 100 q_A(1)], ...
    [make_group('Control', control_1(:,1), mean(control_1(:,2:end),2), std(control_1(:,2:end),0,2), Inf, NaN), ...
     make_group('Anti-PD-1', anti_1(:,1), mean(anti_1(:,2:end),2), std(anti_1(:,2:end),0,2), 1, NaN)]);

control_2 = [
     0 0 0
     3 24.2460487 16.6539513
     6 84.15634344 21.84365656
     9 170.2157198 32.78428022
    12 263.3938488 30.60615122
    15 349.4532251 54.54677488
    18 487.7979496 78.20205041
    21 659.4361384 87.5638616
    24 845.3246476 66.67535241
    27 1123.936352 36.06364801];
anti_2 = [
     0 0 0
     3 21.36075949 7.120253165
     6 85.44303797 18.98734177
     9 140.0316456 21.36075949
    12 182.7531646 23.73417722
    15 244.4620253 33.2278481
    18 310.5696203 35.60126582
    21 401.1075949 52.21518987
    24 562.5 73.57594937
    27 678.7974684 102.056962];
cases(2) = make_case('MC38', 2.4e6, 3:30, [1 0.01 0.001 0.001 0.01 0.1 0.1], ...
    [0.39 2.4e8 5 300 75 q_A(2)], ...
    [make_group('Control', control_2(:,1), control_2(:,2), control_2(:,3), Inf, NaN), ...
     make_group('Anti-PD-1', anti_2(:,1), anti_2(:,2), anti_2(:,3), 6, NaN)]);

control_3 = [1 60.6 40.4; 4 162.0 46.0; 7 296.0 47.0; 12 597.0 54.0; ...
    15 1030.0 100.0; 17 1470.0 160.0; 22 1720.0 130.0];
anti_3 = [1 60.6 20.0; 4 35.0 26.7; 7 76.1 20.0; 12 171.0 46.0; ...
    15 332.0 80.0; 17 386.0 100.0; 22 534.0 140.0];
cases(3) = make_case('A375', 8e6, 1:25, [1 0.01 0.001 0.001 0.01 0.01 0.01], ...
    [0.43 3.4e8 10 250 12 q_A(3)], ...
    [make_group('Control', control_3(:,1), control_3(:,2), control_3(:,3), Inf, NaN), ...
     make_group('Anti-PD-1', anti_3(:,1), anti_3(:,2), anti_3(:,3), 1, NaN)]);

case_1 = [0 0 0; 6 151.85 44.01; 9 323.50 77.02; 12 569.98 112.24; 15 900.09 110.04; 18 1254.40 213.47];
case_2 = [0 0 0; 6 110.04 46.21; 9 253.08 35.21; 12 354.31 61.62; 15 596.39 88.03; 18 853.87 103.43];
case_3 = [0 0 0; 6 110.04 39.61; 9 222.27 114.44; 12 420.33 63.82; 15 721.83 132.04; 18 897.89 176.06];
case_4 = [0 0 0; 6 41.81 17.61; 9 96.83 39.61; 12 145.25 46.21; 15 184.86 66.02; 18 233.27 72.62];
cases(4) = make_case('CT26', 2e6, 0:20, [1 0.01 0.01 0.01 0.01 0.01 0.01], ...
    [0.47 2.8e8 1.4 5 100 q_A(4)], ...
    [make_group('Control', case_1(:,1), case_1(:,2), case_1(:,3), Inf, 100), ...
     make_group('Low PD-L1', case_2(:,1), case_2(:,2), case_2(:,3), Inf, 1), ...
     make_group('Anti-PD-1', case_3(:,1), case_3(:,2), case_3(:,3), 4, 100), ...
     make_group('Low PD-L1 + Anti-PD-1', case_4(:,1), case_4(:,2), case_4(:,3), 4, 1)]);
% Injection days on the experimental time axis; MC38 receives only three doses.
% LLC: Chen et al., corresponding WO2025113697A1 protocol, day 9 then q3d.
% MC38: Zhang et al., days 7,10,13; A375: Rosato et al., weekly.
% CT26: Liu et al., q2d; treatment start retained from the supplied dataset.
intervals=[3,3,7,2];
dose_times={9:3:21,[7 10 13],1:7:22,3:2:19};
for case_index=1:numel(cases)
    cases(case_index).anti_interval_days=intervals(case_index);
    cases(case_index).anti_dose_times=dose_times{case_index};
end

end

function item = make_case(id, c0, time, ratios, manual, groups)
item.id = id;
item.Y0 = c0 .* ratios;
item.time = time;
item.manual.r_C = manual(1);
item.manual.C_max = manual(2);
item.manual.p_1 = manual(3);
item.manual.p_2 = manual(4);
item.manual.epC = manual(5);
item.manual.q_A = manual(6);
item.manual.T_40 = 4.88e8;
item.manual.T_80 = 2.39e8;
item.groups = groups;
end

function group = make_group(name, time, mean_value, sd, treatment_start_row, epC_override)
group.name = name;
group.time = time(:);
group.mean = mean_value(:);
group.sd = sd(:);
group.treatment_start_row = treatment_start_row;
group.epC_override = epC_override;
end
