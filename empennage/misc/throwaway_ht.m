clear; clc; clf;
% Define symbolic variable
syms z L_vs b_VS

% Define Ldash_vs as a symbolic expression
Ldash_vs = (4 * L_vs / (pi * b_VS)) * sqrt(1 - (z / b_VS)^2);

% Perform symbolic integration
shear_result = int(Ldash_vs, z);
moment_result = int(shear_result, z);

% Set values for L and span
L = 89039;
span = 26.41;
n_values = [3.75, -1.5];
colors = ['b', 'g', 'r']; % Colors for different n-values
% panel size and loading
omega = L / span / 2; % (N/m)

% Set the range for z and initialize results
N = 100;
z_vals = linspace(span / 2, 0, N);
force_vals = zeros(N, length(n_values));
shear_vals = zeros(N, length(n_values));
moment_vals = zeros(N, length(n_values));
htailplane_mass = 2981.5;
vtailplane_mass = 2688.32;
htailplane_mass_dis = htailplane_mass / N * ones(1, N);
vtailplane_mass_dis = linspace(0, vtailplane_mass, N);

% Modify the htailplane_mass_dis for points with span < 3
for i = 1:N
    if z_vals(i) < 3
        reduction_factor = -53198.2 * 5 / sum(z_vals < 3);  % Calculate the reduction to be evenly distributed
        htailplane_mass_dis(i) = htailplane_mass_dis(i) - reduction_factor;
    end
end

% Loop over different load factors
for j = 1:length(n_values)
    n = n_values(j);
    force_result_numeric = subs(Ldash_vs, [L_vs, b_VS], [L * n, span / 2]);
    shear_result_numeric = subs(shear_result, [L_vs, b_VS], [L * n, span / 2]);
    moment_result_numeric = subs(moment_result, [L_vs, b_VS], [L * n, span / 2]);
    moment_before = 0;
    
    for i = 1:N
        z_current = z_vals(i);
        force_vals(i, j) = real(double(subs(force_result_numeric, z, span/2))) - real(double(subs(force_result_numeric, z, z_current)));    
        shear_vals(i, j) = real(double(subs(shear_result_numeric, z, span/2))) - real(double(subs(shear_result_numeric, z, z_current))) - htailplane_mass_dis(i);
        if i > 1
            moment_vals(i, j) = (shear_vals(i, j) + shear_vals(i - 1, j))/ 2 + moment_before ;
            moment_before = moment_vals(i, j);
        end
    end
end


[x_generated, y_generated] = generate_interpolated_points(N);

cref = 4.24;
a =  (0.6 + 0.12) * cref / 2 - 0.25 * cref;
a = a / 2; % assume from spreadsheet that a = b1 = half of distance between singbox and load centre
root_c = y_generated;
b1 = a;

AoA = 7.5 - 0.15;
CM_a = -0.5;
CM_0 = 0.13;
Sref = 110.4;
Vau = 651;
M_0 = (AoA / 57.3 * CM_a + CM_0) * Sref * cref * 0.5 * 1.225 * (Vau / 3.6)^2;
ht_mass = 2981;

M_0_2 = (AoA / 57.3 * CM_a + CM_0) * Sref * cref * 0.5 * 1.225 * (500 / 3.6)^2;

lift_moment_dis = -root_c'.*flip(force_vals(:, 1),1)* a;
lift_moment_dis_col2 = -root_c'.*flip(force_vals(:, 2),1)*a;

pitching_moment_dis = zeros(4,1);
pitching_moment_dis_col2 = zeros(4,1);

for i = 1:size(force_vals,1)
    pitching_moment_dis(i,1) = M_0 * force_vals(i, 1) ./ (2*sum(force_vals(:,1))*0.1);
    pitching_moment_dis_col2(i,1) = M_0_2 * force_vals(i, 2) ./ (2*sum(force_vals(:,2))*0.1);
end

pitching_moment_dis = flip(pitching_moment_dis);
pitching_moment_dis_col2 = flip(pitching_moment_dis_col2);

ht_weight_moment_dis(:,1) = ht_mass .* root_c .* b1 .* 9.81 * 3.75 .* span / 2;
ht_weight_moment_dis_col2(:,1) = ht_mass .* root_c .* b1 .* 9.81 * 3.75 .* span / 2;

for i = 1:size(lift_moment_dis,1) - 1
    % Column 1 - Lift, Pitching, and Tailplane Weight Moments
    lift_moment_int(i) = (lift_moment_dis(i+1,1) + lift_moment_dis(i,1)) / (z_vals(i + 1) - z_vals(i));
    pitching_moment_int(i) = (pitching_moment_dis(i+1,1) + pitching_moment_dis(i,1)) / (z_vals(i + 1) - z_vals(i));
    ht_weight_moment_int(i) = (ht_weight_moment_dis(i+1,1) + ht_weight_moment_dis(i,1)) / (z_vals(i + 1) - z_vals(i));
    combined_moment_temp(i) = lift_moment_int(i) + pitching_moment_int(i)*1000 + ht_weight_moment_int(i);
    
    % Column 2 - Lift, Pitching, and Tailplane Weight Moments for the second load factor
    lift_moment_int_col2(i) = (lift_moment_dis_col2(i+1,1) + lift_moment_dis_col2(i,1)) / (z_vals(i + 1) - z_vals(i));
    pitching_moment_int_col2(i) = (pitching_moment_dis_col2(i+1,1) + pitching_moment_dis_col2(i,1)) / (z_vals(i + 1) - z_vals(i));
    ht_weight_moment_int_col2(i) = (ht_weight_moment_dis_col2(i+1,1) + ht_weight_moment_dis_col2(i,1)) / (z_vals(i + 1) - z_vals(i));
    combined_moment_temp_col2(i) = lift_moment_int_col2(i) + pitching_moment_int_col2(i)*1000 + ht_weight_moment_int_col2(i);
end
combined_moment_temp = [flip(cumsum(flip(combined_moment_temp))), 0]./10^3;
combined_moment_temp_col2 = [flip(cumsum(flip(combined_moment_temp_col2))), 0]./10^3;


% Define colors for different load factors
colors = ['b', 'g', 'r']; 
legend_entries = strcat('n = ', string(n_values));

% Create figure with subplots
figure;
set(gcf, 'Position', [100, 100, 800, 600]); % Adjust figure size

% Define the aspect ratio
width = 800;  % Example width (you can adjust as needed)
height = round(width / 2);  % Maintain 2:1 aspect ratio

% Aerodynamic Load vs. Span
figure;
set(gcf, 'Position', [100, 100, width, height]); % Set figure size
hold on;
for j = 1:length(n_values)
    plot(z_vals, force_vals(:, j), strcat(colors(j), '-o'), 'LineWidth', 2);
end
plot(z_vals, -htailplane_mass_dis, 'm-o', 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 2);
title('Aerodynamic Load vs. Span');
xlabel('Span station (m)', 'Interpreter', 'latex');
ylabel('Force (N)', 'Interpreter', 'latex');
grid on;
legend([legend_entries, 'Tailplane Mass', 'y = 0 Reference'], 'Location', 'Best');
% Set global axis properties
set(gca, 'FontSize', 22, 'LineWidth', 1.25, 'Box', 'off');
saveas(gcf, 'aerodynamic_load_vs_span.png'); 

% Shear Force vs. Span
figure;
set(gcf, 'Position', [100, 100, width, height]); % Set figure size
hold on;
for j = 1:length(n_values)
    plot(z_vals, -shear_vals(:, j), strcat(colors(j), '-o'), 'LineWidth', 2);
end
yline(0, 'k--', 'LineWidth', 2);
title('Shear Force vs. Span');
xlabel('Span station (m)', 'Interpreter', 'latex');
ylabel('Shear (N)', 'Interpreter', 'latex');
grid on;
legend([legend_entries, 'y = 0 Reference'], 'Location', 'Best');
% Set global axis properties
set(gca, 'FontSize', 22, 'LineWidth', 1.25, 'Box', 'off');
saveas(gcf, 'shear_vs_span.png'); 


% Bending Moment vs. Span
figure;
set(gcf, 'Position', [100, 100, width, height]); % Set figure size
hold on;
for j = 1:length(n_values)
    plot(z_vals, -moment_vals(:, j), strcat(colors(j), '-o'), 'LineWidth', 2);
end
yline(0, 'k--', 'LineWidth', 2);
title('Bending Moment vs. Span');
xlabel('Span station (m)', 'Interpreter', 'latex');
ylabel('Moment (Nm)', 'Interpreter', 'latex');
grid on;
legend([legend_entries, 'y = 0 Reference'], 'Location', 'Best');
% Set global axis properties
set(gca, 'FontSize', 22, 'LineWidth', 1.25, 'Box', 'off');
saveas(gcf, 'bending_vs_span.png'); 

% Torque vs. Span
figure;
set(gcf, 'Position', [100, 100, width, height]); % Set figure size
plot(z_vals, flip(combined_moment_temp), 'b-o', 'LineWidth', 2);
hold on;
plot(z_vals, -flip(combined_moment_temp_col2), 'g-o', 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 2);
title('Torque vs. Span');
xlabel('Span station (m)', 'Interpreter', 'latex');
ylabel('Torque (Nm)', 'Interpreter', 'latex');
grid on;
legend([legend_entries, 'y = 0 Reference'], 'Location', 'Best');
% Set global axis properties
set(gca, 'FontSize', 22, 'LineWidth', 1.25, 'Box', 'off');
saveas(gcf, 'torque_vs_span.png'); 




%change direction of lift distribution

function [x_generated, y_generated] = generate_interpolated_points(n)
    % Your data (x is the second column and y is the first column)
    x = [
        0.00
        0.50
        1.00
        1.50
        2.00
        2.50
        3.00
        3.50
        4.00
        4.50
        5.00
        5.50
        6.00
        6.50
        7.00
        7.50
        8.00
        8.50
        9.00
        9.50
        10.00
        10.50
        11.00
        11.50
        12.00
        12.50
        13.00
        13.50
        14.00
        14.50
        15.00
        15.50
        16.00
        16.50
        17.00
        17.50
        18.00
        18.50
        19.00
        19.50
        20.00
        20.50
        21.00
        21.50
        22.00
        22.50
        23.00
        23.50
        24.00
        24.50
        25.00
        25.50
        26.00
        26.50
        27.00
        27.50
        28.00
        28.50
        29.00
        29.30
    ];

    y = [
        1.00
        0.979093099
        0.958186199
        0.937279298
        0.916372397
        0.895465496
        0.874558596
        0.853651695
        0.832744794
        0.8132409
        0.801499214
        0.790487652
        0.77947609
        0.768464528
        0.757452965
        0.746441403
        0.735429841
        0.724418279
        0.713406716
        0.702395154
        0.691383592
        0.68037203
        0.669360467
        0.658348905
        0.647337343
        0.636325781
        0.625314219
        0.614302656
        0.603291094
        0.592279532
        0.58126797
        0.570256407
        0.559244845
        0.548233283
        0.537221721
        0.526210158
        0.515198596
        0.504187034
        0.493175472
        0.482163909
        0.471152347
        0.460140785
        0.449129223
        0.438117661
        0.427106098
        0.416094536
        0.405082974
        0.394071412
        0.383059849
        0.372048287
        0.361036725
        0.350025163
        0.3390136
        0.328002038
        0.316990476
        0.305978914
        0.294967351
        0.283955789
        0.272944227
        0.00
    ];

    % Normalize the x-values such that the maximum value is 1
    x_normalized = x / max(x);

    % Create a finer x grid for interpolation (for example, 100 points between 0 and 1)
    x_fine = linspace(min(x_normalized), max(x_normalized), 100);

    % Interpolate y-values at the new normalized x positions
    y_interpolated = interp1(x_normalized, y, x_fine, 'pchip');  % pchip for smooth interpolation

    % To generate 'n' points, for example n = 10:
    x_generated = linspace(min(x_normalized), max(x_normalized), n);
    y_generated = interp1(x_normalized, y, x_generated, 'pchip');

    % Plot the generated points (optional)
    figure;
    plot(x_generated, y_generated, 'rx', 'MarkerSize', 10); % Generated points
    title('Generated Interpolated Points');
    xlabel('Normalized x');
    ylabel('y');
    grid on;
end