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
b = 0.05; % (m)
h = 0.1; % (m)
c = 0.4; % (m)
a = 0.3; % (m)

% Set the range for z and initialize results
N = 4;
z_vals = linspace(span / 2, 0, N);
force_vals = zeros(N, length(n_values));
shear_vals = zeros(N, length(n_values));
moment_vals = zeros(N, length(n_values));
htailplane_mass = 2981.5;
vtailplane_mass = 2688.32;
htailplane_mass_dis = linspace(0, htailplane_mass, N);
vtailplane_mass_dis = linspace(0, vtailplane_mass, N);

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
        shear_vals(i, j) = real(double(subs(shear_result_numeric, z, span/2))) - real(double(subs(shear_result_numeric, z, z_current)));
        if i > 1
            moment_vals(i, j) = (shear_vals(i, j) + shear_vals(i - 1, j)) / (span/2/N) + moment_before;
            moment_before = moment_vals(i, j);
        end
    end
end


[x_generated, y_generated] = generate_interpolated_points(N)

cref = 4.24
a =  (0.6 + 0.12) * cref  / 2 - 0.25 * cref;
a = a / 2; % assume from spreadsheet that a = b1 = half of distance between singbox and load centre
root_c = y_generated;
b1 = a;


AoA = 7.5 - 0.15;
CM_a = -0.5;
CM_0 = 0.13;
Sref = 110.4;
Vau = 651;
M_0 = (AoA / 57.3 * CM_a + CM_0) * Sref * cref * 0.5 * 1.225 * (Vau / 3.6)^2
ht_mass = 2981;

lift_moment_dis = -root_c'.*flip(force_vals(:, 1),1)* a 

pitching_moment_dis = zeros(4,1);

for i = 1:size(force_vals,1)
    pitching_moment_dis(i,1) = M_0 * force_vals(i, 1) ./ (2*sum(force_vals(:,1))*0.1);
end
pitching_moment_dis = flip(pitching_moment_dis)

ht_weight_moment_dis(:,1) = ht_mass .* root_c .* b1 .* 9.81 * 3.75 .* span / 2

for i = 1:size(lift_moment_dis,1) - 1
    lift_moment_int(i) = (lift_moment_dis(i+1,1) + lift_moment_dis(i,1)) / (z_vals(i + 1) - z_vals(i))
    pitching_moment_int(i) = (pitching_moment_dis(i+1,1) + pitching_moment_dis(i,1)) / (z_vals(i + 1) - z_vals(i));
    ht_weight_moment_int(i) = (ht_weight_moment_dis(i+1,1) + ht_weight_moment_dis(i,1)) / (z_vals(i + 1) - z_vals(i));
    combined_moment_temp(i) = lift_moment_int(i) + pitching_moment_int(i)*1000 + ht_weight_moment_int(i);
end
combined_moment_temp = [flip(cumsum(flip(combined_moment_temp))), 0]./10^3;

% Define colors for different load factors
colors = ['b', 'g', 'r']; 
legend_entries = strcat('n = ', string(n_values));

% Create figure with subplots
figure;
set(gcf, 'Position', [100, 100, 800, 600]); % Adjust figure size

% Aerodynamic Load vs. Span
subplot(2, 2, 1);
hold on;
for j = 1:length(n_values)
    plot(z_vals, -force_vals(:, j), strcat(colors(j), '-o'), 'LineWidth', 2);
end
plot(z_vals, -htailplane_mass_dis * 9.81 * 3.75 / 60, 'm-o', 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 2);
title('Aerodynamic Load vs. Span');
xlabel('Span station (m)', 'Interpreter', 'latex');
ylabel('Force (N)', 'Interpreter', 'latex');
grid on;
legend([legend_entries, 'Tailplane Mass', 'y = 0 Reference'], 'Location', 'Best');

% Shear Force vs. Span
subplot(2, 2, 2);
hold on;
for j = 1:length(n_values)
    plot(z_vals, shear_vals(:, j), strcat(colors(j), '-o'), 'LineWidth', 2);
end
yline(0, 'k--', 'LineWidth', 2);
title('Shear Force vs. Span');
xlabel('Span station (m)', 'Interpreter', 'latex');
ylabel('Shear (N)', 'Interpreter', 'latex');
grid on;
legend([legend_entries, 'y = 0 Reference'], 'Location', 'Best');

% Bending Moment vs. Span
subplot(2, 2, 3);
hold on;
for j = 1:length(n_values)
    plot(z_vals, moment_vals(:, j), strcat(colors(j), '-o'), 'LineWidth', 2);
end
yline(0, 'k--', 'LineWidth', 2);
title('Bending Moment vs. Span');
xlabel('Span station (m)', 'Interpreter', 'latex');
ylabel('Moment (Nm)', 'Interpreter', 'latex');
grid on;
legend([legend_entries, 'y = 0 Reference'], 'Location', 'Best');

% Torque vs. Span
subplot(2, 2, 4);
plot(z_vals, flip(combined_moment_temp), 'r-o', 'LineWidth', 2);
hold on;
yline(0, 'k--', 'LineWidth', 2);
title('Torque vs. Span');
xlabel('Span station (m)', 'Interpreter', 'latex');
ylabel('Torque (N·m)', 'Interpreter', 'latex');
grid on;

% Set global axis properties
set(gca, 'FontSize', 11, 'LineWidth', 1.25, 'Box', 'off');




% Define ply types, thickness, and number of plies
plies = [
    struct('ply_type', 'deg0', 'ply_thickness', 0.0005, 'ply_count', 20, 'stack', struct());
    struct('ply_type', 'deg90', 'ply_thickness', 0.0005, 'ply_count', 4, 'stack', struct());
    struct('ply_type', 'deg45', 'ply_thickness', 0.0005, 'ply_count', 8, 'stack', struct());
];


plies_station_I = [
    struct('ply_type','deg0','ply_thickness', 0.0005, 'ply_count', 36, 'span_lower', 0, 'span_upper', 10,'stack', 0),
    struct('ply_type','deg90', 'ply_thickness', 0.0005, 'ply_count', 8, 'span_lower', 0, 'span_upper', 10,'stack', 0),
    struct('ply_type','deg45', 'ply_thickness', 0.0005, 'ply_count', 8, 'span_lower', 0, 'span_upper', 10,'stack', 0) 
    ]

plies_station_II = [...
    struct('ply_type','deg0','ply_thickness', 0.0005, 'ply_count', 14, 'span_lower', 10, 'span_upper', 15,'stack', 0),
    struct('ply_type','deg90', 'ply_thickness', 0.0005, 'ply_count', 8, 'span_lower', 10, 'span_upper', 15,'stack', 0),
    struct('ply_type','deg45', 'ply_thickness', 0.0005, 'ply_count', 8, 'span_lower', 10, 'span_upper', 15,'stack', 0)
    ]

plies_station_III = [
    struct('ply_type','deg0','ply_thickness', 0.0005, 'ply_count', 8, 'span_lower', 15, 'span_upper', 20,'stack', 0),
    struct('ply_type','deg90', 'ply_thickness', 0.0005, 'ply_count', 8, 'span_lower', 15, 'span_upper', 20,'stack', 0), 
    struct('ply_type','deg45', 'ply_thickness', 0.0005, 'ply_count', 10, 'span_lower', 15, 'span_upper', 20,'stack', 0) 
]




% Create the main struct to hold all stations
plies_struct = struct( ...
    'station_I', plies_station_I, ...  % Assign plies data for station I
    'station_II', plies_station_II, ... % Assign plies data for station II
    'station_III', plies_station_III ... % Assign plies data for station III
);

% Get the field names of the struct (i.e., station names)
stations = fieldnames(plies_struct);

% Determine the number of stations
N = size(stations,1); % Example number of elements

% Initialize a struct array to store stack data for each station
stack = repmat(struct('data', []), N, 1); 

% Loop through each station
for l = 1:size(stations,1) 

    % Initialize a struct to track the number of plies in each fiber direction
    ply_count = struct('deg0', 0, 'deg90', 0, 'deg45', 0); 

    % Re-fetch the station names in case of updates
    stations = fieldnames(plies_struct); 

    % Retrieve the current station's plies data
    current_station = plies_struct.(stations{l}); 

    % Initialize a temporary stack variable to store plies for the current station
    temp_stack = []; 

    % Track the remaining number of plies to be added in each direction
    remaining_ply_count = [ ...
        plies_struct.(stations{l})(1).ply_count, ... % Number of 0-degree plies
        plies_struct.(stations{l})(2).ply_count, ... % Number of 90-degree plies
        plies_struct.(stations{l})(3).ply_count  ... % Number of 45-degree plies
    ];

    % Define the symmetry condition based on the total number of plies
    total_layers = sum(remaining_ply_count); 

    % Extract the span lower values for each ply in the station
    span_lower_values = [current_station.span_lower]; 

    % Extract the span upper values for each ply in the station
    span_upper_values = [current_station.span_upper]; 

    % Loop until all plies are added
    while remaining_ply_count(1) > 0 || remaining_ply_count(2) > 0 || remaining_ply_count(3) > 0  
        % Continue looping until all plies are placed

        added_ply = false;  % Initialize a flag to track if a ply was added in the current iteration

        % Add 45-degree plies (ensuring symmetry)
        if remaining_ply_count(3) > 0  % Check if there are any 45-degree plies left to add
            
            if size(temp_stack) == 0  % If the stack is empty, place the first ply symmetrically
                % Create a 45-degree ply structure
                ply = struct('ply_type', 'deg45', 'ply_thickness', 0.0005, ...
                            'lower_span', span_lower_values(3), 'upper_span', span_upper_values(3));
                
                % Insert ply at both ends of the stack to maintain symmetry
                temp_stack = [ply, temp_stack, ply];  
                
                % Update the count of 45-degree plies added
                ply_count.deg45 = ply_count.deg45 + 2;
                
                % Reduce the number of remaining 45-degree plies
                remaining_ply_count(3) = remaining_ply_count(3) - 2;
                
                added_ply = true;  % Mark that a ply was added in this iteration

            elseif remaining_ply_count(3) > 0 && length(temp_stack) >= 2 && ...
                   strcmp(temp_stack(end-1).ply_type, 'deg45') && strcmp(temp_stack(end).ply_type, 'deg45')
                % If at least two plies exist in the stack and both are 45-degree, add a 90-degree ply

                % Create a 90-degree ply structure
                ply = struct('ply_type', 'deg90', 'ply_thickness', 0.0005, ...
                            'lower_span', span_lower_values(2), 'upper_span', span_upper_values(2));

                % Insert ply symmetrically at both ends of the stack
                temp_stack = [ply, temp_stack, ply];  
                
                if remaining_ply_count(3) > 2  % If more than 2 plies remain
                    ply_count.deg90 = ply_count.deg90 + 2;  % Update the 90-degree ply count
                    remaining_ply_count(2) = remaining_ply_count(2) - 2;  % Reduce 90-degree ply count
                    added_ply = true;
                else  % If 2 or fewer 45-degree plies remain
                    ply_count.deg90 = ply_count.deg90 + 2;  % Update the 90-degree ply count
                    added_ply = true;
                end 

            elseif remaining_ply_count(3) >= 4  % If there are at least 4 remaining 45-degree plies
                % Create a 45-degree ply structure
                ply = struct('ply_type', 'deg45', 'ply_thickness', 0.0005, ...
                            'lower_span', span_lower_values(3), 'upper_span', span_upper_values(3));

                % Insert two plies symmetrically at both ends
                temp_stack = [ply, ply, temp_stack, ply, ply];  
                
                % Update the count of 45-degree plies added
                ply_count.deg45 = ply_count.deg45 + 4;
                
                % Reduce the number of remaining 45-degree plies
                remaining_ply_count(3) = remaining_ply_count(3) - 4;
                
                added_ply = true;  

            else  % If there are fewer than 4 remaining 45-degree plies
                remaining_ply_count(3) == 0;  % Ensure that no more 45-degree plies remain
                
                % Create a 45-degree ply structure
                ply = struct('ply_type', 'deg45', 'ply_thickness', 0.0005, ...
                            'lower_span', span_lower_values(3), 'upper_span', span_upper_values(3));

                % Insert the last remaining 45-degree ply symmetrically
                temp_stack = [ply, temp_stack, ply];  
                
                % Update the count of 45-degree plies added
                ply_count.deg45 = ply_count.deg45 + 2;
                
                % Reduce the number of remaining 45-degree plies
                remaining_ply_count(3) = remaining_ply_count(3) - 2;
                
                added_ply = true;
            end
        end

    
        % Add 90 degree plies (1 at a time for symmetry)
        if remaining_ply_count(2) > 0  % Check if there are 90-degree plies left to add
            % Place ply symmetrically
            ply = struct('ply_type', 'deg90', 'ply_thickness', 0.0005, 'lower_span', span_lower_values(2), 'upper_span', span_upper_values(2));
            temp_stack = [ply, temp_stack, ply];  % Insert at both ends for symmetry
            ply_count.deg90 = ply_count.deg90 + 2;  % Update the count of 90-degree plies
            remaining_ply_count(2) = remaining_ply_count(2) - 2;  % Decrease the count of remaining 90-degree plies
            added_ply = true;  % Flag that a ply was added
    
            % Always add a deg0 ply between deg90 and deg45 if needed
            if remaining_ply_count(3) > 0 && length(temp_stack) >= 2 && strcmp(temp_stack(end-1).ply_type, 'deg90') && strcmp(temp_stack(end).ply_type, 'deg45')
                % Create a 0-degree ply to be placed between 90-degree and 45-degree plies
                ply = struct('ply_type', 'deg0', 'ply_thickness', 0.0005, 'lower_span', span_lower_values(1), 'upper_span', span_upper_values(1));
                temp_stack = [ply, temp_stack, ply];  % Insert at both ends for symmetry
                ply_count.deg0 = ply_count.deg0 + 2;  % Update the count of 0-degree plies
                remaining_ply_count(1) = remaining_ply_count(1) - 2;  % Decrease the count of remaining 0-degree plies
                added_ply = true;  % Flag that a ply was added
            end
        end
    
        % Add 0 degree plies (1 at a time for symmetry)
        while remaining_ply_count(1) > 0  % Continue while there are 0-degree plies left to add
            % Check if adding another 0° ply will create more than 2 consecutive 0° plies
            if length(temp_stack) >= 2 && strcmp(temp_stack(end-1).ply_type, 'deg0') && strcmp(temp_stack(end).ply_type, 'deg0')
                if remaining_ply_count(3) == 0  % If no 45-degree plies are left
                    ply = struct('ply_type', 'deg45', 'ply_thickness', 0.0005, 'lower_span', span_lower_values(2), 'upper_span', span_upper_values(2));
                    temp_stack = [ply, temp_stack, ply];  % Insert at both ends for symmetry
                    ply_count.deg45 = ply_count.deg45 + 2;  % Update the count of 90-degree plies
                    % remaining_ply_count(2) = remaining_ply_count(2) - 2; % Commented out, possibly intentional
                    added_ply = true;  % Flag that a ply was added
                else 
                    ply = struct('ply_type', 'deg45', 'ply_thickness', 0.0005, 'lower_span', span_lower_values(2), 'upper_span', span_upper_values(2));
                    temp_stack = [ply, temp_stack, ply];  % Insert at both ends for symmetry
                    ply_count.deg45 = ply_count.deg45 + 2;  % Update the count of 90-degree plies
                    remaining_ply_count(2) = remaining_ply_count(2) - 2;  % Decrease the count of remaining 90-degree plies
                    added_ply = true;  % Flag that a ply was added
                end
            else
                % Place a 0-degree ply symmetrically
                ply = struct('ply_type', 'deg0', 'ply_thickness', 0.0005, 'lower_span', span_lower_values(1), 'upper_span', span_upper_values(1));
                temp_stack = [ply, temp_stack, ply];  % Insert at both ends for symmetry
                ply_count.deg0 = ply_count.deg0 + 2;  % Update the count of 0-degree plies
                remaining_ply_count(1) = remaining_ply_count(1) - 2;  % Decrease the count of remaining 0-degree plies
                added_ply = true;  % Flag that a ply was added
            end
        end
            
        % If no ply was added, break the loop to avoid an infinite loop
        if ~added_ply
            break
        end
    end

    % Assign the temporary stack to the correct index
    stack(l).data = temp_stack; % Store the final 1x40 struct array in the stack

    % Display results
    disp('Final stack of plies (symmetrical):');
    for i = 1:length(temp_stack)  % Loop through each ply in the stack
        disp(['Ply ', num2str(i), ': ', temp_stack(i).ply_type, ' Thickness: ', num2str(temp_stack(i).ply_thickness)]);
    end

    % Display the final count of plies by direction
    disp('Ply counts by direction:');
    disp(ply_count);
end  


% Define colors for each ply orientation
colors = struct('deg0', [0, 0.5, 0], ...   % Green for 0°
                'deg45', [1, 0, 0], ...    % Red for 45°
                'deg90', [0, 0, 1]);       % Blue for 90°

% Initialize bottom of the stack
bottom = 0;
total_thickness = 0; % To keep track of total thickness



% Iterate through the stack and plot each ply
for i = 1:length(stack)
    % Create figure
    figure;
    hold on;
    title('Ply Stack Visualization');
    xlabel('Span (m)');
    ylabel('Thickness (m)');
    legendHandles = []; % Store patch handles for legend
    legendEntries = {}; % Store labels for legend
    y = [];
    x = [];
    h = [];
    bottom = 0;
    total_thickness = 0;
    for j = 1:length(stack(i).data)
        ply_type = stack(i).data(j).ply_type; % Access ply_type using dot notation
        ply_thickness = stack(i).data(j).ply_thickness; % Access ply_thickness using dot notation
        ply_lowerspan = stack(i).data(j).lower_span;
        ply_upperspan = stack(i).data(j).upper_span;

        % Define x-coordinates for the ply
        x = [ply_lowerspan ply_upperspan ply_upperspan ply_lowerspan];

        % Define y-coordinates for the ply
        y = [bottom, bottom, bottom + ply_thickness, bottom + ply_thickness];
    
        % Fill rectangle with corresponding color
        h = fill(x, y, colors.(ply_type), 'EdgeColor', 'k');
    
        % Update bottom position for next ply
        bottom = bottom + ply_thickness;
        total_thickness = total_thickness + ply_thickness; % Update total thickness
    
        % Store legend entry if not already added
        if ~ismember(ply_type, legendEntries)
            legendHandles(end+1) = patch(NaN, NaN, colors.(ply_type)); % Create invisible patch
            legendEntries{end+1} = ply_type;
        end  
    end

    % Check if total thickness is even and plot a midway line if true
    if mod(length(stack(i).data), 2) == 0
        midway_height = total_thickness / 2;  % Midway point
        plot([0, 20], [midway_height, midway_height], 'k--', 'LineWidth', 2); % Midway line
    end
    
    % Add legend with correct colors
    legend(legendHandles, legendEntries, 'Location', 'northeast');
    
    % Adjust plot limits
    ylim([0, bottom]);
    xlim([ply_lowerspan, ply_upperspan]);
    hold off;
end

% Loop through each column and perform the division
for i = 1:size(Nx, 1)
    D11 = 5.14e6;
    D22 = 4.4e6;
    D12 = 9.62e5;
    D33 = 1.8e6;
    Ncr_compressive = Nx(i) * 10^3 / (2 * pi^2 / b^2 * ( (D11 * D22)^0.5 + D12 + 2 * D33)) ;
    
    % Display the result
    fprintf('Ncr_compressive for column %d = %.4f\n', i, Ncr_compressive);

    % Check if it exceeds 1
    if Ncr_compressive > 1
        warning('Warning! The structure will buckle under compression.');
        fprintf('Ncr_compressive for column %d = %.4f\n', i, Ncr_compressive);
    else
        fprintf('Ncr_compressive for column %d = %.4f, will survive loading!\n', i, Ncr_compressive);
    end
end