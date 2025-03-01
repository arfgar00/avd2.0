clc; clear; close all;

%% -------------------- USER INPUTS -------------------- %%
% Define spanwise locations
spanwise_stations = linspace(0, 30, 50); % 30m span, 50 spanwise points

% Thickness exaggeration for better visibility
thickness_exaggeration = 5;  

% Define ply stackings for different regions (ensuring symmetry about midline)
% Format: {Span range, Ply angles}
PlyStackingRegions = {...
    [0, 10], [0, 45, -45, 90, 90, -45, 45, 0];  % 8 plies (Symmetric)
    [10, 20], [0, 90, 45, -45, 90, -45, 45, 0]; % 8 plies (Symmetric, Corrected)
    [20, 30], [0, 45, -45, 90, 90, -45, 45, 0, 0, 0]}; % 10 plies (Symmetric)

% Global stress state at each spanwise location (dummy values)
sigma_x_span = linspace(120, 30, length(spanwise_stations)); % Stress reducing spanwise
sigma_y_span = linspace(50, 15, length(spanwise_stations));
tau_xy_span  = linspace(25, 5, length(spanwise_stations));

% ---- Define Material Strengths by Orientation ----
Material_0deg  = [500, 250, 40, 100, 30];  % 0° Ply Properties
Material_45deg = [450, 220, 50, 120, 35];  % 45° and -45° Ply Properties
Material_90deg = [300, 150, 60, 130, 25];  % 90° Ply Properties

%% ----------------- TSAI-HILL FAILURE ANALYSIS ----------------- %%
FailureMatrix = {};
Y_PlyStacks = {};
Span_PlyStacks = {};

for j = 1:length(spanwise_stations) % Iterate over spanwise locations
    sigma_x = sigma_x_span(j);
    sigma_y = sigma_y_span(j);
    tau_xy = tau_xy_span(j);
    
    % Select appropriate ply stacking sequence for this span location
    for k = 1:length(PlyStackingRegions)
        region = PlyStackingRegions{k,1};
        if spanwise_stations(j) >= region(1) && spanwise_stations(j) < region(2)
            theta = PlyStackingRegions{k,2};
            break;
        end
    end
    
    % Ensure the stack is symmetric by mirroring the ply sequence if needed
    nPlies = length(theta);
    if mod(nPlies, 2) ~= 0
        error('Ply stack must have an even number of plies for symmetry.');
    end
    
    % Assign material properties based on ply orientation
    MaterialProps = zeros(nPlies, 5);
    for i = 1:nPlies
        if theta(i) == 0
            MaterialProps(i, :) = Material_0deg;
        elseif abs(theta(i)) == 45
            MaterialProps(i, :) = Material_45deg;
        elseif theta(i) == 90
            MaterialProps(i, :) = Material_90deg;
        end
    end
    
    % Compute Tsai-Hill for each ply
    FailureIndex = zeros(nPlies, 1);
    for i = 1:nPlies
        % Convert angle to radians
        theta_rad = deg2rad(theta(i));

        % Transformation matrix for stress (T-matrix)
        c = cos(theta_rad);
        s = sin(theta_rad);

        T = [ c^2, s^2,  2*c*s;
              s^2, c^2, -2*c*s;
             -c*s, c*s, c^2 - s^2 ];

        % Compute stress in the local 1-2 system
        sigma_global = [sigma_x; sigma_y; tau_xy]; 
        sigma_local = T * sigma_global;  % Transform stress

        % Extract transformed stress components
        sigma_1  = sigma_local(1); % Stress along fiber direction
        sigma_2  = sigma_local(2); % Transverse stress
        tau_12   = sigma_local(3); % In-plane shear

        % Get material strengths for this ply
        sigma_1t = MaterialProps(i,1);
        sigma_1c = MaterialProps(i,2);
        sigma_2t = MaterialProps(i,3);
        sigma_2c = MaterialProps(i,4);
        tau_12s  = MaterialProps(i,5);

        % Use appropriate strengths based on stress sign
        sigma_1_max = sigma_1t * (sigma_1 >= 0) + sigma_1c * (sigma_1 < 0);
        sigma_2_max = sigma_2t * (sigma_2 >= 0) + sigma_2c * (sigma_2 < 0);

        % ----------------- TSAI-HILL FAILURE CRITERION ----------------- %
        TsaiHillIndex = (sigma_1 / sigma_1_max)^2 ...
                      - (sigma_1 * sigma_2) / sigma_1_max^2 ...
                      + (sigma_2 / sigma_2_max)^2 ...
                      + (tau_12 / tau_12s)^2;

        % Store failure index
        FailureIndex(i) = TsaiHillIndex;
    end
    
    % Compute correct symmetric Y positions
    total_thickness = nPlies * 0.25 * thickness_exaggeration;
    Y_stack = linspace(-total_thickness/2, total_thickness/2, nPlies + 1);

    % Store results for this spanwise section
    FailureMatrix{j} = FailureIndex;
    Y_PlyStacks{j} = Y_stack;
    Span_PlyStacks{j} = spanwise_stations(j) * ones(size(Y_stack));
end

%% ----------------- PLOTTING FAILURE INDEX ----------------- %%
figure;
hold on;

% Get color map based on max failure value
maxFI = max(cellfun(@max, FailureMatrix)); % Maximum failure index for color scaling
minFI = 0; % Failure index lower bound
cmap = jet(256); % Gradient from blue to red
caxis([minFI, maxFI]); % Normalize colors

% Iterate through each spanwise location and plot
for j = 1:length(spanwise_stations)-1
    FailureIndex = FailureMatrix{j};
    Y_stack = Y_PlyStacks{j};
    X_stack = [spanwise_stations(j), spanwise_stations(j+1)];

    for i = 1:length(FailureIndex)
        % Normalize Tsai-Hill index to a color
        failureValue = FailureIndex(i);
        colorIdx = round((failureValue - minFI) / (maxFI - minFI) * 255) + 1;
        colorIdx = max(1, min(256, colorIdx)); % Ensure valid index
        color = cmap(colorIdx, :);

        % Create a rectangle for this ply section
        fill([X_stack(1), X_stack(2), X_stack(2), X_stack(1)], ...
             [Y_stack(i), Y_stack(i), Y_stack(i+1), Y_stack(i+1)], color, 'EdgeColor', 'k');
    end
end

% Labels and Formatting
xlabel('Span Location (m)');
ylabel('Thickness (mm) (Exaggerated)');
title('Tsai-Hill Failure Index Across Span and Thickness (Symmetric)');
colormap(jet);
colorbar; % Show color legend
grid on;
axis equal; % Keeps aspect ratio realistic
ylim([-max(cellfun(@max, Y_PlyStacks)), max(cellfun(@max, Y_PlyStacks))]); % Symmetric thickness
xlim([0, max(spanwise_stations)]); % Matches span range

hold off;
