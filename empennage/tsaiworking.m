
clc; close all;

% [0, actual * 1/5], repmat([90, 90, (-45), (-45), 0, 45, 45, 90, (-45), (-45), 0, 45, 45, 90, (-45), (-45), 0, 45, 45, 90, (-45)], 1, 1); % Section 1
% [actual * 1/5, actual * 4/5], repmat([90, (-45), (-45), 0, 45, 45, 90, (-45), (-45), 0, 45, 45, 90, (-45), (-45), 0, 45, 45, 90, (-45)], 1, 1); % Section 2
% [actual * 4/5, actual], repmat([(-45), (-45), 0, 45, 45, 90, (-45), (-45), 0, 45, 45, 90, (-45), (-45), 0, 45, 45, 90, (-45)], 1, 1)}; % Section 3
%rear spar web 

% [0, actual * 1/5], repmat([90, 90, 45, 0, 0, -45, 90, 45, repmat([0, 0, 45, 0, 0,-45], 1, 5), 45, 90, 90, -45], 1, 2); % Section 1
% [actual * 1/5, actual * 4/5], repmat([90, 45, 0, 0, -45, 90, 45, repmat([0, 0, 45, 0, 0, -45], 1, 5), 45, 90, 90, -45], 1, 2); % Section 2
% [actual * 4/5, actual], repmat([45, 0, 0, -45, 90, 45, repmat([0, 0, 45, 0, 0,-45], 1, 5), 45, 90, 90, -45], 1, 2)}; % Section 3
%rear spar cap


% [0, actual * 1/5], [90, 90, 45, 0, 0, -45, 90, 45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 45, 90, -45]; % Section 1
% [actual * 1/5, actual * 4/5], [90, 45, 0, 0, -45, 90, 45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 45, 90, -45]; % Section 2
% [actual * 4/5, actual], [45, 0, 0, -45, 90, 45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 45, 90, -45]; % Section 3
% front spar cap

% [0, actual * 1/5], [90, 90, -45, -45, 90, 45, 45, 0, -45, -45, 90, 45, 45, 0, -45, -45, 90, -45, -45]; % Section 1
% [actual * 1/5, actual * 4/5], [90, 90, -45, -45, 90, 45, 45, 0, -45, -45, 90, 45, 45, 0, -45, -45, 90, -45, -45]; % Section 2
% [actual * 4/5, actual], [90, 90, -45, -45, 90, 45, 45, 0, -45, -45, 90, 45, 45, 0, -45, -45, 90, -45, -45]; % Section 3
% front spar web 

% ^ horizontal

% [0, actual * 1/5], [90, 90, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, 45]; % Section 1
% [actual * 1/5, actual * 4/5], [90, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, 45]; % Section 2
% [actual * 4/5, actual], [-45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, 45]; % Section 3
% rear spar cap

% [0, actual * 1/5], [90, 90, -45, -45, 0, 45, 45, 90, -45, -45, 0, 45, 45, 90, -45, -45, 0, 45, 45, 90, -45]; % Section 1
% [actual * 1/5, actual * 4/5], [90, -45, -45, 0, 45, 45, 90, -45, -45, 0, 45, 45, 90, -45, -45, 0, 45, 45, 90, -45]; % Section 2
% [actual * 4/5, actual], [-45, -45, 0, 45, 45, 90, -45, -45, 0, 45, 45, 90, -45, -45, 0, 45, 45, 90, -45]; % Section 3
% rear spar web

% [0, actual * 1/5], [90, 90, 45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 90, 90, -45]; % Section 1
% [actual * 1/5, actual * 4/5], [90, 45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 90, 90, -45]; % Section 2
% [actual * 4/5, actual], [45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 90, 90, -45]}; % Section 3
% front spar cap

% [0, actual * 1/5], [90, 90, -45, -45, 0, 45, 45, 90, -45, -45, 0, 45, 45, 90, -45, -45, 0, 45, 45, 90, -45]; % Section 1
% [actual * 1/5, actual * 4/5], [90, -45, -45, 0, 45, 45, 90, -45, -45, 0, 45, 45, 90, -45, -45, 0, 45, 45, 90, -45]; % Section 2
% [actual * 4/5, actual], [-45, -45, 0, 45, 45, 90, -45, -45, 0, 45, 45, 90, -45, -45, 0, 45, 45, 90, -45]}; % Section 3
% front spar web

%% -------------------- USER INPUTS -------------------- %%
% Define spanwise locations
spanwise_stations = linspace(0, 10.68 * 2 / 2, 10); % 30m span, 50 spanwise points

actual = 26.41 / 2;
% Thickness exaggeration for better visibility
thickness_exaggeration = 2;  

% Define ply stackings for different regions (ensuring symmetry about midline)
% Format: {Span range, Ply angles}
PlyStackingRegions = {...
[0, actual * 1/5], [90, 90, 45, 0, 0, -45, 90, 45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 45, 90, -45]; % Section 1
[actual * 1/5, actual * 4/5], [90, 45, 0, 0, -45, 90, 45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 45, 90, -45]; % Section 2
[actual * 4/5, actual], [45, 0, 0, -45, 90, 45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 0, 0, 45, 0, 0, -45, 45, 90, -45]}; % Section 3



 

% Global stress state at each spanwise location (dummy values)
sigma_x_span = flip(abs(Nx(:,1))); % Stress reducing spanwise
sigma_y_span = flip(abs(zeros(10)));
tau_xy_span  = flip(abs(qFS(:,1)));

% ---- Define Material Strengths by Orientation ----
Material_0deg  = [1.62e11, 1.689e9, 4.67e9, 1.00e9, 3.00e8];  % 0° Ply Properties
Material_45deg = [4.50e10, 2.20e9, 5.00e9, 1.20e9, 3.50e8];  % 45° and -45° Ply Properties
Material_90deg = [9.6e9, 6.23e7, 4.67e9, 1.30e9, 2.50e8];  % 90° Ply Properties

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
    
    % Assign material properties based on ply orientation
    MaterialProps = zeros(length(theta), 5);
    for i = 1:length(theta)
        if theta(i) == 0
            MaterialProps(i, :) = Material_0deg;
        elseif abs(theta(i)) == 45
            MaterialProps(i, :) = Material_45deg;
        elseif theta(i) == 90
            MaterialProps(i, :) = Material_90deg;
        end
    end
    
    % Compute Tsai-Hill for each ply
    FailureIndex = zeros(length(theta), 1);
    for i = 1:length(theta)
        theta_rad = deg2rad(theta(i));
        c = cos(theta_rad);
        s = sin(theta_rad);
        T = [ c^2, s^2,  2*c*s;
              s^2, c^2, -2*c*s;
             -c*s, c*s, c^2 - s^2 ];
        sigma_global = [sigma_x; sigma_y; tau_xy]; 
        sigma_local = T * sigma_global; 
        sigma_1  = sigma_local(1); 
        sigma_2  = sigma_local(2);
        tau_12   = sigma_local(3);
        sigma_1t = MaterialProps(i,1);
        sigma_1c = MaterialProps(i,2);
        sigma_2t = MaterialProps(i,3);
        sigma_2c = MaterialProps(i,4);
        tau_12s  = MaterialProps(i,5);
        sigma_1_max = sigma_1t * (sigma_1 >= 0) + sigma_1c * (sigma_1 < 0);
        sigma_2_max = sigma_2t * (sigma_2 >= 0) + sigma_2c * (sigma_2 < 0);
        TsaiHillIndex = 1e4 * ((sigma_1 / sigma_1_max)^2 ...
                      - (sigma_1 * sigma_2) / sigma_1_max^2 ...
                      + (sigma_2 / sigma_2_max)^2 ...
                      + (tau_12 / tau_12s))^2;
        FailureIndex(i) = min(TsaiHillIndex, 1.1); % Cap values above 1
    end
    
    total_thickness = length(theta) * 0.25 * thickness_exaggeration;
    Y_stack = linspace(0, total_thickness/2, length(theta) + 1); % Only positive side
    FailureMatrix{j} = FailureIndex;
    Y_PlyStacks{j} = Y_stack;
    Span_PlyStacks{j} = spanwise_stations(j) * ones(size(Y_stack));
end

%% ----------------- PLOTTING FAILURE INDEX ----------------- %%
figure;
set(gca, 'DataAspectRatio', [1 2 1]);
hold on;
cmap = parula(256);
caxis([0, max(cellfun(@max, FailureMatrix))]);

for j = 1:length(spanwise_stations)-1
    FailureIndex = FailureMatrix{j};
    Y_stack = Y_PlyStacks{j};
    X_stack = [spanwise_stations(j), spanwise_stations(j+1)];
    for i = 1:length(FailureIndex)
        if FailureIndex(i) >= 1
            color = [1, 0, 0]; % Red for failure
        else
            colorIdx = round(FailureIndex(i) / max(cellfun(@max, FailureMatrix)) * 255) + 1;
            colorIdx = max(1, min(256, colorIdx));
            color = cmap(colorIdx, :);
        end
        fill([X_stack(1), X_stack(2), X_stack(2), X_stack(1)], ...
             [Y_stack(i), Y_stack(i), Y_stack(i+1), Y_stack(i+1)], color, 'EdgeColor', 'k');
    end
end

% LaTeX formatting and boldened axes
xlabel('\textbf{Span Location (m)}', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('\textbf{Thickness (mm)}', 'Interpreter', 'latex', 'FontSize', 14);
% title('\textbf{Tsai-Hill Failure Index Across Span and Thickness, Rear Caps (One Half)}', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontWeight', 'bold', 'LineWidth', 1.5);
colormap(parula);
colorbar;
grid on;
xlim([0, max(spanwise_stations)]);
% exportgraphics(gcf, 'frontsparwebverticaltsai.png', 'Resolution', 400)
hold off;
