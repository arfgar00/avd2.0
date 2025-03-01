% NACA 63-015A Airfoil with Wing Box - Scalable Plot

clear; clc; close all;

%% User Input for Chord Length
% chord_length = input('Enter desired chord length: ');
chord_length = 1

%% Define Original Airfoil Data (Normalized to Chord = 1)
x_original = [0.000000 0.005000 0.007500 0.012500 0.025000 0.050000 0.075000 ...
              0.100000 0.150000 0.200000 0.250000 0.300000 0.350000 0.400000 ...
              0.450000 0.500000 0.550000 0.600000 0.650000 0.700000 0.750000 ...
              0.800000 0.850000 0.900000 0.950000 1.000000];

y_upper =   [0.000000 0.012030 0.014480 0.018440 0.025790 0.036180 0.043820 ...
              0.049970 0.059420 0.066190 0.070910 0.073840 0.074960 0.074350 ...
              0.072150 0.068580 0.063870 0.058200 0.051730 0.044680 0.037310 ...
              0.029910 0.022520 0.015120 0.007720 0.000320];

y_lower =   [0.000000 -0.012030 -0.014480 -0.018440 -0.025790 -0.036180 -0.043820 ...
             -0.049970 -0.059420 -0.066190 -0.070910 -0.073840 -0.074960 -0.074350 ...
             -0.072150 -0.068580 -0.063870 -0.058200 -0.051730 -0.044680 -0.037310 ...
             -0.029910 -0.022520 -0.015120 -0.007720 -0.000320];

%% Scale the Airfoil
x_scaled = x_original * chord_length;
y_upper_scaled = y_upper * chord_length;
y_lower_scaled = y_lower * chord_length;

%% Define Wing Box Dimensions (as percentage of chord)
% chord for horizontal == 4.24 to 2.12
% chord for vertical == 7.63 to 3.81

box_x_front = 0.12 * chord_length;   % Front spar at 12% chord
box_x_rear = 0.60 * chord_length;    % Rear spar at 60% chord
box_y_top = 0.05 * chord_length;    % Top plate (8% of chord)
box_y_bottom = -0.05 * chord_length; % Bottom plate (-8% of chord)

%% Compute Important Points
aero_center_x = 0.25 * chord_length;   % Aerodynamic center at 25% chord
aero_center_y = 0;                     % Typically at zero lift line

shear_center_x = (box_x_front + box_x_rear) / 2; % Middle of wing box
shear_center_y = (box_y_top + box_y_bottom) / 2; % Middle of wing box height

%% Plot Airfoil
figure;
hold on;
plot(x_scaled, y_upper_scaled, 'b', 'LineWidth', 2); % Upper Surface
plot(x_scaled, y_lower_scaled, 'r', 'LineWidth', 2); % Lower Surface

%% Plot Wing Box
% Bottom plate
plot([box_x_front, box_x_rear], [box_y_bottom, box_y_bottom], 'k-', 'LineWidth', 2);
% Top plate
plot([box_x_front, box_x_rear], [box_y_top, box_y_top], 'k-', 'LineWidth', 2);
% Front spar4
plot([box_x_front, box_x_front], [box_y_bottom, box_y_top], 'k-', 'LineWidth', 2);
% Rear spar
plot([box_x_rear, box_x_rear], [box_y_bottom, box_y_top], 'k-', 'LineWidth', 2);

%% Plot Aerodynamic and Shear Centers
plot(aero_center_x, aero_center_y, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'Aerodynamic Center (c/4)');
plot(shear_center_x, shear_center_y, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'DisplayName', 'Shear Center');



%% Final Formatting
grid on;
axis equal;
xlabel('\textbf{Chord Length (x)}', 'Interpreter', 'latex');
ylabel('\textbf{Thickness (y)}', 'Interpreter', 'latex');
title(sprintf('\\textbf{NACA 63-015A Airfoil with Wing Box (Chord = %.2f)}', chord_length), 'Interpreter', 'latex');
legend({'\textbf{Upper Surface}', '\textbf{Lower Surface}', '\textbf{Wing Box}'}, ...
       'Location', 'northeast', 'Interpreter', 'latex');

xlim([0 1.1*chord_length]);
disp('Airfoil and wing box plotted successfully.');

% Customize axis appearance
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.LineWidth = 1.25;
ax.Box = 'off';

fig = gcf;  % Get handle of the current figure (instead of creating a new one)
pbaspect([2.5 1 1]); % Sets the aspect ratio of the axes to 2:1

% Save the existing figure
print(fig, 'airfoil_plot.png', '-dpng', '-r600');  % High-res PNG
print(fig, 'airfoil_plot.pdf', '-dpdf', '-bestfit');  % High-quality PDF
