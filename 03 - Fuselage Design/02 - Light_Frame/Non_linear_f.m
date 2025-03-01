function [c, ceq] = Non_linear_f(x)
% This function is used to define the non-linear constraint values
%
% output:  c      - inequality constriant
%          ceq    - equality constraint
% inputs:  x      - (L_f, t_f, b_f, h_f)
%          t_f    - frame thickness (m)
%          b_f    - frame width (m)
%          h_f    - frame height (m)

% Assign x values
t_f = x(1);   % frame thickness (m)
b_f = x(2);   % frame width (m)
h_f = x(3);   % frame height (m)

% Recall the initial design parameters inside the function
D = 6.38;                          % Fuselage diameter (m)
r = D ./ 2;                        % Fuselage radius (m)
% C = 2 * pi * r;                  % Fuselage circumference (m)
L_fuse = 69.1;                     % Fuselage length (m)

E_f = 71700000000;                 % Young's Modulus (N/m^2)
% rho_f = 2850;                    % Density (kg/m^3)

M = 5.067912648168767e+06;        % Maximum bending moment (Nm)

sigma_crit = 6.053057743137613e+07; % Critical buckling stress (N/m^2)

L_f = 0.508;    % frame spacing (m)

% Number of frame
n = ceil((L_fuse - 2 * L_f) / L_f);

% Adjust frame spacing
L_f = (L_fuse - 2 * L_f) / n;

% Find second moment of inertia
[I_x, I_y, A] = I_f(t_f, b_f, h_f);

% Define the inequality constriant
% 
c(1) = (h_f / t_f) - pi * sqrt(E_f/683000000);
%
c(2) = (pi ^ 2 * E_f * I_y) / (A * L_f ^ 2) - sigma_crit;

% Define the equality constraint
%
% Minimum sizing requirement
ceq (1) = E_f * I_x - ((1/16000) * M * D ^ 2) / L_f;
%
% Panel instability equal to general instability
ceq (2) = (pi ^ 2 * E_f * I_y) / (A * L_f ^ 2) - sigma_crit;

end