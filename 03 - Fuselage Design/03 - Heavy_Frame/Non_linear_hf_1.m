function [c, ceq] = Non_linear_hf_1(x)
% This function is used to define the non-linear constraint values
%
% output:  c      - inequality constriant
%          ceq    - equality constraint
% inputs:  x      - (L_f, t_f, b_f, h_f)
%          t_f    - frame thickness (m)
%          b_f    - frame width (m)
%          h_f    - frame height (m)

% Recall the maximum loads
S_max = 1.040490021754786e+06;
N_max = 9.035042617803515e+05;
M_max = 1.584783676014650e+06;

% Assign x values
t_web = x(1); % web thickness (m)
t_fl = x(2);  % flange thickness (m)
b_f = x(3);   % frame width (m)
h_f = x(4);   % frame height (m)

% Recall the initial design parameters inside the function
E_f = 67000000000;                 % Young's Modulus (N/m^2)
sigma_y = 614000000;               % Tensile yield stress (Pa)
tau_y = sigma_y ./ sqrt(3);        % Shear yield stress (Pa)

% Find first moment of inertia
Q_x = b_f * t_fl * (h_f/2 - t_fl/2);

% Find second moment of inertia
[I_x, A] = I_hf(t_web, t_fl, b_f, h_f);

% Define the inequality constriant
%
% Slenderness ratio
c(1) = (h_f / t_web) - 1.2 * sqrt(E_f/sigma_y);
c(2) = (b_f / t_fl) - 0.7 * sqrt(E_f/sigma_y);
%
% Maximum direct stress
c(3) = N_max / A - sigma_y;
%
% Maximum shear stress
c(4) = (S_max * Q_x) / (I_x * t_web) - tau_y;
% c(4) = S_max / A - tau_y;
%
% Slenderness ratio
c(5) = 3 - (h_f / t_web);
c(6) = 8 - (b_f / t_fl);
c(7) = 1.5 - (t_fl / t_web);
c(8) = 0.3 - (b_f / h_f);
c(9) = (b_f / h_f) - 0.8; 
% Maximum bending stress
% c(10) = (M_max * h_f) / (2 * I_x) - sigma_y;

% Define the equality constraint
ceq = [];
% Maximum  bending stress
ceq(1) = (M_max * h_f) / (2 * I_x) - sigma_y;

end