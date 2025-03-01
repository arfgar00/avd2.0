function m_f = mass_f(t_f, b_f, h_f)
% This function is used to calculate the light frames total mass,
% according to the the 'Fuselage Design' Excel sheet
%
% output:  m_f    - total frame mass (kg)
% inputs:  t_f    - frame thickness (m)
%          b_f    - frame width (m)
%          h_f    - frame height (m)

% Recall the initial design parameters inside the function
D = 6.38;                          % Fuselage diameter (m)
r = D ./ 2;                        % Fuselage radius (m)
C = 2 * pi * r;                    % Fuselage circumference (m)
L_fuse = 69.1;                     % Fuselage length (m)

E_f = 71700000000;                 % Young's Modulus (N/m^2)
rho_f = 2850;                      % Density (kg/m^3)

L_f = 0.508;    % frame spacing (m)

% Number of frame
n = ceil((L_fuse - 2 * L_f) / L_f);

% sigma_crit = 6.053057743137613e+07; % Critical buckling stress (N/m^2)

% Find second moment of inertia
[I_x, I_y, A] = I_f(t_f, b_f, h_f);

% Single frame volume
V = A * C;

% Calcualte total frame mass
m_f = n * (rho_f * V);

% m_f = n * (rho_f * V) + (E_f * I_x) / (rho_f * V);

end