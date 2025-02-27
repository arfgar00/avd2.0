function m_f = mass_hf(t_web, t_fl, b_f, h_f)
% This function is used to calculate the heavy frame mass
%
% output:  m_hf   - frame mass (kg)
% inputs:  x      - (L_f, t_f, b_f, h_f)
%          t_f    - frame thickness (m)
%          b_f    - frame width (m)
%          h_f    - frame height (m)

% Recall the initial design parameters inside the function
D = 6.38;                          % Fuselage diameter (m)
r = D ./ 2;                        % Fuselage radius (m)
C = 2 * pi * r;                    % Fuselage circumference (m)
% E_f = 67000000000;                 % Young's Modulus (N/m^2)
rho_f = 2860;                      % Density (kg/m^3)

% Find second moment of inertia
[I_x, A] = I_hf(t_web, t_fl, b_f, h_f);

% Calcualte frame mass
m_f = rho_f * A * C;

end