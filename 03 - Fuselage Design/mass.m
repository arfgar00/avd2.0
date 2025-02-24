function m_l = mass(t_skin, s_st, t_st, L_st)
% This function is used to calculate the skin-stringer panel mass per
% unit length, according to the the 'Fuselage Design' Excel sheet
%
% output:  m_l    - skin-stringer panel mass per unit length (kg/m)
% inputs:  t_skin - skin thickness (m)
%          s_st   - stringer spacing (m)
%          t_st   - stringer thickness (m)
%          L_st   - stringer flange (m)

% Recall the initial design parameters inside the function
D = 6.38;
r = D / 2;
C = 2 * pi * r;
rho_skin = 2780;
rho_st = 2850;
h_st = L_st / 0.3;

% Calcualte Z-shaped stringer area
B_st = 2 * (L_st * t_st) + (h_st - 2 * t_st) * t_st;

% Calcualte total Z-shaped stringer mass per unit length (kg/m)
m_st = rho_st * B_st * (C / s_st);

% Calcualte skin mass per unit length (kg/m)
m_skin = rho_skin * (C * t_skin);

% Mass per unit length (kg/m)
% m_l = m_st + m_skin;

% Apply penalty if the number of stringers is excessive
min_spacing = 0.16;
penalty_factor = 500^2;

if s_st < min_spacing
    spacing_penalty = penalty_factor * (min_spacing - s_st)^2;
else
    spacing_penalty = 0;
end

% Define cost scaling factor
k = 0.01;  % Scaling coefficient for non-standard flange-width ratios

% Compute cost factor
cost_factor = 1 + k * ((L_st / t_st) - 10)^2;

% Apply penalty if ratio is outside the 8 to 12 range
if (L_st / t_st) > 12 || (L_st / t_st) < 8
    cost_factor = cost_factor * 1.2;  % Apply mild penalty if out of range
end

% Apply cost factor to mass function
m_l = cost_factor * (m_st + m_skin) + spacing_penalty;

end