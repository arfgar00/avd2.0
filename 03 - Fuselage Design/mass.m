function m_l = mass(t_skin, s_st, t_st, L_st, h_st)
% This function is used to calculate the skin-stringer panel mass per
% unit length, according to the the 'Fuselage Design' Excel sheet
%
% output:  m_l    - skin-stringer panel mass per unit length (kg/m)
% inputs:  t_skin - skin thickness (m)
%          s_st   - stringer spacing (m)
%          t_st   - stringer thickness (m)
%          L_st   - stringer flange (m)
%          h_st   - stringer height (m)

% Recall the initial design parameters inside the function
D = 6.38;
r = D / 2;
C = 2 * pi * r;
rho_skin = 2780;
rho_st = 2850;

% Calcualte Z-shaped stringer area
B_st = 2 * (L_st * t_st) + (h_st - 2 * t_st) * t_st;

% Calcualte total Z-shaped stringer mass per unit length (kg/m)
m_st = rho_st * B_st * (C / s_st);

% Calcualte skin mass per unit length (kg/m)
m_skin = rho_skin * (C * t_skin); 

% Mass per unit length (kg/m)
m_l = m_st + m_skin;

end