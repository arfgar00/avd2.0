function [I_1, I_2, A] = I_str(t_st, L_st)
% This function is used to calculate the second-moment of inertia of
% Z-shaped stringer
%
% output:  I_1    - principal axes Ixx (m^4)
%          I_2    - principal axes Iyy (m^4)
%          A      - stringer cross-sectional area (m^4)
% inputs:  t_st   - stringer thickness (m)
%          L_st   - stringer flange (m)

h_st = L_st / 0.3; % stringer height (m)

% Compute Sectional Areas
A1 = L_st * t_st; % Top flange
A2 = (h_st - 2 * t_st) * t_st; % Web
A3 = L_st * t_st; % Bottom flange

% Compute Centroids
y1 = h_st - t_st / 2; % Top flange
y2 = (h_st - 2 * t_st) / 2; % Web
y3 = t_st / 2; % Bottom flange

x1 = L_st / 2; % Top flange (center)
x2 = 0; % Web (on Y-axis)
x3 = -L_st / 2; % Bottom flange

% Compute Moments of Inertia (about own centroid)
Ix1 = (1/12) * L_st * t_st^3; % Top flange
Ix2 = (1/12) * (h_st - 2*t_st) * t_st^3; % Web
Ix3 = (1/12) * L_st * t_st^3; % Bottom flange

Iy1 = (1/12) * t_st * L_st^3; % Top flange
Iy2 = (1/12) * t_st * (h_st - 2*t_st)^3; % Web
Iy3 = (1/12) * t_st * L_st^3; % Bottom flange

% Compute Parallel Axis Theorem Contributions
% Ix = Ix1 + A1 * y1^2 + Ix2 + A2 * y2^2 + Ix3 + A3 * y3^2;
% Iy = Iy1 + A1 * x1^2 + Iy2 + A2 * x2^2 + Iy3 + A3 * x3^2;
Ix = Ix1 + A1 * y1^2 + Ix2 + A2 * y2^2 + Ix3 + A3 * y3^2;
Iy = Iy1 + A1 * x1^2 + Iy2 + A2 * x2^2 + Iy3 + A3 * x3^2;
Ixy = A1 * x1 * y1 + A2 * x2 * y2 + A3 * x3 * y3;

% Compute Principal Moments of Inertia
I_1 = (Ix + Iy) / 2 + sqrt(((Ix - Iy) / 2)^2 + Ixy^2);
I_2 = (Ix + Iy) / 2 - sqrt(((Ix - Iy) / 2)^2 + Ixy^2);
A = A1 + A2 + A3;

end