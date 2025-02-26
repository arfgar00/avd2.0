function [I_x, I_y, A] = I_f(t_f, b_f, h_f)
% This function calculates the second moment of inertia (I_x and I_y) 
% and the cross-sectional area (A) for a C-shaped frame.
%
% Outputs: I_x    - second moment of area about x-axis (m^4)
%          I_y    - second moment of area about y-axis (m^4)
%          A      - frame cross-sectional area (m^2)
% Inputs:  t_f    - frame thickness (m)
%          b_f    - frame width (m)
%          h_f    - Frame height (m)

% Cross-sectional area
A = (t_f * h_f) + 2 * (b_f * t_f);

% Compute y_bar (vertical centroid, measured from the bottom)
y_bar = 0 / A;

% Compute x_bar (horizontal centroid, measured from the left edge)
x_bar = (t_f * h_f * (b_f / 2)) / A;

% Second moment of area about x-axis (Ix)
I_x_web = (1/12) * t_f * h_f^3 + t_f * h_f * (y_bar)^2;
I_x_flange = (1/12) * b_f * t_f^3 + b_f * t_f * (y_bar - h_f/2)^2;
I_x = I_x_web + 2 * I_x_flange;

% Second moment of area about y-axis (Iy)
I_y_web = (1/12) * h_f * t_f^3 + h_f * t_f * (x_bar - b_f/2)^2;
I_y_flange = (1/12) * b_f^3 * t_f + b_f * t_f * (x_bar - 0)^2;
I_y = I_y_web + 2 * I_y_flange;

end