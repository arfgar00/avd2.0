function [I_x, A] = I_hf(t_web, t_fl, b_f, h_f)
% This function calculates the second moment of inertia (I_x and I_y) 
% and the cross-sectional area (A) for a I-shaped frame.
%
% Outputs: I_x    - second moment of area about x-axis (m^4)
%          I_y    - second moment of area about y-axis (m^4)
%          A      - frame cross-sectional area (m^2)
% Inputs:  t_web  - web thickness (m)
%          t_fl   - flange thickness (m)
%          b_f    - frame width (m)
%          h_f    - Frame height (m)

% Cross-sectional area
A = (t_web * h_f) + 2 * (b_f * t_fl);

% Second moment of area about x-axis (Ix)
I_x_web = (1/12) * t_web * h_f^3;
I_x_flange = (1/12) * b_f * t_fl^3 + b_f * t_fl * (h_f/2)^2;
I_x = I_x_web + 2 * I_x_flange;

end