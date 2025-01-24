function [V_S] = Stall_Speed(CL, h, W, Sref, Sweep)
% This function is used to calculate the stall speed due to the
% compressibility effects use iterative method
%
% [V_S] = Stall_Speed(CL, V, h, Sweep)
%
% outputs: V_S     - equivalent stall speed (m/s)
% inputs:  CL      - lift coefficient without compressibility effect
%          h       - altitude
%          W       - weight (N)
%          Sred    - Reference wing area (m^2)
%          Sweep   - wing 1/4 chord sweep (rad)

% Basic parameters

[T0, a0, P0, rho0] = atmosisa(0);

[T, a, P, rho] = atmosisa(h);

% Iterative method

error = 1;

V_S = 10;

while abs(error) > 1e-5

CL_C = Comp_Effects(CL, V_S, h, Sweep);

V_S_new = sqrt(2 .* W ./ (Sref .* rho0 .* CL_C));

error = V_S_new - V_S;

V_S = V_S + 0.5 * error;

end

end