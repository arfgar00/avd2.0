function [V_C] = C2E(V, h)
% This function is used to calculate the calibrated airspeed
% corrected for compressibility
%
% [V_C] = C2E(V, h)
%
% outputs: V_C     - calibrated airspeed
% inputs:  V       - equivalent airspeed (m/s)
%          h       - altitude

% Calculate True Airspeed

[T0, a0, P0, rho0] = atmosisa(0);

[T, a, P, rho] = atmosisa(h);

V_TAS = V .* sqrt(rho0 ./ rho);

% Calculate Mach number

M = V_TAS / a;

% Calcualte the calibrated airspeed

V_C = V * sqrt((1 + 0.2 .* M .^ 2) ./ (1 - M .^ 2));

end