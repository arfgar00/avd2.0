function [CL_C] = Comp_Effects(CL, V, h, Sweep)
% This function is used to calculate the lift coefficient
% due to the compressibility effect using the Prandtl-Glauert Rule
% for subsonic Mach numbers
%
% [CL_C] = Comp_Effects(CL, V, Sweep)
%
% outputs: CL_C    - lift coefficient
% inputs:  CL      - lift coefficient without compressibility effect
%          V       - equivalent airspeed (m/s)
%          h       - altitude
%          Sweep   - wing 1/4 chord sweep (rad)

% Calculate True Airspeed

[T0, a0, P0, rho0] = atmosisa(0);

[T, a, P, rho] = atmosisa(h);

V_TAS = V .* sqrt(rho0 ./ rho);

% Calculate Mach number

M = V_TAS / a;

% Calcualte the lift coefficient

M_eff = M * cos(Sweep);

if M_eff <= 0.3
    
    CL_C = CL;

elseif M_eff <= 0.7
    
    CL_C = CL ./ sqrt(1 - M_eff .^ 2);

else
    
    CL_C = (CL ./ sqrt(1 - M_eff .^ 2)) ...
        * (1 + (M_eff ./ (1 + sqrt(1 - M_eff .^ 2))));

end

end