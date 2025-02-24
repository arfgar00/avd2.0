function [c, ceq] = Non_linear(x)
% This function is used to define the non-linear constraint values
%
% output:  c      - inequality constriant
%          ceq    - equality constraint
% inputs:  x      - (t_skin, s_st, t_st, L_st, h_st)
%          t_skin - skin thickness (m)
%          s_st   - stringer spacing (m)
%          t_st   - stringer thickness (m)
%          L_st   - stringer flange (m)
%          h_st   - stringer height (m)

% Assign x values
t_skin = x(1); % skin thickness (m)
s_st = x(2);   % stringer spacing (m)
t_st = x(3);   % stringer thickness (m)
L_st = x(4);   % stringer flange (m)
h_st = x(5);   % stringer height (m)

% Recall the initial design parameters inside the function
D = 6.38;
r = D / 2;
C = 2 * pi * r;
E = 71300000000;
q = s_st - 15 * t_skin;

% From ESDU datasheet 
% (see Structure 3 Lecture Notes for k to K calculation)
Kc = 3.69; % 3.62 for v = 0.3
Ks = 4.93; % 4.83 for v = 0.3

% Calcualte critical buckling stress
shear_crit = Ks * E * (t_skin / q) ^ 2;
sigma_crit = Kc * E * (t_skin / q) ^ 2;

% Recall the initial design parameters inside the function
% Q = 5.951582721709430e+05;  % Maximum shear force (N)
% P = 0;                      % Maximum tangential load (N)
M = -5.067912648168767e+06;   % Maximum bending moment (Nm)
% T = 1.263599895277778e+06;  % Maximum torque (Nm)

% Recall the maximum shear flow
q = 7.914986682604700e+04;    % Maximum shear flow (Nm)

% Find the number of the stringer
n = ceil(C / s_st);           % Round-up for safety
s_st = C / n;                 % real stringer spacing (m)

% Compute the stringer section area
B_st = 2 * (L_st * t_st) + (h_st - 2 * t_st) * t_st;

% Find the y-location for each case
y = zeros(1, n);

for i = 1 : n

    y(i) = r * cos(2 * pi * (i - 1) / n); % Evenly distributed

end

% Compute the skin collaborative area
B_i = zeros(1, n);

B_i(1) = (t_skin * s_st / 6) * (2 + y(end) / y(1)) + ...
    (t_skin * s_st / 6) * (2 + y(2) / y(1)); % First stringer
B_i(n) = (t_skin * s_st / 6) * (2 + y(n-1) / y(n)) + ...
    (t_skin * s_st / 6) * (2 + y(1) / y(n)); % Last stringer

for i = 2 : n - 1

    B_i(i) = (t_skin * s_st / 6) * (2 + y(i-1) / y(i)) + ... 
        (t_skin * s_st / 6) * (2 + y(i+1) / y(i));

end

% Compute total boom area for each stringer
B_total = zeros(1, n);

for i = 1 : n

    B_total(i) = B_st + B_i(i);

end

% Compute the second moment of area, Ix
I_x = 0;

for i = 1 : n

    I_x = I_x + B_total(i) * (y(i) ^ 2);
end

% Find the shear stress distribution
sigma = zeros(1; n)

for i = 1 : n

    sigma = (M / I_x) * y(i);

end

% Find the maximum tensile stress
sigma_ten = abs(min(sigma));

% Find the maximum bending & compression stress
sigma_ben = abs(max(sigma));

% Find the maximum shear stress
shear = q / t_skin;

% Define the inequality constriant
c(1) = sigma_ten - 

% Define the equality constraint
ceq = [];

end