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

% Recall the initial design parameters inside the function
D = 6.38;
r = D / 2;
C = 2 * pi * r;
h_st = L_st / 0.3;
E = 71300000000;
E_st = 71700000000;
q_st = s_st;
% q_st = s_st - L_st;
sigma_skin = 455000000;

% From ESDU datasheet 
% (see Structure 3 Lecture Notes for k to K calculation)
Kc = 3.69; % 3.62 for v = 0.3
Ks = 4.93; % 4.83 for v = 0.3

% Recall the initial design parameters inside the function
% Q = 5.951582721709430e+05;  % Maximum shear force (N)
% P = 0;                      % Maximum tangential load (N)
M = -5.067912648168767e+06;   % Maximum bending moment (Nm)
% T = 1.263599895277778e+06;  % Maximum torque (Nm)

% Recall the maximum shear flow
q = 7.914986682604700e+04;    % Maximum shear flow (Nm)

% Find the number of the stringer
n = ceil(C / s_st);           % Round-up for safety
s_st = C / n;

% Compute the stringer section area
B_st = 2 * (L_st * t_st) + (h_st - 2 * t_st) * t_st;

% Find the y-location for each case
y = zeros(1, n);

for i = 1 : n

    y(i) = r * sin((2 * pi * (i - 1)) / n);

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

B_i(1) = B_i(n);

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
sigma = zeros(1, n);

for i = 1 : n

    sigma(i) = (M / I_x) * y(i);

end

% Find the maximum tensile stress
sigma_ten = abs(min(sigma));

% Find the maximum bending & compression stress
sigma_ben = abs(max(sigma));

% Find the maximum shear stress
shear = q / t_skin;

% Critical buckling stress of the skin
shear_crit = Ks * E * ((t_skin / q_st) ^ 2);
sigma_crit = Kc * E * ((t_skin / q_st) ^ 2);

% Critical buckling stress of the stringer
sigma_crit_st = Kc * E_st * ((t_st / h_st) ^ 2);

% Find Farrar's efficiency
As_bt = B_st / (s_st * t_skin);
ts_t = t_st / t_skin;
Farrar = interpF(As_bt, ts_t);

% Define the inequality constriant
%
% tensile stress
c(1) = sigma_ten - sigma_skin;

% compression-shear interaction
c(2) = (shear / shear_crit) ^ 2 + sigma_ben / sigma_crit - 0.9999;
% c(3) = 0.985 - ((shear / shear_crit) ^ 2 + sigma_ben / sigma_crit);

% Critical buckling stress of the stringer
c(3) = sigma_ben - sigma_crit_st;

% Stringer design
c(4) = 2 * t_st - h_st;

% Farrar's efficiency
c(5) = As_bt - 1.99;
c(6) = ts_t - 1.79;
c(7) = 0.65 - Farrar;
% c(7) = 0.75 - Farrar;

% Define the equality constraint
%
% stringer must be uniformly spacing
% ceq = mod(C, s_st) / C;
ceq = [];

end