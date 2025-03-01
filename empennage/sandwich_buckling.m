% Given Material and Panel Properties
K_F = 0;          % Buckling coefficient from face bending (assumed negligible)
K_M = 3;          % Buckling coefficient from bending shear (typically 3-6)
K = K_F + K_M;    % Total buckling coefficient

b_initial = 1.5;  % Panel width (m)
D = 445440;       % Flexural rigidity (Nm) (Typical for aerospace panels)

h = linspace(0.01, 0.10, 10);  % Core thickness (m) (Nomex honeycomb)
t_c = 0.002;      % Face sheet thickness (m) (Carbon fiber skin)
G_c = 30e6;       % Core shear modulus (N/m²) (Nomex honeycomb)

a = 2.12 * (0.7 - 0.12); % (m)
b = a * 0.9; % Reassigning panel width based on `a`

% Compute critical load per unit width:
N_cr = K * (pi^2 ./ b^2) * D; % Scalar

% Compute transverse shear stiffness:
U = (h.^2 ./ t_c) * G_c; % Vector

% Compute bending shear rigidity:
V = (pi^2 * D) ./ (b^2 .* U) % Vector

% Display results
fprintf('Critical Load per Unit Width (N_cr): %.2f N/m\n', N_cr);

% Displaying U and V as vectors
fprintf('Transverse Shear Stiffness (U) for each h value:\n');
disp(U);

fprintf('Bending Shear Rigidity (V) for each h value:\n');
disp(V);

% Loop through Nx and check for buckling failure
for i = 1:length(Nx)
    if Nx(i) > N_cr
        fprintf('Element %d exceeds critical load: Nx = %.2f N/m > N_cr = %.2f N/m (Buckling Failure)\n', i, Nx(i), N_cr);
    else
        fprintf('Element %d is within safe limits: Nx = %.2f N/m <= N_cr = %.2f N/m\n', i, Nx(i), N_cr);
    end
end
