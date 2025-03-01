load = linspace(0.1,5,20); % Load values
t = 1; % Avoid division by zero
sigma_c = 20; 
b = linspace(0.1,10,20); % Width values
h = 3; % Fixed height

% Compute initial number of covers
n0covers = ceil(load ./ (t * sigma_c));

% Introduce jagged behavior by using ceil and mod
n0coversactual = n0covers + ceil(0.125 * (n0covers)) * 2 + 2;

% Calculate mass
mass = n0coversactual .* b * h;

% Plot mass to observe jagged pattern
figure;
plot(b, mass, 'b', 'LineWidth', 1.5);
xlabel('Load');
ylabel('Mass (n0coversactual * b * h)');
title('Jagged Mass Distribution');
grid on;
