load =  1562500;
t = linspace(0,2e-4,10);  
t(1) = 1e-6;  % Avoid division by zero
sigmac = 1e9;

% Compute n0covers (element-wise division)
n0covers = ceil(load ./ (t * sigmac));  % Use ./ for element-wise division
n0covers = n0covers + 2 * ceil(0.125 * n0covers / 2) + 2;

% Define linspace for plotting
b = linspace(0,10,10);

% Define missing variables (example values)
h = 0.1;
c = 0.4;       
% Compute mass
mass = b .* t;

% Plot the results
figure;
plot( n0covers,mass, '-o');  % Directly plot n0covers (it's a vector)
xlabel('NumberofCovers');
ylabel('Mass');
title('Plot of n0covers vs. Mass');
grid on;

