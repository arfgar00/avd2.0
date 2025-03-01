function plotDcell(R1, R2, b2)
    
    % System of equations for:
    % - Upper arc center (h1, k1)
    % - Lower arc center (h2, k2)
    % - Junction point (x_p, y_p)
    fun = @(vars) [
        (0 - vars(1))^2 + (b2 - vars(2))^2 - R1^2;        % Upper arc passes through (0, b2)
        (vars(5) - vars(1))^2 + (vars(6) - vars(2))^2 - R1^2; % Upper arc passes through (x_p, y_p)
        (0 - vars(3))^2 + (0 - vars(4))^2 - R2^2;          % Lower arc passes through (0, 0)
        (vars(5) - vars(3))^2 + (vars(6) - vars(4))^2 - R2^2; % Lower arc passes through (x_p, y_p)
        (vars(6) - vars(2)) * (vars(5) - vars(3)) - (vars(6) - vars(4)) * (vars(5) - vars(1)) % Tangency
    ];
    
    % Initial guess: Centers to the LEFT of the spar (h1 < 0, h2 < 0)
    initial_guess = [-R1, b2 + R1, -R2, -R2, 0.5, b2/2]; 
    
    % Solve numerically
    options = optimoptions('fsolve', 'Display', 'off', 'Algorithm', 'levenberg-marquardt');
    solution = fsolve(fun, initial_guess, options);
    
    % Extract variables
    h1 = solution(1);   % Upper arc center x
    k1 = solution(2);   % Upper arc center y
    h2 = solution(3);   % Lower arc center x
    k2 = solution(4);   % Lower arc center y
    x_p = solution(5);  % Junction x-coordinate
    y_p = solution(6);  % Junction y-coordinate
    
    % Implicitly compute b (horizontal protrusion)
    b = x_p;
    
    % Generate upper arc (R1: curves downward and outward)
    theta1_start = atan2(b2 - k1, 0 - h1);   % Angle at (0, b2)
    theta1_end = atan2(y_p - k1, x_p - h1);  % Angle at (x_p, y_p)
    theta1 = linspace(theta1_start, theta1_end, 100);
    x_upper = h1 + R1 * cos(theta1);
    y_upper = k1 + R1 * sin(theta1);
    
    % Generate lower arc (R2: curves upward and outward)
    theta2_start = atan2(0 - k2, 0 - h2);    % Angle at (0, 0)
    theta2_end = atan2(y_p - k2, x_p - h2);  % Angle at (x_p, y_p)
    theta2 = linspace(theta2_start, theta2_end, 100);
    x_lower = h2 + R2 * cos(theta2);
    y_lower = k2 + R2 * sin(theta2);
    
    % Vertical front spar
    x_spar = zeros(1, 50);
    y_spar = linspace(0, b2, 50);
    offset = max(max(x_upper) , max(x_lower));
    % Plot
    figure;
    plot(x_spar+offset, y_spar, 'k-', 'LineWidth', 2); hold on; % Front spar
    plot(-x_upper+offset, y_upper, 'b-', 'LineWidth', 2);        % Upper arc (R1)
    plot(-x_lower+offset, y_lower, 'r-', 'LineWidth', 2);        % Lower arc (R2)
    
end
