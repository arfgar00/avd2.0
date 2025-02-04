classdef wingbox
    properties
        t1 %thickness of skin on D-cell
        t2 %thickness of skin on main cell
        tF %thickness of front spar
        tR %thickness of rear spar
        c_c  %total width of main cell, normalized by chord length
        b1_c  %separation between stringers on main cell, normalized by chord length
        R1_c  %radius of curvature of upper arc, normalized by chord length
        R2_c  %radius of curvature of lower arc, normalized by chord length
        a %separation between ribs
        b2_c % height of front spar, normalized by chord length
        ts %thickness of stringers
        ds %flange width of Z stringer (of one - on the Z)
        hs %height of Z stinger
        As %area of one stringer
        sr_s0 %sigma_r/sigma_0, looked up from catchpole diagram, array
        N %number of panels on main cell
    end

    methods
        function obj = calcb1(obj)
            obj.b1_c = obj.c_c/(obj.N+1);
        end

        function obj = calcStringerArea(obj)
            obj.As = (obj.hs + 2 * obj.ds)*obj.ts;
        end
        
        function x = plotDcell(obj)
            % System of equations for:
            % - Upper arc center (h1, k1)
            % - Lower arc center (h2, k2)
            % - Junction point (x_p, y_p)
            R1 = obj.R1_c;
            R2 = obj.R2_c;
            b2 = obj.b2_c;
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
            plot(x_spar+offset, y_spar-b2/2, 'k-', 'LineWidth', 2); hold on; % Front spar
            plot(-x_upper+offset, y_upper-b2/2, 'k-', 'LineWidth', 2);        % Upper arc (R1)
            plot(-x_lower+offset, y_lower-b2/2, 'k-', 'LineWidth', 2);        % Lower arc (R2)
            x = mean(-x_spar+offset);
        end


        function drawwithAirfoil(obj, airfoil)
            % Plot airfoil shape
            airfoil.plotShape();
            hold on;
            
            % Plot D-cell and get its x-position
            x = obj.plotDcell();
            hold on;
            
            % Plot wing box boundaries
            % Bottom plate
            plot([x, x+obj.c_c], [-obj.b2_c/2, -obj.b2_c/2], 'k-', 'LineWidth', 2) 
            hold on;
            % Front spar
            plot([x, x], [-obj.b2_c/2, obj.b2_c/2], 'k-', 'LineWidth', 2)
            hold on;
            % Rear spar
            plot([x+obj.c_c, x+obj.c_c], [-obj.b2_c/2, obj.b2_c/2], 'k-', 'LineWidth', 2)
            hold on;
            % Top plate
            plot([x, x+obj.c_c], [obj.b2_c/2, obj.b2_c/2], 'k-', 'LineWidth', 2) 
            hold on;
        
            % Draw Z-stringers on top and bottom plates (inside the wing box)
            if obj.N > 0
                % Calculate stringer spacing along chord
                spacing = obj.c_c / (obj.N + 1);
                
                % Define Z-stringer geometry (proportional to wing box height)
                z_size = obj.b2_c / 10;  % Size of Z-stringer (scaled to box height)
                
                % Z-shape coordinates (normalized)
                z_shape = [0, 0.5, 0.5, 1;  % X-coordinates (normalized)
                           0.5, 0.5, -0.5, -0.5]; % Y-coordinates (normalized)
                
                % Scale Z-shape to fit inside the wing box
                z_shape = z_shape .* [spacing / 2; z_size];  % Scale to fit spacing and height
                
                % Draw stringers on both plates
                for i = 1:obj.N
                    % Calculate chordwise position
                    x_pos = x + i * spacing;
                    
                    % Top plate stringer (oriented downward)
                    top_z = z_shape + [x_pos; obj.b2_c/2 - z_size];  % Offset inside the box
                    plot(top_z(1,:), top_z(2,:), 'b-', 'LineWidth', 1.5)
                    
                    % Bottom plate stringer (oriented upward, flipped vertically)
                    bot_z = z_shape .* [1; -1] + [x_pos; -obj.b2_c/2 + z_size];  % Offset inside the box
                    plot(bot_z(1,:), bot_z(2,:), 'b-', 'LineWidth', 1.5)
                end
            end
            
            axis equal;
        end
    end
end