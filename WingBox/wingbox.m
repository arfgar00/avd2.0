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
        F %Farror factor, an array in spanwise
        L %rib spacing, an array in spanwise
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

        function obj = interpsr_s0(obj, digitized_data, wing)
            ts_t = obj.ts / obj.t2;  % Current ts/t ratio
            As_bt = obj.As ./ (obj.b1_c .* wing.cn) ./ obj.t2;  % Current As/bt values
        
            % Extract and sort digitized_data by ts_t_ratio
            ts_t_ratios = [digitized_data.ts_t_ratio];
            [sorted_ts, sort_idx] = sort(ts_t_ratios);
            digitized_data_sorted = digitized_data(sort_idx);
        
            % Create interpolants for each curve
            for k = 1:length(digitized_data_sorted)
                [x_sorted, idx] = sort(digitized_data_sorted(k).x_data);  % Sort x_data
                y_sorted = digitized_data_sorted(k).y_data(idx);          % Align y_data
                % Store interpolant in the struct
                digitized_data_sorted(k).interpolant = griddedInterpolant(x_sorted, y_sorted, 'linear', 'none');
            end
        
            % Find bounding ts/t ratios for interpolation
            idx_lower = find(sorted_ts <= ts_t, 1, 'last');  % Nearest lower ts/t
            idx_upper = find(sorted_ts >= ts_t, 1, 'first'); % Nearest upper ts/t
        
            % Handle edge cases (extrapolation)
            if isempty(idx_lower)
                idx_lower = 1;  % Use first curve if ts_t < min(ts_t_ratio)
                idx_upper = 1;
            elseif isempty(idx_upper)
                idx_lower = length(sorted_ts);  % Use last curve if ts_t > max(ts_t_ratio)
                idx_upper = length(sorted_ts);
            end
        
            % Get interpolated y-values for each As_bt
            for i = 1:length(As_bt)
                x_query = As_bt(i);  % Current As/bt value to query
        
                % Get y-values from bounding ts/t curves
                y_lower = digitized_data_sorted(idx_lower).interpolant(x_query);  % From lower ts/t
                y_upper = digitized_data_sorted(idx_upper).interpolant(x_query);  % From upper ts/t
        
                % Linear interpolation between ts/t curves
                t_lower = sorted_ts(idx_lower);
                t_upper = sorted_ts(idx_upper);
                if t_upper == t_lower
                    obj.sr_s0(i) = y_lower;  % Exact match, no interpolation needed
                else
                    alpha = (ts_t - t_lower) / (t_upper - t_lower);  % Weighting factor
                    obj.sr_s0(i) = y_lower + alpha * (y_upper - y_lower);  % Final interpolation
                end
            end
        end
        
        function obj = interpretF(obj, F_data, wing)
            % F_data = 
            %   1×8 struct array with fields:
            %     F, interpolated value
            %     x_data, query value: As_bt = obj.As ./ (obj.b1_c .*
            %     wing.cn) ./ obj.t2,
            %     y_data, query value: ts_t = obj.ts/obj.t2, a constant
            As_bt = obj.As ./ (obj.b1_c .* wing.cn) ./ obj.t2;
            ts_t = obj.ts/obj.t2;
            obj.F = zeros(size(As_bt));
            
            xq = As_bt;
            yq = ts_t;
            Xall = [];
            Yall = [];
            Fall = [];
            for i = 1:length(F_data)
                % Each contour line is F_data(i).x_data vs F_data(i).y_data 
                % with a *constant* value F_data(i).F on that line.
                x_i = F_data(i).x_data(:);
                y_i = F_data(i).y_data(:);
                Xall = [Xall; x_i];
                Yall = [Yall; y_i];
                Fall = [Fall; repmat(F_data(i).F, size(x_i))];
            end
            Finterp = scatteredInterpolant(Xall, Yall, Fall, 'natural', 'nearest');
            for i = 1:length(xq)
                obj.F(i) = Finterp(xq(i), yq);
            end
        end

    end
end