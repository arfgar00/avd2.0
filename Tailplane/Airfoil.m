classdef Airfoil
    properties
        alpha          % angle of attack original, an array
        Cl             % Cl polar original, an array
        Cd             % Cd polar original, an array
        Cm             % Cm polar original, an array
        uppershape     % coordinate of shape, original
        lowershape     % coordinate of shape, original
        upperShapefun  % polyfit of upper surface shape
        lowerShapefun  % polyfit of lower surface shape
        x_c            % linspace(0,1,500), chordwise location
        t_cfun         % thickness/chord function from interpolation
        t_c            % maximum thickness/chord
        x_cm           % x_c location at maximum thickness
        a0             % CL = a0*alpha + b
        b              % CL = a0*alpha + b
        alphaSPos      % positive stall angle of attack
        alphaSNeg      % negative stall angle of attack
        CM0            % CM at zero lift angle of attack
        Clmax          % maximum Cl
    end
    
    methods
        % Read polar data
        function obj = readPolar(obj, filename)
            tab = readtable(filename);
            obj.Cl = tab.Cl;
            obj.alpha = tab.Alpha .* pi / 180;
            obj.Cd = tab.Cd;
            obj.Cm = tab.Cm;
            obj.Clmax = max(obj.Cl);
            % Define the range in radians
            alpha_min_deg = -7;  % Minimum angle in degrees
            alpha_max_deg = 7;   % Maximum angle in degrees
            alpha_min = alpha_min_deg * pi / 180;  % Convert to radians
            alpha_max = alpha_max_deg * pi / 180;  % Convert to radians

            % Find indices within the specified range
            idx_linear = (obj.alpha >= alpha_min) & (obj.alpha <= alpha_max);

            % Extract the linear region data
            alpha_linear = obj.alpha(idx_linear);
            Cl_linear = obj.Cl(idx_linear);

            % Linear regression for Cl
            p = polyfit(alpha_linear, Cl_linear, 1);
            obj.a0 = p(1);  % Lift curve slope (Cl per radian)
            obj.b = p(2);   % Intercept (Cl at alpha = 0)
            [~, maxIndex] = max(obj.Cl);
            obj.alphaSPos = obj.alpha(maxIndex);
            [~, minIndex] = min(obj.Cl);
            obj.alphaSNeg = obj.alpha(minIndex);
            alpha0L = -obj.b/obj.a0;
            obj.CM0 = interp1(obj.alpha,obj.Cm,alpha0L);
        end

        % Plot polar data (Cl, Cd, and Cl/Cd)
        function plotPolar(obj)
            subplot(2,2,1)
            plot(obj.alpha.*180/pi, obj.Cl)
            hold on
            plot(obj.alpha.*180/pi, obj.alpha*obj.a0 + obj.b,'--')
            xlabel("alpha")
            ylabel("Cl")
            subplot(2,2,2)
            plot(obj.alpha.*180/pi, obj.Cd)
            xlabel("alpha")
            ylabel("Cd")
            subplot(2,2,3)
            plot(obj.alpha.*180/pi, obj.Cl./obj.Cd)
            xlabel("alpha")
            ylabel("Cl/Cd")
        end

        % Interpolate Cl and Cd for a given angle of attack
        function [Cl, Cd] = interpPolar(obj, alpha)
            Cl = interp1(obj.alpha, obj.Cl, alpha);
            Cd = interp1(obj.alpha, obj.Cd, alpha);
        end

        % Read airfoil shape data
        function obj = readShape(obj, filename)
            obj.x_c = linspace(0, 1, 500);
            
            fid = fopen(filename, 'r');
            if fid == -1
                error('Cannot open the file: %s', filename);
            end

            headerLine1 = fgetl(fid); % Airfoil name
            headerLine2 = fgetl(fid); % Numbers (e.g., 103. 103.)

            coordinates = [];
            while ~feof(fid)
                line = fgetl(fid);
                if isempty(line)
                    continue;
                end
                data = sscanf(line, '%f %f');
                if length(data) == 2
                    coordinates = [coordinates; data'];
                end
            end
            fclose(fid);

            x_diff = diff(coordinates(:,1));
            idx_reset = find(x_diff < 0, 1, 'first');
            
            if isempty(idx_reset)
                error('Could not find the splitting point between upper and lower surfaces.');
            end

            upperSurface = coordinates(1:idx_reset, :);
            lowerSurface = coordinates(idx_reset+1:end, :);

            if lowerSurface(1,1) > lowerSurface(end,1)
                lowerSurface = flipud(lowerSurface);
            end
            obj.uppershape = upperSurface;
            obj.lowershape = lowerSurface;
        end

        % Interpolate shape using polynomial fitting
        function obj = interpShape(obj, n)
            p_upper = polyfit(obj.uppershape(:,1), obj.uppershape(:,2), n);
            p_lower = polyfit(obj.lowershape(:,1), obj.lowershape(:,2), n);
            obj.upperShapefun = @(x_c) polyval(p_upper, x_c);
            obj.lowerShapefun = @(x_c) polyval(p_lower, x_c);
            
            thickness = obj.upperShapefun(obj.x_c) - obj.lowerShapefun(obj.x_c);
            p_thickness = polyfit(obj.x_c, thickness, n);
            obj.t_cfun = @(x_c) polyval(p_thickness, x_c);
            
            thickness_values = obj.t_cfun(obj.x_c);
            [obj.t_c, max_index] = max(thickness_values);
            obj.x_cm = obj.x_c(max_index);
        end

        % Calculate L/D ratio at a given Cl
        function L_Ddes = L_DatCL(obj, CLdes)
            if CLdes > min(obj.Cl) && CLdes < max(obj.Cl)
                alphades = (CLdes - obj.b)/obj.a0;
                CDdes = interp1(obj.alpha, obj.Cd, alphades);
                L_Ddes = CLdes / CDdes;
            else
                error("Input CL out of range")
            end
        end

        % Displace the shape relative to quarter chord point
        function obj = displaceShape(obj, dx, dy)
            % Quarter chord point (x = 0.25, y at x=0.25 from upper surface)
            quarterChordX = 0.25;
            quarterChordY_upper = obj.upperShapefun(quarterChordX);
            quarterChordY_lower = obj.lowerShapefun(quarterChordX);
            quarterChordY = (quarterChordY_upper + quarterChordY_lower) / 2;

            % Calculate displacement
            dispX = dx - quarterChordX;
            dispY = dy - quarterChordY;

            % Apply displacement to both upper and lower surfaces
            obj.uppershape(:, 1) = obj.uppershape(:, 1) + dispX;
            obj.uppershape(:, 2) = obj.uppershape(:, 2) + dispY;
            obj.lowershape(:, 1) = obj.lowershape(:, 1) + dispX;
            obj.lowershape(:, 2) = obj.lowershape(:, 2) + dispY;
        end

        % Plot the airfoil shape with optional scaling factor and mark quarter chord point
        function plotShape(obj, scaleFactor)
            if nargin < 2
                scaleFactor = 1;  % Default scale factor is 1 (no scaling)
            end

            % Scale the airfoil shapes
            scaledUpper = obj.uppershape * scaleFactor;
            scaledLower = obj.lowershape * scaleFactor;

            % Plot the scaled shapes
            hold on;
            plot(scaledUpper(:,1), scaledUpper(:,2), 'b--', 'LineWidth', 2);
            hold on;
            plot(scaledLower(:,1), scaledLower(:,2), 'r--', 'LineWidth', 2);

            % Mark the quarter-chord point
            quarterChordX = 0.25;
            quarterChordY_upper = obj.upperShapefun(quarterChordX);
            quarterChordY_lower = obj.lowerShapefun(quarterChordX);
            quarterChordY = (quarterChordY_upper + quarterChordY_lower) / 2;
            plot(quarterChordX * scaleFactor, quarterChordY * scaleFactor, 'go', 'MarkerSize', 4, 'MarkerFaceColor', 'g');
            %text(quarterChordX * scaleFactor - 5, quarterChordY * scaleFactor + 5, '  Quarter Chord', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        end
    end
end
