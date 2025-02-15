classdef wingbox
    properties
        t_D %thickness of skin on D-cell
        t_Upper %thickness of upper panel
        t_Lower %thickness of lower panel
        tF %thickness of front spar
        tR %thickness of rear spar
        c_c  %total width of main cell, normalized by chord length
        b1_c  %separation between stringers on main cell, normalized by chord length
        b_c   %length of D Cell
        R_c  %radius of curvature of D Cell
        a %separation between ribs
        b2_c % height of front spar, normalized by chord length
        ts %thickness of stringers
        ds %flange width of Z stringer (of one - on the Z)
        hs %height of Z stinger
        As %area of one stringer
        sr_s0 %sigma_r/sigma_0, looked up from catchpole diagram, array
        N %number of stringers
        F %Farror factor, an array in spanwise
        

        t1_spar %spar flange thickness
        t2_spar %spar web thickness
        b_spar  %spar flange width

        An   %area of the D cell
        Ar   %area of the main cell

        Sn   %circumference of D cell
        Sr   %circumference of main cell

        L %rib spacing, an array in spanwise, same length as wing.cn
        Nrib % number of ribs
        yrib % y location of each rib
    end

    methods
        function obj = calcb1(obj,D,wing)
            %1. find wingbox width
            theta = wing.Lambdax_c(0,wing.b/2);
            fs = obj.b_c; %front spar percentage
            rs = obj.b_c + obj.c_c; %rear spar percentage
            R = D/2;
            froot = wing.LEx(R) + wing.c_at_y(R)*fs;
            rroot = wing.LEx(R) + wing.c_at_y(R)*rs;
            ftip = wing.LEx(wing.b/2) + wing.c_at_y(wing.b/2)*fs;
            rtip = wing.LEx(wing.b/2) + wing.c_at_y(wing.b/2)*rs;
            
            xf = @(y) froot + (ftip - froot)/(wing.b/2 - R)*(abs(y)-R);
            xr = @(y) rroot + (rtip - rroot)/(wing.b/2 - R)*(abs(y)-R);
            for i = 1:length(wing.stripy)
                y = wing.stripy(i);
                obj.b1_c(i) = (xr(y) - xf(y))/(obj.N(i)+1)/wing.c_at_y(y);
            end
            
        end

        function obj = calcR(obj)
            obj.R_c = (obj.b2_c^2/4  + obj.b_c^2)/(obj.b2_c);
        end

        function obj = calcStringerArea(obj)
            obj.As = (obj.hs + 2 * obj.ds)*obj.ts;
        end

        function obj = calcCellArea(obj,wing)
            obj.Ar = obj.b2_c * obj.c_c .*wing.cn.^2;
            obj.Sr = wing.cn.*(obj.b2_c + 2*obj.c_c);
            if obj.R_c == []
                obj = obj.calcR();
            else
                theta = atan(obj.b_c/(obj.R_c - obj.b2_c/2));
                obj.An = wing.cn.^2.*(theta*obj.R_c^2 - obj.b_c*(obj.R_c - obj.b2_c/2));
                obj.Sr = 2*theta*obj.R_c.*wing.cn;
            end
        end
        
        function b = plotDcell(obj)
            R = obj.R_c;
            b2 = obj.b2_c;
            b = obj.b_c;
            x_upper = linspace(0,b,30);
            y_upper = sqrt(R^2 - (x_upper - b).^2) + b2/2-R;

            y_lower = -y_upper;
            % Plot
            plot([b, b],[b2/2, -b2/2], 'k-', 'LineWidth', 2); hold on; % Front spar
            hold on
            plot(x_upper,y_upper, 'k-', 'LineWidth', 2)
            hold on
            plot(x_upper,y_lower, 'k-', 'LineWidth', 2)
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

        function obj = interpsr_s0(obj, sr_s0, wing)
            As_bt = obj.As ./ (obj.b1_c .* wing.cn) ./ obj.t_Upper;
            ts_t = obj.ts / obj.t_Upper;
            obj.sr_s0 = zeros(size(As_bt));
            
            xq = As_bt;
            zq = ts_t .* ones(size(xq)); % Expand ts_t to match xq dimensions
            
            % Aggregate data from sr_s0 (unsorted)
            X_all = [];
            Z_all = [];
            Y_all = [];
            for i = 1:numel(sr_s0)
                ts_t_ratio = sr_s0(i).ts_t_ratio; % Use actual field name from your struct
                x_vals = sr_s0(i).x_data(:);
                y_vals = sr_s0(i).y_data(:);
                X_all = [X_all; x_vals];
                Z_all = [Z_all; repmat(ts_t_ratio, numel(x_vals), 1)]; % Match x_vals size
                Y_all = [Y_all; y_vals];
            end
            
            % Build interpolant and evaluate
            sr_interpolant = scatteredInterpolant(X_all, Z_all, Y_all, 'linear', 'none');
            obj.sr_s0 = sr_interpolant(xq, zq);
        end
                
        function obj = interpretF(obj, F_data, wing)
            As_bt = obj.As ./ (obj.b1_c .* wing.cn) ./ obj.t_Upper;
            ts_t = obj.ts / obj.t_Upper;
            obj.F = zeros(size(As_bt));
            
            xq = As_bt;
            yq = ts_t.*ones(size(xq));
            
            % Collect all (x, y) points and their corresponding F values from contours
            X_all = [];
            Y_all = [];
            F_all = [];
            for i = 1:numel(F_data)
                F_val = F_data(i).F;
                x_vals = F_data(i).x_data(:); % Ensure column vector
                y_vals = F_data(i).y_data(:); 
                % Append data
                X_all = [X_all; x_vals];
                Y_all = [Y_all; y_vals];
                F_all = [F_all; repmat(F_val, numel(x_vals), 1)];
            end
            
            % Create interpolant (linear interpolation, no extrapolation)
            F_interpolant = scatteredInterpolant(X_all, Y_all, F_all, 'linear', 'none');
            
            % Interpolate at query points and assign to obj.F
            obj.F = F_interpolant(xq, yq);
        end

        function obj = createRibs(obj, Nrib, wing, D, distribution, alphaExp)
            % D is diameter of fuselage
            b = wing.b;
            R = D/2;
            obj.Nrib = Nrib;
        
            if mod(Nrib,2) ~= 0
                error("Number of ribs must be even")
            end
            
            switch distribution
                case "uniform"
                    % Linear (uniform) spacing
                    y = linspace(R, b/2, Nrib/2);
                    
                case "exponential"
                    % You can pick either approach:
                    % (a) standard logspace:
                    % y = logspace(log10(R), log10(b/2), Nrib/2);
        
                    % (b) power-logspace approach:
                    % y = customLogspace(R, b/2, Nrib/2, alphaExp);
        
                    % (c) direct geometric interpolation with alpha:
                    y = R .* ((b/2) / R) .^ (linspace(0,1,Nrib/2).^alphaExp);
        
                otherwise
                    error("Unknown distribution type")
            end
            
            obj.yrib = [-flip(y), y];
        end

        function plotRibs(obj,wing)
            hold on
            for i = 1:length(obj.yrib)
                thisy = obj.yrib(i);
                if thisy > wing.b/2
                    thisy = wing.b/2;
                elseif thisy < -wing.b/2
                    thisy = -wing.b/2;
                end
                xstart = wing.LEx(thisy);
                chord = wing.c_at_y(thisy);
                plot([xstart xstart + chord],[obj.yrib(i) obj.yrib(i)],'r-')
                hold on
            end

        end
        
        function obj = mapRibSpacing(obj, wing, D)
            % variable to assign: obj.L
            % y_wing   = wing.stripy;  % sample points on the wing
            % yrib     = obj.yrib;     % location of ribs, size(obj.Nrib)
            %
            % We want L to be size(y_wing).
            % For each y_wing(i):
            %   1) If y_wing(i) is between -D/2 and D/2, set L(i) = NaN.
            %   2) Else, find the two adjacent ribs in yrib that bracket y_wing(i).
            %      L(i) = y_rib(k+1) - y_rib(k).
        
            y_wing = wing.stripy;
            yrib   = obj.yrib;
        
            % Sort yrib if not already sorted
            yrib = sort(yrib);
        
            % Pre-allocate L
            L = zeros(size(y_wing));
        
            for i = 1:length(y_wing)
                y_i = y_wing(i);
        
                % 1) Check if y_i is within [-D/2, D/2]
                if (y_i >= -D/2) && (y_i <= D/2)
                    % If in the fuselage region, set L to NaN
                    L(i) = NaN;
                else
                    % 2) Otherwise, find the adjacent ribs
                    k = find(yrib <= y_i, 1, 'last');
        
                    if isempty(k)
                        % y_i < yrib(1)
                        L(i) = yrib(2) - yrib(1);
                    elseif k == length(yrib)
                        % y_i >= yrib(end)
                        L(i) = yrib(end) - yrib(end-1);
                    else
                        % y_i lies between yrib(k) and yrib(k+1)
                        L(i) = yrib(k+1) - yrib(k);
                    end
                end
            end
        
            obj.L = L;
        end

        function plotWingBox(obj,wing,D)
            fs = obj.b_c; %front spar percentage
            rs = obj.b_c + obj.c_c; %rear spar percentage
            R = D/2;
            froot = wing.LEx(R) + wing.c_at_y(R)*fs;
            rroot = wing.LEx(R) + wing.c_at_y(R)*rs;
            plot([froot rroot],[R R],'b')
            hold on
            ftip = wing.LEx(wing.b/2) + wing.c_at_y(wing.b/2)*fs;
            rtip = wing.LEx(wing.b/2) + wing.c_at_y(wing.b/2)*rs;
            plot([ftip rtip],[wing.b/2 wing.b/2],'b')
            hold on
            plot([ftip froot],[wing.b/2 R],'b')
            hold on
            plot([rtip rroot],[wing.b/2 R],'b')
        end

        function plotStringer(obj,wing,D)
            theta = wing.Lambdax_c(0,wing.b/2);
            fs = obj.b_c; %front spar percentage
            rs = obj.b_c + obj.c_c; %rear spar percentage
            R = D/2;
            froot = wing.LEx(R) + wing.c_at_y(R)*fs;
            rroot = wing.LEx(R) + wing.c_at_y(R)*rs;
            ftip = wing.LEx(wing.b/2) + wing.c_at_y(wing.b/2)*fs;
            rtip = wing.LEx(wing.b/2) + wing.c_at_y(wing.b/2)*rs;
            
            xf = @(y) froot + (ftip - froot)/(wing.b/2 - R)*(abs(y)-R);
            xr = @(y) rroot + (rtip - rroot)/(wing.b/2 - R)*(abs(y)-R);
            for i = 1:length(wing.stripy)-1
                y = wing.stripy(i);
                if y >= R
                    y2 = wing.stripy(i+1);
                    N = obj.N(i);
                    for j = 1:N
                        xs1 = xf(y) + j/(N+1)*(xr(y) - xf(y));
                        xs2 = xf(y2) + j/(N+1)*(xr(y2) - xf(y2));

                        ys1 = y;
                        ys2 = y2;

                        plot([xs1 xs2], [ys1 ys2], 'k')
                        hold on
                    end
                end
            end
        end

    end
end