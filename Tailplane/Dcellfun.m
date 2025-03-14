function t = Dcellfun(wb,wing,span,material,t1,n_rib)
    %return the thickness of D cell
    data1 = load('curvedPanel.mat');
    curvedPanel = data1.curvedPanel;
    data2 = load('curvedPanel2.mat');
    curvedPanel2 = data2.curvedPanel2;

    bref = span;
    a = wb.a;
    b = wb.b_c.*wing.cn;
    R1 = wb.R_c.*wing.cn;

    % Number of elements in the struct
    numEntries = length(curvedPanel2);

    % Initialize a struct to store polyfit coefficients
    fitResults = struct();
    
    % Initialize figure
    % figure; hold on;
    % colors = lines(numEntries); % Different colors for each line
    
    for i = 1:numEntries
        % Extract x_data and y_data
        x = curvedPanel2(i).x_data;
        y = curvedPanel2(i).y_data;
        b_a = curvedPanel2(i).b_a; % Get the b_a value
    
        % Fit an 8th-degree polynomial
        p = polyfit(x, y, 8);
    
        % Store the polyfit coefficients in the struct
        fitResults(i).b_a = b_a;
        fitResults(i).coeffs = p;  % Store all polynomial coefficients
    
        % Generate fitted y-values
        y_fit = polyval(p, x);
        % 
        % % Plot original data points
        % scatter(x, y, 50, colors(i, :), 'filled');
        % 
        % % Plot fitted polynomial curve
        % plot(x, y_fit, 'Color', colors(i, :), 'LineWidth', 2, 'DisplayName', sprintf('b_a = %.1f', b_a));
        % 
        % %Display polynomial equation in console
        % fprintf('b_a = %.1f: Polynomial Coefficients = [', b_a);
        % fprintf('%.4f ', p);
        % fprintf(']\n');
    end
    
    % % Display a legend for the plot
    % legend('show');
    % % Formatting plot
    % xlabel('x data');
    % ylabel('y data');
    % title('8th-Degree Polynomial Fit for Each b_a Value');
    % legend('show');
    % grid on;
    % hold off;

    % % Save the fitted results for later use
    % save('polyfit_results_curvedPanel2.mat', 'fitResults');



    % figure; hold on;

    
    % Number of elements in the struct
    numEntries = length(curvedPanel);
    
    % Initialize a struct to store polyfit coefficients
    fitResultsa_b = struct();
    
    % Initialize figure
    % figure; hold on;
    % colors = lines(numEntries); % Different colors for each line
    
    for i = 1:numEntries
        % Extract x_data and y_data
        x = curvedPanel(i).x_data;
        y = curvedPanel(i).y_data;
        a_b = curvedPanel(i).a_b; % Get the b_a value
    
        % Fit an 8th-degree polynomial
        p = polyfit(x, y, 8);
    
        % Store the polyfit coefficients in the struct
        fitResultsa_b(i).a_b = a_b;
        fitResultsa_b(i).coeffs = p;  % Store all polynomial coefficients
    
        % Generate fitted y-values
        y_fit = polyval(p, x);
    
        % % Plot original data points
        % scatter(x, y, 50, colors(i, :), 'filled');
        % 
        % % Plot fitted polynomial curve
        % plot(x, y_fit, 'Color', colors(i, :), 'LineWidth', 2, 'DisplayName', sprintf('a_b = %.3g', a_b));
        % 
        % %Display polynomial equation in console
        % fprintf('a_b = %.1f: Polynomial Coefficients = [', a_b);
        % fprintf('%.4f ', p);
        % fprintf(']\n');
    end
    
    % % Display a legend for the plot
    % legend('show');
    % % Formatting plot
    % xlabel('x data');
    % ylabel('y data');
    % title('8th-Degree Polynomial Fit for Each a/b Value');
    % legend('show');
    % grid on;
    % hold off;
    
    % % Save the fitted results for later use
    % save('polyfit_results_curvedPanela_b.mat', 'fitResults');
    
    n = size(a,1);
    
    Ks = [];
    display(["size of a:", size(a)])
    for i = 1:n
        if b(i) > a(i)
            Ks(i) = compute_Ks_b_a(a(i), b(i), R1(i), t1, fitResults);
        elseif a(i) == b(i)
            error("lmao")
        else 
            Ks(i) = compute_Ks_a_b(a(i), b(i), R1(i), t1, fitResultsa_b);
        end
    end

    E = material.E;
    display(["size of Ks:", size(Ks)])
    tau_cr = Ks.*E.*(t1./b).^2;

    tau_tres = 125*1e6;
    t = sqrt(tau_tres.*b.^2./Ks./E);

    % % weight check
    % cellArea1 = 11000;
    % Al_rho = 2.7;
    % P_rib = 0.5 * cellArea1 * n_rib * Al_rho * 1000 / 10^9;
    % skin = t / 1000 * b * 2 * bref * Al_rho * 1000;

    function Ks = compute_Ks_b_a(a, b, R1, t1, fitResults)
        % Ensure b > a, otherwise return an error
        if b <= a
            error('Invalid input: b must be greater than a.');
        end
    
        % Compute x_val
        x_val = a / sqrt(R1 * t1);
        %disp(['x_val = ', num2str(x_val)]);  % Debugging output
    
        if x_val < 0 | x_val > 20
            error('Invalid input: Outside of ESDU scope');
        end
    
        % Extract all available b/a values from fitResults
        b_a_values = [fitResults.b_a];
    
        % Determine the appropriate polynomial fit based on b/a conditions
        b_a_query = b / a;  % Compute b/a ratio
        %disp(b_a_query)
        %disp(['b/a = ', num2str(b_a_query)]);  % Debugging output
    
        if b_a_query < 1
            error('b/a < 1 is not supported.');
        elseif b_a_query > 5
            % Use the polynomial fit where b/a = 5
            [~, idx] = min(abs(b_a_values - 5));
            Ks = polyval(fitResults(idx).coeffs, x_val);
        else
            % Find the two closest b/a values surrounding b_a_query
            lower_idx = find(b_a_values <= b_a_query, 1, 'last'); % Lower bound
            upper_idx = find(b_a_values >= b_a_query, 1, 'first'); % Upper bound
            %disp(['lower_idx = ', num2str(lower_idx), ', upper_idx = ', num2str(upper_idx)])  % Debugging output
    
            % If exact match found, use it directly
            if b_a_values(lower_idx) == b_a_query
                Ks = polyval(fitResults(lower_idx).coeffs, x_val);
            else
                % Evaluate both polynomial fits at x_val
                Ks_lower = polyval(fitResults(lower_idx).coeffs, x_val);
                Ks_upper = polyval(fitResults(upper_idx).coeffs, x_val);
    
                % Perform linear interpolation between the two polynomial outputs
                w = (b_a_query - b_a_values(lower_idx)) / (b_a_values(upper_idx) - b_a_values(lower_idx));
                Ks = (1 - w) * Ks_lower + w * Ks_upper;
            end
        end
    
        % Debugging output to ensure Ks is returned
        %disp(['Computed Ks = ', num2str(Ks)]);
    end


    function Ks = compute_Ks_a_b(a, b, R1, t1, fitResultsa_b)
        % Ensure b > a, otherwise return an error
        if a <= b
            error('Invalid input: b must be greater than a.');
        end
    
        % Compute x_val
        x_val = b / sqrt(R1 * t1);
        %disp(['x_val = ', num2str(x_val)]);  % Debugging output
    
        if x_val < 0 || x_val > 20
            error('Invalid input: Outside of ESDU scope');
        end
    
        % Extract all available a/b values from fitResultsa_b
        a_b_values = [fitResultsa_b.a_b];
    
        % Determine the appropriate polynomial fit based on a/b conditions
        a_b_query = a / b;  % Compute b/a ratio
        disp(['a/b = ', num2str(a_b_query)]);  % Debugging output
    
        if a_b_query < 1
            error('a/b < 1 is not supported.');
        elseif a_b_query > 1e100
            % Use the polynomial fit where a/b = 5
            [~, idx] = min(abs(a_b_values - 5));
            Ks = polyval(fitResultsa_b(idx).coeffs, x_val);
        else
            % Find the two closest a/b values surrounding b_a_query
            lower_idx = find(a_b_values <= a_b_query, 1, 'last'); % Lower bound
            upper_idx = find(a_b_values >= a_b_query, 1, 'first'); % Upper bound
            disp(['lower_idx = ', num2str(lower_idx), ', upper_idx = ', num2str(upper_idx)]);  % Debugging output
    
            % If exact match found, use it directly
            if a_b_values(lower_idx) == a_b_query
                Ks = polyval(fitResults(lower_idx).coeffs, x_val);
            else
                % Evaluate both polynomial fits at x_val
                Ks_lower = polyval(fitResultsa_b(lower_idx).coeffs, x_val);
                Ks_upper = polyval(fitResultsa_b(upper_idx).coeffs, x_val);
    
                % Perform linear interpolation between the two polynomial outputs
                w = (a_b_query - a_b_values(lower_idx)) / (a_b_values(upper_idx) - a_b_values(lower_idx));
                Ks = (1 - w) * Ks_lower + w * Ks_upper;
            end
        end
    
        % Debugging output to ensure Ks is returned
        disp(['Computed Ks = ', num2str(Ks)]);
    end
end