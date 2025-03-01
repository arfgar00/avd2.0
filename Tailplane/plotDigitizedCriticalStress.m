figure(1)
clf;
load("scr_s0.mat")

for i = 1:length(digitized_data)
    [sorted_x, sort_order] = sort(digitized_data(i).x_data);
    sorted_y = digitized_data(i).y_data(sort_order);
    % Define a dense grid of x-values for smooth plotting
    xi = linspace(min(sorted_x), max(sorted_x), 1000); % 1000 points for smoothness
    
    % Interpolate using cubic spline
    yi = interp1(sorted_x, sorted_y, xi, 'spline');
    % Compute spline coefficients
    pp = spline(sorted_x, sorted_y); % Returns piecewise polynomial (pp) form
    
    % Evaluate spline at desired x-values
    xi = linspace(min(sorted_x), max(sorted_x), 1000);
    yi = ppval(pp, xi);
    plot(xi, yi, 'LineWidth', 2, 'DisplayName', sprintf('ts/t = %.1f (Spline)', digitized_data(i).ts_t_ratio));
    hold on;
end

figure(2)
clf;
load("Farror.mat")

for i = 1:length(F_data)
    plot(F_data(i).x_data, F_data(i).y_data); hold on
    
    % Pick the midpoint to place label
    midIdx = floor(numel(F_data(i).x_data)/2);
    xm = F_data(i).x_data(midIdx);
    ym = F_data(i).y_data(midIdx);
    
    % Create a label using F_data(i).F
    txt = sprintf('F=%.2f', F_data(i).F);
    text(xm, ym, txt, ...
        'FontSize',8, ...
        'BackgroundColor','w', ...   % white background for visibility
        'Margin',1);
end
xlabel("$A_s/bt$","Interpreter","latex")
ylabel("$t_s/t_2$","Interpreter","latex")