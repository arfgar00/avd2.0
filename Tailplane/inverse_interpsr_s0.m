function [As_bt, ts_t] = inverse_interpsr_s0(sr_target)
    % Load the data and build the interpolant (same as in interpsr_s0)
    data = load("scr_s0.mat");
    data = data.scr_s0;

    % Aggregate data from scr_s0.mat
    X_all = [];
    Z_all = [];
    Y_all = [];
    for i = 1:numel(data)
        ts_t_ratio = data(i).ts_t_ratio;
        x_vals = data(i).x_data(:);
        y_vals = data(i).y_data(:);
        X_all = [X_all; x_vals];
        Z_all = [Z_all; repmat(ts_t_ratio, numel(x_vals), 1)];
        Y_all = [Y_all; y_vals];
    end

    % Create the scattered interpolant
    sr_interpolant = scatteredInterpolant(X_all, Z_all, Y_all, 'linear', 'none');

    % Define the objective function for optimization.
    % x(1) corresponds to As_bt and x(2) corresponds to ts_t.
    obj = @(x) abs(sr_interpolant(x(1), x(2)) - sr_target);

    % Choose an initial guess, e.g., the mean values of the data.
    x0 = [mean(X_all), mean(Z_all)];

    % Use fminsearch to find the pair (As_bt, ts_t) that minimizes the error.
    options = optimset('Display', 'off');
    sol = fminsearch(obj, x0, options);

    % Return the found solution.
    As_bt = sol(1);
    ts_t = sol(2);
end
