function sr_s0 = interpsr_s0(As_bt, ts_t)
    data = load("scr_s0.mat");
    data = data.scr_s0;

    xq = As_bt;
    zq = ts_t .* ones(size(xq)); % Expand ts_t to match xq dimensions

    % Aggregate data from sr_s0 (unsorted)
    X_all = [];
    Z_all = [];
    Y_all = [];
    for i = 1:numel(data)
        ts_t_ratio = data(i).ts_t_ratio;
        x_vals = data(i).x_data(:);
        y_vals = data(i).y_data(:);
        X_all = [X_all; x_vals];
        Z_all = [Z_all; repmat(ts_t_ratio, numel(x_vals), 1)]; % Match x_vals size
        Y_all = [Y_all; y_vals];
    end
    
    % Build interpolant and evaluate
    sr_interpolant = scatteredInterpolant(X_all, Z_all, Y_all, 'linear', 'none');
    sr_s0 = sr_interpolant(xq, zq);
end