function sr_s0 = interpKsSkinStringer(h_b, ts_t)
    data = load("K_s.mat");
    data = data.scr_s0;

    xq = h_b;
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
    for i = 1:length(sr_s0)
        if isnan(sr_s0(i))
            % Find the nearest non-NaN value
            distances = abs(sr_s0 - sr_s0(i));
            distances(isnan(sr_s0)) = Inf; % Ignore NaN values
            [~, nearest_idx] = min(distances);
            sr_s0(i) = sr_s0(nearest_idx);
        end
    end
end