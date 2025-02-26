function F = interpF(As_bt, ts_t)
    data = load("Farror.mat");
    data = data.F_data;
    xq = As_bt;
    yq = ts_t.*ones(size(xq));
    
    % Collect all (x, y) points and their corresponding F values from contours
    X_all = [];
    Y_all = [];
    F_all = [];
    for i = 1:numel(data)
        F_val = data(i).F;
        x_vals = data(i).x_data(:); % Ensure column vector
        y_vals = data(i).y_data(:); 
        % Append data
        X_all = [X_all; x_vals];
        Y_all = [Y_all; y_vals];
        F_all = [F_all; repmat(F_val, numel(x_vals), 1)];
    end
    
    % Create interpolant (linear interpolation, no extrapolation)
    F_interpolant = scatteredInterpolant(X_all, Y_all, F_all, 'linear', 'none');
    
    % Interpolate at query points and assign to obj.F
    F = F_interpolant(xq, yq);
end