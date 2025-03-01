function [ts_t, As_bt] = inverse_interpF(F_target)
    data = load("Farror.mat");
    data = data.F_data;
    
    % Collect all points and F values as in original function
    X_all = [];
    Y_all = [];
    F_all = [];
    for i = 1:numel(data)
        F_val = data(i).F;
        x_vals = data(i).x_data(:);
        y_vals = data(i).y_data(:);
        X_all = [X_all; x_vals];
        Y_all = [Y_all; y_vals];
        F_all = [F_all; repmat(F_val, numel(x_vals), 1)];
    end
    
    % Create interpolant
    F_interpolant = scatteredInterpolant(X_all, Y_all, F_all, 'linear', 'none');
    
    % Check if F_target exists in data with a tolerance
    F_vals = [data.F];
    tol = 1e-6;
    idx = find(abs(F_vals - F_target) < tol);
    if ~isempty(idx)
        % Return all points from the matching F_val contours
        As_bt = cell(1, numel(idx));
        ts_t = cell(1, numel(idx));
        for k = 1:numel(idx)
            As_bt{k} = data(idx(k)).x_data(:)';
            ts_t{k} = data(idx(k)).y_data(:)';
        end
        % If only one contour, return as vectors
        if numel(idx) == 1
            As_bt = As_bt{1};
            ts_t = ts_t{1};
        end
        return;
    end
    
    % Determine grid range
    x_min = min(X_all);
    x_max = max(X_all);
    y_min = min(Y_all);
    y_max = max(Y_all);
    
    % Create a grid with a reasonable resolution
    num_points = 500; % Increase for higher accuracy
    x_grid = linspace(x_min, x_max, num_points);
    y_grid = linspace(y_min, y_max, num_points);
    [X_grid, Y_grid] = meshgrid(x_grid, y_grid);
    
    % Evaluate interpolant on the grid
    F_grid = F_interpolant(X_grid, Y_grid);
    
    % Compute contour at F_target
    contour_level = F_target;
    C = contourc(x_grid, y_grid, F_grid, [contour_level, contour_level]);
    
    % Parse contourc output into segments
    segments = {};
    idx = 1;
    while idx < size(C, 2)
        level = C(1, idx);
        num_points_segment = C(2, idx);
        next_idx = idx + num_points_segment + 1;
        if next_idx > size(C, 2)
            break;
        end
        x_segment = C(1, idx+1:idx+num_points_segment);
        y_segment = C(2, idx+1:idx+num_points_segment);
        segments{end+1} = [x_segment; y_segment];
        idx = next_idx;
    end
    
    % Prepare output based on segments
    if isempty(segments)
        As_bt = [];
        ts_t = [];
    else
        As_bt = cell(1, length(segments));
        ts_t = cell(1, length(segments));
        for i = 1:length(segments)
            As_bt{i} = segments{i}(1, :);
            ts_t{i} = segments{i}(2, :);
        end
        % Return as vectors if single segment
        if length(segments) == 1
            As_bt = As_bt{1};
            ts_t = ts_t{1};
        end
    end
end