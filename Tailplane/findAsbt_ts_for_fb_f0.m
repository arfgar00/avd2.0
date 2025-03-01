function [As_bt_vals, ts_t_vals] = findAsbt_ts_for_fb_f0(targetFBF0, nAsbt, nTs_t)
    % FINDASBT_TS_FOR_FB_F0  Find all (As/bt, ts/t) combinations giving fb/f0 = targetFBF0
    %
    %   [As_bt_vals, ts_t_vals] = findAsbt_ts_for_fb_f0(targetFBF0, nAsbt, nTs_t)
    %
    % Inputs:
    %   targetFBF0  - Desired fb/f0 value (e.g. 1.0).
    %   nAsbt       - Number of samples in As/bt direction for mesh (e.g. 100).
    %   nTs_t       - Number of samples in ts/t direction for mesh (e.g. 100).
    %
    % Outputs:
    %   As_bt_vals  - All As/bt points on the fb/f0 = targetFBF0 contour.
    %   ts_t_vals   - Matching ts/t points on that same contour.
    %
    % Example:
    %    % Suppose you want fb/f0 = 1.0 and a 100x100 sampling:
    %    [As_bt_vals, ts_t_vals] = findAsbt_ts_for_fb_f0(1.0, 100, 100);
    %    plot(As_bt_vals, ts_t_vals, 'o');
    %

    %% 1. Load raw data (same as in interpsr_s0)
    rawData = load("scr_s0.mat");
    rawData = rawData.scr_s0;

    % Collect all raw (X, Z, Y) points for building the interpolant
    X_all = [];  % This will store all As/bt
    Z_all = [];  % This will store all ts/t
    Y_all = [];  % This will store fb/f0
    for i = 1:numel(rawData)
        ts_t_ratio = rawData(i).ts_t_ratio;
        x_vals     = rawData(i).x_data(:);  % As/bt
        y_vals     = rawData(i).y_data(:);  % fb/f0
        nPts       = numel(x_vals);

        X_all = [X_all; x_vals];
        Z_all = [Z_all; repmat(ts_t_ratio, nPts, 1)];
        Y_all = [Y_all; y_vals];
    end

    % Build the scatteredInterpolant (as in interpsr_s0)
    sr_interpolant = scatteredInterpolant(X_all, Z_all, Y_all, 'linear', 'none');

    %% 2. Define a mesh over the valid range of (As/bt, ts/t)
    x_min = min(X_all);  x_max = max(X_all);
    z_min = min(Z_all);  z_max = max(Z_all);

    % Generate uniform grids
    As_bt_grid = linspace(x_min, x_max, nAsbt);
    ts_t_grid  = linspace(z_min, z_max, nTs_t);

    [Xq, Zq] = meshgrid(As_bt_grid, ts_t_grid);

    %% 3. Evaluate fb/f0 on the mesh
    fb_f0_vals = sr_interpolant(Xq, Zq);

    %% 4. Extract the contour for fb/f0 == targetFBF0
    % The contour command returns a matrix C with all contour segments.
    % Each segment has a "header" with the level and # points, followed by
    % columns of (x, y). We parse that here.

    % For clarity, let's treat Xq as "As/bt" and Zq as "ts/t" for plotting:
    [C, ~] = contour(Xq, Zq, fb_f0_vals, [targetFBF0, targetFBF0]);

    % Initialize arrays for solutions
    As_bt_vals = [];
    ts_t_vals  = [];

    % Parse the contour matrix C
    idx = 1;
    while idx < size(C, 2)
        thisLevel  = C(1, idx);
        numPoints  = C(2, idx);
        contourPts = C(:, (idx+1):(idx + numPoints));  % columns of [x; y]

        idx = idx + numPoints + 1;  % move to next segment

        if abs(thisLevel - targetFBF0) < 1e-10
            % This segment corresponds to the desired fb/f0 = targetFBF0
            As_bt_vals = [As_bt_vals, contourPts(1, :)];
            ts_t_vals  = [ts_t_vals,  contourPts(2, :)];
        end
    end

    % That’s it! As_bt_vals and ts_t_vals are all the solutions on the fb/f0=targetFBF0 curve.
end