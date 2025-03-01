function Ks = interpSingleKs(R,t , a, b)
    if a/b > 1
        data = load("curvedPanel.mat");
        data = data.curvedPanel;
        zq = a/b;
        xq = b/sqrt(R*t);
        load1 = 1;
    else
        data = load("curvedPanel2.mat");
        data = data.curvedPanel2;
        zq = b/a;
        xq = a/sqrt(R*t);
        load1 = 0;
    end

    X_all = [];
    Z_all = [];
    Y_all = [];
    for i = 1:numel(data)
        if load1 == 1
            Z_val = data(i).a_b;
        else
            Z_val = data(i).b_a;
        end
        x_vals = data(i).x_data(:);
        y_vals = data(i).y_data(:);
        X_all = [X_all; x_vals];
        Z_all = [Z_all; repmat(Z_val, numel(x_vals), 1)];
        Y_all = [Y_all; y_vals];
    end

    pts = [X_all, Z_all];
    [uniquePts, ~, ic] = unique(pts, 'rows');
    Y_unique = accumarray(ic, Y_all, [], @mean);

    interpolant = scatteredInterpolant(uniquePts(:,1), uniquePts(:,2), Y_unique, 'linear', 'none');
    Ks = interpolant(xq, zq);
end