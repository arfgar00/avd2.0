function err = objective_Upanel(x, ci, N, E)
    % Extract variables
    ts = x(1);
    t = x(2);
    b1_c = x(3);
    hs = x(4);
    L = x(5);
    d_h = 0.3;%aspect ratio
    b1 = b1_c*ci;

    %interpolate from the graph
    As_bt = ts/t*hs/b1*(1 + 2/d_h);
    ts_t = ts/t;
    scr_s0 = interpsr_s0(As_bt, ts_t);

    F = interpF(As_bt, ts_t);

    
    %calc err
    t = scr_s0/F*sqrt(N*L/E);
    err = (1/(b1_c) - 1/ci^2)*ts*hs*(1 + 2*d_h) + t/ci;
end