function sigma_cr = findCricStress(wb,M,wing)
    %inputs: wb: wingbox object
    %        M:  moment, an array varies in y direction
    %        wing: wingGeometry object
    %using the method in excel skin

    sr_s0 = wb.sr_s0; % sigma_critical/sigma_0, an array varies in y direction
    c = wb.c_c.*wing.cn;
    b2 = wb.b2_c.*wing.cn;
    N = M./c./b2;
    sigma0 = N./wb.t2;
    sigma_cr = sigma0.*sr_s0;
end