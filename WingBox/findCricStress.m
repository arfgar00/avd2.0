function sigma_cr = findCricStress(wb,M,wing,material)
    %inputs: wb: wingbox object
    %        M:  moment, an array varies in y direction
    %        wing: wingGeometry object
    %using the method in excel skin

    sr_s0 = wb.sr_s0; % sigma_critical/sigma_0, an array varies in y direction
    c = wb.c_c.*wing.cn;
    b2 = wb.b2_c.*wing.cn;
    N = M./c./b2;
    E = material.E;
    %sigma0 = N./wb.t2;
    sigma0 = 3.62*E*(wb.t2./b2).^2;
    sigma_cr = sigma0.*sr_s0;
    
end