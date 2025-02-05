function sigma_cr = findFStress(wb,M,wing,material)
    %inputs: wb: wingbox object
    %        M:  moment, an array varies in y direction
    %        wing: wingGeometry object
    %using the method in excel skin

    F = wb.F; % sigma_critical/sigma_0, an array varies in y direction
    c = wb.c_c.*wing.cn;
    b2 = wb.b2_c.*wing.cn;
    N = M./c./b2;
    Et = material.E;
    L = wb.L;
    sigma_cr = F.*sqrt(N.*Et./L);
end