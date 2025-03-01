function [thickness_buckling, thickness_yield]= Ribfun(t_e,M,wing,wb,material,tr0)
    
    %tr0 guess rib thickness
    E = material.E;
    chord = wing.cn;
    h_c = wb.b2_c.*wing.cn;
    s = wb.a;
    sigma_y = material.sigma_y;

    I = chord .* (t_e)^3 / 12 + chord .* (t_e) .* (h_c./2).^2;
    C = (M.^2 .* s .* h_c .* ...
        (t_e) .* chord ./ 2) ./ ...
        (E .* I.^2);
    CS = (C ./ tr0) .* chord;
    
    
    thickness_buckling = ((C ./ chord) ./ (3.62 * E) .* h_c.^2).^(1/3);
    BS = (C ./ thickness_buckling) ./ chord;
    
    thickness_yield = (C ./ sigma_y) ./ chord;
end