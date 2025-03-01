function [M] = ribMass(wb, wing_opt, rho)
    crib = [];
    for i = 1:length(wb.yrib)
        crib(i) = wing_opt.c_at_y(wb.yrib(i));
    end
    Ar = (wb.c_c .* wb.b2_c .*crib.^2);
    V_ribs = sum(wb.tr.* Ar');
    M = V_ribs*rho;
end