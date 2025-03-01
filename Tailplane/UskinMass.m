function [M_skin_total,M_Uskin,M_Uskin_stringer,b1_c] = UskinMass(rho_Uskin, wing_opt, wb,wingidx)
    V_Uskin = 2*trapz(wing_opt.stripy(wingidx), wb.t_Upper(wingidx).*wb.c_c.*wing_opt.cn(wingidx)); %plate
    As_Upper = wb.ts_Upper.*wb.hs.*(1 + 2*0.3);
    b1_c = wb.bs_Upper./(wing_opt.cn);
    AsN = As_Upper(wingidx)./(b1_c(wingidx)./wb.c_c);
    AsN(isnan(AsN)) = 0; %condition for integration
    V_Uskin_stringer = 2*trapz(wing_opt.stripy(wingidx), AsN);

    V_Uskin_total = V_Uskin + V_Uskin_stringer;
    M_skin_total = rho_Uskin*V_Uskin_total;

    M_Uskin = V_Uskin*rho_Uskin;
    M_Uskin_stringer=  V_Uskin_stringer * rho_Uskin;
end