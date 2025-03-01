function [M_skin_total,M_Lskin,M_Lskin_stringer,b1_c_Lower] = LskinMass(rho_Lskin, wing_opt, wb, wingidx)
    V_Lskin = 2*trapz(wing_opt.stripy(wingidx), wb.t_Lower(wingidx).*wb.c_c.*wing_opt.cn(wingidx)); %plate
    As_Lower = wb.ts_Lower.*wb.hs_Lower.*(1 + 2*0.3);
    b1_c_Lower = wb.bs_Lower./(wing_opt.cn);
    AsN = As_Lower(wingidx)./(b1_c_Lower(wingidx)./wb.c_c);
    AsN(isnan(AsN)) = 0; %condition for integration
    V_Lskin_stringer = 2*trapz(wing_opt.stripy(wingidx), AsN);
    V_Lskin_total = V_Lskin + V_Lskin_stringer;
    M_skin_total = rho_Lskin*V_Lskin_total;

    M_Lskin = V_Lskin*rho_Lskin;
    M_Lskin_stringer=  V_Lskin_stringer * rho_Lskin;
end