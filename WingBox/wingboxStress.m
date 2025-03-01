classdef wingboxStress
    properties
        s_cr_local_skin
        s_cr_global_skin
        wb
        s_0_global_skin
    end

    methods
        function obj = findGlobalStress(obj, M,wing,material, t_panel)
            %inputs: wb: wingbox object
            %        M:  moment, an array varies in y direction
            %        wing: wingGeometry object
            %using the method in excel skin
            wb = obj.wb;
            sr_s0 = wb.sr_s0; % sigma_critical/sigma_0, an array varies in y direction
            c = wb.c_c.*wing.cn;
            b2 = wb.b2_c.*wing.cn;
            N = M./c./b2;
            E = material.E;
            sigma_0 = wb.Ks_Upper.*E.*(wb.t_Upper./(wb.b1_c.*wing.cn)).^2;
            sigma_cr = sigma_0.*sr_s0;
            obj.s_0_global_skin = sigma_0;
            obj.s_cr_global_skin = sigma_cr;
        end
        

        function obj = findLocalStress(obj,M,wing,material)
            %inputs: wb: wingbox object
            %        M:  moment, an array varies in y direction
            %        wing: wingGeometry object
            %using the method in excel skin
            wb = obj.wb;
            F = wb.F; % sigma_critical/sigma_0, an array varies in y direction
            c = wb.c_c.*wing.cn;
            b2 = wb.b2_c.*wing.cn;
            N = M./c./b2;
            Et = material.E;
            L = wb.L;
            sigma_cr = F.*sqrt(N.*Et./L);
            obj.s_cr_local_skin = sigma_cr;
        end
    end
end