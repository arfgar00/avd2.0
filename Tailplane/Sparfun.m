function [t_FS, t_RS] = Sparfun(wb, wing,material,V,T,wbStress)
    data = load('webBuckling.mat'); % Load as struct
    webBuckling = data.webBuckling; % Extract the specific variable
    %shear load: V
    %torque: T
    
    a = wb.a;
    c = wb.c_c.*wing.cn;
    b2 = wb.b2_c.*wing.cn;
    
    ratio = a./b2;
    
    Ks = interp1(webBuckling.x_data,webBuckling.y_data,ratio); % to edit ie digiscan and look up value for a/b2

    E = material.E;

    q0 = T./2./b2./c;
    q2 = V./2./b2;
    q_FS = q0 + q2; %front spar
    q_RS = q2 - q0; %rear spar

    t_FS = (q_FS .* b2 ./ Ks ./ E).^(1/3);%thickness front spar
    t_RS = (q_FS .* b2 ./ Ks ./ E).^(1/3);%thickness

    tau_F = q_FS./t_FS;%shear stress front
    tau_R = q_RS./t_RS;

    %combined loading
    t2 = wb.t2;
    tau_0 = q0./t2;
    b1 = wb.b1_c .* wing.cn;

    tau_cr = 8.1*E.*(t2./b1).^2;
    tau_cr_tresca = 125*1e6;

    sigma_0 = wbStress.s_0_global_skin;
    sigma_cr = wbStress.s_cr_global_skin;

    R_c = sigma_0/sigma_cr;
    R_s = tau_cr_tresca/tau_0;
    R_b = 0;

    if R_s^2 + R_c > 0.99
        error('Combined loading restriction not fulfilled: R_s^2 + R_c exceeds 0.99.');
    elseif R_s^2 + R_b^2 > 0.99
        error('Combined loading restriction not fulfilled: R_s^2 + R_c exceeds 0.99.');
    elseif R_b^1.75 + R_c > 0.99
        error('Combined loading restriction not fulfilled: R_s^2 + R_c exceeds 0.99.');
    elseif R_s^2 + (R_c + R_b)^2 > 0.99
        error('Combined loading restriction not fulfilled: R_s^2 + R_c exceeds 0.99.');
    else
    end


end
