%shear load & torque
a = 2;
c = 1;
b2 = 0.2;

Ks = 0; % to edit ie digiscan and look up value for a/b2

E_Al2024 = 70;
V = 20000;
T = 10000;
q0 = T/2/b2/c/1000;
q2 = V/2/b2/1000;
q_FS = q0 + q2;
q_RS = q2 - q0;

t_FS = (q_FS * 1000 * b2 / Ks / (E_Al2024 * 10^9))^(1/3)*1000;
t_RS = (q_FS * 1000 * b2 / Ks / (E_Al2024 * 10^9))^(1/3)*1000;

tau_F = q_FS/t_FS;
tau_R = q_RS/t_RS;

%combined loading
t2 = 1.05;
q0 = 25;
tau_0 = q0/t2;
b1 = 50;

tau_cr = 8.1*E_Al2024*10^9*(t2/b1)^2/10^6;
tau_cr_tresca = 125;

sigma_0 = 111.18;
sigma_cr = 133.4105;

R_c = sigma_0/sigma_cr;
R_s = tau_cr_tresca/tau_0;
R_b = 0;

if R_s^2 + R_c > 0.99
    warning('Combined loading restriction not fulfilled: R_s^2 + R_c exceeds 0.99.');
end

if R_s^2 + R_b^2 > 0.99
    warning('Combined loading restriction not fulfilled: R_s^2 + R_c exceeds 0.99.');
end

if R_b^1.75 + R_c > 0.99
    warning('Combined loading restriction not fulfilled: R_s^2 + R_c exceeds 0.99.');
end

if R_s^2 + (R_c + R_b)^2 > 0.99
    warning('Combined loading restriction not fulfilled: R_s^2 + R_c exceeds 0.99.');
end


