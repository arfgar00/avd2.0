function [tau_n, tau_r, tau_w, sigmax] = boxfun(wb, wing, P, M, material)
    D = wb.b2_c.*wing.cn;
    t1 = wb.t1_spar;
    t2 = wb.t2_spar;
    b = wb.b_spar;
    c = D./2;
    
    Ixx = 1/6 .* b .* t1.^3 + 1/2 .* b .* t1 .* (D - t1).^2 + (1/12) .* t2 .* (D - 2.*t1).^3;

    sigmax = M .* c ./ Ixx;
    
    An = wb.An; %n for nose, D cell
    Ar = wb.Ar; %r for rear, main cell
    t_w = wb.tF; %confusing
    tn = wb.t1;
    tr = wb.t2;
    d = (0.5-wb.b_c).*wing.cn;      %distance to stress centre, for front spar
    Sn = wb.Sn;
    Sr = wb.Sr;
    h_w = wb.b2_c.*wing.cn;
    G = material.G;

    torque = P .* d;
    tau_n = zeros(size(An));
    tau_r = zeros(size(An));
    tau_w = zeros(size(An));
    for i = 1:length(An)
        matrix1 = [2*An(i) 2*Ar(i) 0; 
                   h_w(i) -h_w(i) -h_w(i); 
            1/2/An(i)/G*Sn(i)/tn  -1/2/Ar(i)/G*Sr(i)/tr   h_w(i)/t_w*1/2/An(i)/G+h_w(i)/t_w*1/2/Ar(i)/G];
        vector1 = [torque(i); P(i); 0];
        solution = matrix1 \ vector1;
        tau_n(i) = solution(1) / tn;
        tau_r(i) = solution(2) / tr;
        tau_w(i) = solution(3) / t_w;
    end

    % %add checks here
    % P = (solution(1) - solution(2) - solution(3))*torsion_parameters.h_w;
    % T = 2*torsion_parameters.An*solution(1) + 2*torsion_parameters.Ar*solution(2);
    % theta_xn = 1/2/torsion_parameters.An/torsion_parameters.G * (torsion_parameters.Sn/torsion_parameters.tn * solution(1) + torsion_parameters.h_w/torsion_parameters.t_w * solution(3));
    % theta_xr = 1/2/torsion_parameters.Ar/torsion_parameters.G * (torsion_parameters.Sr/torsion_parameters.tr * solution(2) - torsion_parameters.h_w/torsion_parameters.t_w * solution(3));
    % 
    % 
    
end