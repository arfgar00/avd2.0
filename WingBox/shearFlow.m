function  [q_wb_num, q_D_num] = shearFlow(chord, wb, Tau, tU, tL, tF, tR, tD)

    % Define symbolic variables
    syms q_wb q_D                   % Unknowns
    syms A_wb A_D c_wb t_U t_L h_wb t_F t_R b_D t_D tau % Knowns
    
    % Define equations
    % Equation 1: Torque equilibrium
    eq1 = tau == 2*A_wb*q_wb + 2*A_D*q_D;
    
    % Equation 2: Compatibility of twist rate
    eq2 = (1/(2*A_wb)) * (c_wb*(q_wb/t_U + q_wb/t_L) + h_wb*((q_wb - q_D)/t_F + q_wb/t_R)) ...
           == ...
           (1/(2*A_D)) * (h_wb*((q_D - q_wb)/t_F) + b_D*q_D/t_D);
    
    % Solve the system symbolically
    sol = solve([eq1, eq2], [q_wb, q_D]);
    q_wb_sol = sol.q_wb;
    q_D_sol = sol.q_D;
    
    % Substitute numerical values (replace with your actual values)
    % Example values:
    A_wb_val = wb.c_c.*wb.b2_c.*chord^2;     % Wingbox enclosed area [m²]
    A_D_val = pi*wb.R_c^2*chord^2/2;      % D-cell enclosed area [m²]
    c_wb_val = chord*wb.c_c;     % Wingbox chord length [m]
    t_U_val = tU;    % Upper skin thickness [m]
    t_L_val = tL;    % Lower skin thickness [m]
    h_wb_val = wb.b2_c*chord;    % Wingbox height [m]
    t_F_val = tF;    % Front spar thickness [m]
    t_R_val = tR;    % Rear spar thickness [m]
    b_D_val = pi*wb.R_c*chord;      % D-cell circumference [m]
    t_D_val = tD;    % D-cell skin thickness [m]
    tau_val = Tau;      % Applied torque [Nm]
    
    q_wb_num = subs(q_wb_sol, ...
        {A_wb, A_D, c_wb, t_U, t_L, h_wb, t_F, t_R, b_D, t_D, tau}, ...
        {A_wb_val, A_D_val, c_wb_val, t_U_val, t_L_val, h_wb_val, ...
         t_F_val, t_R_val, b_D_val, t_D_val, tau_val});
    
    q_D_num = subs(q_D_sol, ...
        {A_wb, A_D, c_wb, t_U, t_L, h_wb, t_F, t_R, b_D, t_D, tau}, ...
        {A_wb_val, A_D_val, c_wb_val, t_U_val, t_L_val, h_wb_val, ...
         t_F_val, t_R_val, b_D_val, t_D_val, tau_val});
end