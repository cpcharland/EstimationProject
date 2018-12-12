function dxdt = RHS(t,x,u_,w)
    % Params
    Gamma = [zeros(1,3);1,0,0;zeros(1,3);0,1,0;zeros(1,3);0,0,1];
    rho = 0.02;
    g = 3.711;
    beta = pi/4;
    C_D = 0.2;
    m_f = 390;
    w_f = 1;
    h_f = 0.5;
    d_f = 1;
    m_b = 1510;
    w_b = 3.2;
    h_b = 2.5;
    d_b = 2.9;
    h_cm = .9421;
    w_cm = w_b/2;
    A_side = h_b*d_b+h_f*d_f;
    A_bot = w_b*d_b+w_f*d_f;
    
    T_1 = u_(1);
    T_2 = u_(2);
    
    % Alpha pre-processing
    theta = x(5);
    alpha = atan(x(4)/x(2));
    if isnan(alpha)
        alpha  = 0;
    end
    alpha = -abs(alpha);
    
    % Drag forces and Moment of Inertia
    F_D_xi = 1/2*C_D*rho*(A_side*cos(theta-alpha)+A_bot*sin(theta-alpha))*x(2)*sqrt(x(2)^2+x(4)^2);
    F_D_z = 1/2*C_D*rho*(A_side*cos(theta-alpha)+A_bot*sin(theta-alpha))*x(4)*sqrt(x(2)^2+x(4)^2);
    I_n = 1/12*(m_b*(w_b^2+h_b^2)+m_f*(w_f^2+h_f^2));
    
    % d/dt vector from nonlinear dyamics
    dxdt(1) = x(2);
    dxdt(2) = 1/(m_b+m_f)*(T_1*(cos(beta)*sin(theta)+sin(beta)*cos(theta))+T_2*(cos(beta)*sin(theta)-sin(beta)*cos(theta))-F_D_xi);
    dxdt(3) = x(4);
    dxdt(4) = 1/(m_b+m_f)*(T_1*(cos(beta)*cos(theta)-sin(beta)*sin(theta))+T_2*(cos(beta)*cos(theta)+sin(beta)*sin(theta))-F_D_z)-g;
    dxdt(5) = x(6);
    dxdt(6) = 1/I_n*((T_1-T_2)*cos(beta)*w_b/2+(T_2-T_1)*sin(beta)*h_cm);
    dxdt = dxdt'+Gamma*w;
end

