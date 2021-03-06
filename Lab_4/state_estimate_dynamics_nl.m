function xdot_gen = state_estimate_dynamics_nl(t,x_gen,parameters)
    % extracting parameters
    m = parameters.m;
    M = parameters.M;
    g = parameters.g;
    l = parameters.l;
    A = parameters.A;
    B = parameters.B;
    C = parameters.C;
    K = parameters.K;
    L = parameters.L;
    
    x = x_gen(1:4);
    xhat = x_gen(5:end);

    % parameters used for the nonlinear model
    u = K*x;
    z = x(1);
    z_dot = x(2);
    theta = x(3);
    theta_dot = x(4);
    
    % obtain nonlinear derivatives
    z_ddot=(-m*l*sin(theta)*theta_dot^2 + m*g*sin(theta)*cos(theta)+ u)/...
            (M + m*sin(theta)^2);
    theta_ddot=(-m*l*sin(theta)*cos(theta)*theta_dot^2 + ...
            (M+m)*g*sin(theta) + u*cos(theta))/(l*(M+m*sin(theta)^2));
    xdot = [z_dot; z_ddot; theta_dot; theta_ddot];
    
    y = C*x;
    xhatdot = (A + L*C)*xhat + B*K*x - L*y;
    
    xdot_gen = [xdot; xhatdot];
end