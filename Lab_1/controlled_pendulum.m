function xdot = controlled_pendulum(t,x,parameters)

    % extracting parameters
    M = parameters.M;
    g = parameters.g;
    l = parameters.l;
    F = parameters.F;
    G = parameters.G;
    H = parameters.H;
    L = parameters.L;
    
    % Defining x1, x2 and z value
    x1 = x(1);
    x2 = x(2);
    z = x(3);
    
    % Solving for controller states and extacting the u value
    z_dot = F*z - G*x1;
    u = H*z - L*x1;

    % Solving for pendulum states with control input u
    x1_dot = x2;
    x2_dot = -g/l*sin(x1) - 1/(M*l)*cos(x1)*u;
    
    % concatenating all pendulum and controller states in xdot vector
    xdot = [x1_dot; x2_dot; z_dot;];
end
