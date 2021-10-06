function xdot = controlled_pendulum(t,x,parameters)

    % extracting parameters
    M = parameters.M;
    g = parameters.g;
    l = parameters.l;
    F = parameters.F;
    G = parameters.G;
    H = parameters.H;
    L = parameters.L;
    
    x1 = x(1);
    x2 = x(2);
    z = x(3);
    
    %zero input
    z_dot = F*z - G*x1;
    u = H*z - L*x1;

    x1_dot = x2;
    x2_dot = -g/l*sin(x1) - 1/(M*l)*cos(x1)*u;
    
    %compute xdot = Ax + Bu
    xdot = [x1_dot; x2_dot; z_dot;];
end
