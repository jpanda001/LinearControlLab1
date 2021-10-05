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
    
    %setting state-space matrices
    A = [0 1; -g/l 0];
    B = [0; -1/(M*l)];
    
    %setting state variables
    x_cur = [x1; x2];
    
    %compute xdot = Ax + Bu
    xdot = [A*x_cur + B*u; z_dot;];
end
