function xdot_gen = state_estimate_dynamics_linear(t,x_gen,parameters)

    % extracting parameters
    A = parameters.A;
    B = parameters.B;
    C = parameters.C;
    K = parameters.K;
    L = parameters.L;
    
    x = x_gen(1:4);
    xhat = x_gen(5:end);
    
    xdot = (A + B*K)*x;
    y = C*x;
    xhatdot = (A + L*C)*xhat + B*K*x - L*y;
    
    xdot_gen = [xdot; xhatdot];
end