clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms M m l g u y y_dot theta theta_dot t q1 q2 k1 k2 k3 k4

% define the system state
x = [y; y_dot; theta; theta_dot];

% define system dynamics
y_ddot=(-m*l*sin(theta)*theta_dot^2 + m*g*sin(theta)*cos(theta) + u)/...
    (M + m*sin(theta)^2);
theta_ddot=(-m*l*sin(theta)*cos(theta)*theta_dot^2 + (M+m)*g*sin(theta) +...
    u*cos(theta))/...
    (l*(M+m*sin(theta)^2));

xdot = [y_dot; y_ddot; theta_dot;theta_ddot]


% define the set of parameters specific to our ODE
parameters.M = 1.0731;
parameters.m = 0.2300;
parameters.l= 0.3302;
parameters.g = 9.8;

%linearization starts here

% equilibrium point
x_bar = [0; 0; 0; 0];
u_bar = 0;

% find jacobian for any point x
A_raw = jacobian(xdot,x);
B_raw = jacobian(xdot,u);

% substitute parameters as defined above
A_specific = subs(subs(subs(subs(A_raw,m,parameters.m),M,parameters.M), l, parameters.l), g, parameters.g);
B_specific = subs(subs(subs(subs(B_raw,m,parameters.m),M,parameters.M), l, parameters.l), g, parameters.g);

% find jacobian for equilibrium point
A = subs(subs(A_specific,x,x_bar), u, u_bar)
B = subs(subs(B_specific,x,x_bar), u, u_bar)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_double = double(A);
B_double = double(B);
K = [k1 k2 k3 k4]

C = (A_double+B_double*K)
[V,D] = eig(C)