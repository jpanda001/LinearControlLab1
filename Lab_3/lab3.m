clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms M m l g u y y_dot theta theta_dot t

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
poles_1 = [-1, -2, -3, -4];
poles_2 = [-1, -2, -3, -20];
negK_1 = place(A_double, B_double,poles_1);
negK_2 = place(A_double, B_double,poles_2);
K_1 = -negK_1               % to convert to the convention used in class
K_2 = -negK_2               % to convert to the convention used in class

%sanity check of closed loop poles_1 -> all poles_1 are at desired locations: [-1, -2, -3, -4]
A_cl_1 = A+B*K_1;           % close-loop system
double(eig(A_cl_1))

%sanity check of closed loop poles_2 -> all poles_2 are at desired locations: [-1, -2, -3, -20]
A_cl_2 = A+B*K_2;           % close-loop system
double(eig(A_cl_2))

% create the linearized close-loop system
Xdot_linearized_controlled_1 = A_cl_1 * x;
inverted_pen_1 = matlabFunction(Xdot_linearized_controlled_1,'Vars',{t, x});

Xdot_linearized_controlled_2 = A_cl_2 * x;
inverted_pen_2 = matlabFunction(Xdot_linearized_controlled_2,'Vars',{t, x});


options = odeset('RelTol',1e-7,'AbsTol',1e-7);  
Tspan = linspace(0,10,1e3);

x0 = [-0.5; 0; -pi/4; 0];
[t_val_1,x_val_1]=ode45(inverted_pen_1,Tspan,x0,options);
[t_val_2,x_val_2]=ode45(inverted_pen_2,Tspan,x0,options);


%plotting figures
subplot(3,2,1)              %plotting y
plot(t_val_1, x_val_1(:,1))
hold on                     
plot(t_val_2, x_val_2(:,1))
legend('[-1,-2,-3,-4]', '[-1,-2,-3,-20]')
xlabel("t (sec)")
ylabel("y (m)")
subplot(3,2,2)              %plotting y_dot
plot(t_val_1, x_val_1(:,2))
hold on
plot(t_val_2, x_val_2(:,2))
legend('[-1,-2,-3,-4]', '[-1,-2,-3,-20]')
xlabel("t (sec)")
ylabel("$\dot{y}$ (m/s)", 'Interpreter','latex')
subplot(3,2,3)              %plotting theta
plot(t_val_1, x_val_1(:,3))
hold on
plot(t_val_2, x_val_2(:,3))
legend('[-1,-2,-3,-4]', '[-1,-2,-3,-20]')
xlabel("t (sec)")
ylabel("\theta (rad)")
subplot(3,2,4)              %plotting theta_dot
plot(t_val_1, x_val_1(:,4))
hold on
plot(t_val_2, x_val_2(:,4))
legend('[-1,-2,-3,-4]', '[-1,-2,-3,-20]')
xlabel("t (sec)")
ylabel("$\dot{\theta}$ (rad/s)", 'Interpreter','latex')
subplot(3,2,5)              %plotting control input u
u_1 = -negK_1*x_val_1';                 % input u = K_1x; convert rows of x_val_1 to become columns of x's
u_2 = -negK_2*x_val_2';                 % input u = K_1x; convert rows of x_val_1 to become columns of x's
u_n_by_1_1 = u_1';                       % convert u to Nx1 as requested in lab manual
u_n_by_1_2 = u_2';                       % convert u to Nx1 as requested in lab manual

plot(t_val_1, u_n_by_1_1)   
hold on
plot(t_val_1, u_n_by_1_2)
legend('[-1,-2,-3,-4]', '[-1,-2,-3,-20]')
xlabel("t (sec)")
ylabel("conroller input (u)")



