clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%% Numerical Integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define the set of parameters specific to our ODE
parameters.M = 0.2;
parameters.l= 0.15;
parameters.g = 9.81;


% define numerical integration parameters: relativa and absolute tolerance
options = odeset('RelTol',1e-7,'AbsTol',1e-7);

% define ode solver parameters
Tspan = linspace(0,10,1e3);


% initial condition 1
x0 = [0; sqrt(parameters.g/parameters.l)];

% numerically solve ode using ode45
[t,x]=ode45(@pendulum,Tspan,x0,options,parameters);

% extracting the solution to the ODE
x1 = x(:,1);
x2 = x(:,2);

% plot x1 vs t
figure
subplot(211)
plot(t, x1)
xlabel("t")
ylabel("x1")

% plot x2 vs t
subplot(212)
plot(t, x2)
xlabel("t")
ylabel("x2")

sgtitle('initial condition 1: [0; sqrt(g/l)]') 

% plot orbit x2 vs x1
figure
plot(x1, x2)
xlabel("x1")
ylabel("x2")
title("initial condition 1: [0; sqrt(g/l)]")

% initial condition 2
x0 = [0; 1.99*sqrt(parameters.g/parameters.l)];

% numerically solve ode using ode45
[t,x]=ode45(@pendulum,Tspan,x0,options,parameters);

% extracting the solution to the ODE
x1 = x(:,1);
x2 = x(:,2);

% plot x1 vs t
figure
subplot(211)
plot(t, x1)
xlabel("t")
ylabel("x1")

% plot x2 vs t
subplot(212)
plot(t, x2)
xlabel("t")
ylabel("x2")

sgtitle('initial condition 2: [0; 1.99*sqrt(g/l)]') 

% plot orbit x2 vs x1
figure
plot(x1, x2)
xlabel("x1")
ylabel("x2")
title("initial condition 2: [0; 1.99*sqrt(g/l)]")

%%%%%%%%%%%%%%%%%%%%%%%% Symbolic Linearization %%%%%%%%%%%%%%%%%%%%%%%%%%%

% define symbolic variables
syms x1 x2 t u m l g theta_bar real

% define state vector as a composite of symbolic state variables
x = [x1; x2];

% define system dynamics
xdot = [x2; -g/l * sin(x1) - 1/(m*l) * cos(x1) * u];
y = x1;

% find jacobian
A_raw = jacobian(xdot,x);
B_raw = jacobian(xdot,u);
C_raw = jacobian(y,x);
D_raw = jacobian(y,u);


% equilibrium condition 1
x_bar_1 = [0; 0];
u_bar_1 = 0;

% linearize around equilibirum conditions 1
A = subs(subs(A_raw,x,x_bar_1), u, u_bar_1)
B = subs(subs(B_raw,x,x_bar_1), u, u_bar_1)
C = subs(subs(C_raw,x,x_bar_1), u, u_bar_1)
D = subs(subs(D_raw,x,x_bar_1), u, u_bar_1)


% equilibrium condition 2
x_bar_2 = [theta_bar; 0];
u_bar_2 = -m*g*tan(theta_bar);

% linearize around equilibirum conditions 2
A = subs(subs(A_raw,x,x_bar_2), u, u_bar_2)
B = subs(subs(B_raw,x,x_bar_2), u, u_bar_2)
C = subs(subs(C_raw,x,x_bar_2), u, u_bar_2)
D = subs(subs(D_raw,x,x_bar_2), u, u_bar_2)

%%%%%%%%%%%% From Symbolic Expression to Numerical Integration %%%%%%%%%%%%

% relinearize the system of equations for the appropriate equilibrium
% condition
x_bar = [0; 0];
u_bar = 0;

%linearize around equilibirum conditions 1
A = subs(subs(A_raw,x,x_bar), u, u_bar)
B = subs(subs(B_raw,x,x_bar), u, u_bar)
C = subs(subs(C_raw,x,x_bar), u, u_bar)
D = subs(subs(D_raw,x,x_bar), u, u_bar)

% assign create vector
syms z1 z2
z = [z1; z2];
zdot = A*z + B*u;

Xdot = [xdot; zdot];

%substitute numerical values to non-state-variables
g_val = 9.81;
m_val = 0.2;
l_val = 0.15;
u_val = 0;

Xdot = subs(Xdot, g, g_val);
Xdot = subs(Xdot, m, m_val);
Xdot = subs(Xdot, l, l_val);
Xdot = subs(Xdot, u, u_val);

% verify only state variables remain as variables
symvar(Xdot)

% create a Matlab function inline
augmented_pen = matlabFunction(Xdot,'Vars',{t,[x;z]})


% configure parameter of ode45, assume to be the same as in part 3
options = odeset('RelTol',1e-7,'AbsTol',1e-7);  
Tspan = linspace(0,10,1e3);

% initial condition 1 (note:set initial condition for all 4 variables)
x0 = [0; sqrt(g_val/l_val); 0; sqrt(g_val/l_val)];
[t,x]=ode45(augmented_pen,Tspan,x0,options);

% plot with 2 subplots
figure

subplot(211)
plot(t, x(:,1))
hold on
plot(t, x(:,3))
hold off
legend('x1', 'z1')
xlabel("t")
ylabel("\theta")

%plot x2 vs t
subplot(212)
plot(t, x(:,2))
hold on
plot(t, x(:,4))
hold off
legend('x2', 'z2')
xlabel("t")
ylabel('$\dot{\theta}$', 'Interpreter','latex')

sgtitle('initial condition 1: [0; sqrt(g/l)]') 

%plot orbit x2 vs x1 and z2 vs z1
figure
plot(x(:,1),  x(:,2))
hold on
plot(x(:,3),  x(:,4))
hold off
xlabel("\theta")
ylabel('$\dot{\theta}$', 'Interpreter','latex')
legend('x2 vs x1', 'z2 vs z1')
title("initial condition 1: [0; sqrt(g/l)]")


% initial condition 2 (note:set initial condition for all 4 variables)
x0 = [0; 1.99*sqrt(g_val/l_val); 0; 1.99*sqrt(g_val/l_val)];
[t,x]=ode45(augmented_pen,Tspan,x0,options);

% plot with 2 subplots
figure

subplot(211)
plot(t, x(:,1))
hold on
plot(t, x(:,3))
hold off
legend('x1', 'z1')
xlabel("t")
ylabel("\theta")

%plot x2 vs t
subplot(212)
plot(t, x(:,2))
hold on
plot(t, x(:,4))
hold off
legend('x2', 'z2')
xlabel("t")
ylabel('$\dot{\theta}$', 'Interpreter','latex')

sgtitle('initial condition 2: [0; 1.99*sqrt(g/l)]') 

%plot orbit x2 vs x1 and z2 vs z1
figure
plot(x(:,1),  x(:,2))
hold on
plot(x(:,3),  x(:,4))
hold off
xlabel("\theta")
ylabel('$\dot{\theta}$', 'Interpreter','latex')
legend('x2 vs x1', 'z2 vs z1')
title('initial condition 2: [0; 1.99*sqrt(g/l)]')

%%%%%%%%%%%%%%% LTI Representation %%%%%%%%%%%%%%%%
syms s

% Converting symbolic matrices to double matrices
doubleA = double(subs(A, {m, g, l}, {m_val, g_val, l_val}));
doubleB = double(subs(B, {m, g, l}, {m_val, g_val, l_val}));
doubleC = double(subs(C, {m, g, l}, {m_val, g_val, l_val}));
doubleD = double(subs(D, {m, g, l}, {m_val, g_val, l_val}));

% Creating an LTI object
sys = ss(doubleA, doubleB, doubleC, doubleD);

% Converting state representation to TF representation
% TF is the same as the one in page 58 of textbook
zpk(sys)

% Eigenvalues of A
eig(A)

%%%%%%%%%%%%% Pendulum Stabilization %%%%%%%%%%%%%%
close all
% Can be further optimized
K_d = 17;
K_p = 1.9;
K_i = 7;
n = 55;

% PD Compensator TF
C = -1*tf([K_p+K_d, K_p*n],[1, n]);
% C_2 = -0.5*tf([K_p + K_d, K_i + K_p*n, K_i*n],[1, n, 0]);
% controlSystemDesigner(sys);

%check root locus to ensure that locus does not have region on right hand
%plane

% check that there is no zero on the right hand plane
% zpk(1+C*sys)

% Creating bode plot of Controller TF to check how good it is
% figure
% bode(C)

% Finding Zeros and Poles of Controller TF
% [num, den] = tfdata(1 + C*sys)
% [z, p, k] = tf2zpk(cell2mat(num),cell2mat(den))

% Plotting zeros and poles of Controller TF
% fvtool(cell2mat(num),cell2mat(den),'polezero')
% text(real(z)-0.1,imag(z)-0.1,'\bfZeros','color',[0 0.4 0])
% text(real(p)-0.1,imag(p)-0.1,'\bfPoles','color',[0.6 0 0])

% Estimating to controller state-representation matrices 
[F, G, H, L] = ssdata(C)
% 
% % Solving with ODE45
parameters = struct('M', m_val, 'g', g_val, 'l', l_val, 'F', F, 'G', G, 'H', H, 'L', L);
% 
x0 = [0 1.99*sqrt(g_val/l_val) 0];
% 
Tspan = linspace(0,10,1e3);
[t,x] = ode45(@controlled_pendulum, Tspan, x0, options, parameters);
% 
% % Plotting x1 and x2 against time
x1 = x(:,1);
x2 = x(:,2);

figure
subplot(211)
plot(t, x1)
xlabel("t")
ylabel("x1")

subplot(212)
plot(t, x2)
xlabel("t")
ylabel("x2")