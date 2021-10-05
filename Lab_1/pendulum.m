function xdot = pendulum(t,x,parameters)

% extracting parameters
M = parameters.M;
g = parameters.g;
l = parameters.l;

x1 = x(1);
x2 = x(2);

%zero input
u = 0;

%setting state-space matrices
A = [0 1; -g/l 0];
B = [0; -1/(M*l)];

%setting state variables
x_cur = [x1; x2];

%compute xdot = Ax + Bu
xdot = A*x_cur + B * u;
end
