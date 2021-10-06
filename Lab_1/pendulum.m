function xdot = pendulum(t,x,parameters)

% extracting parameters
M = parameters.M;
g = parameters.g;
l = parameters.l;

x1 = x(1);
x2 = x(2);

%zero input
u = 0;

%compute xdot using nonlinear equation
xdot = [x2; -g/l*sin(x1) - 1/(M*l)*cos(x1)*u];
end
