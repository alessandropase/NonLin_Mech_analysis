function  [tt, x] = time_integration(y_0, sys, dt, t_end, omega)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to time-integrate the equation of motion using the
% ODE45 algorithm.
% y_0   = initial conditions
% sys   = mechanical system - type: struct
% dt    = time interval for the time integration
% t_end = time span of integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nDof = length(sys.M);
L = [zeros(nDof), eye(nDof);
    -inv(sys.M)*sys.K, -inv(sys.M)*sys.C];

g_nl =@(x, x_dot) [zeros(nDof, 1);
                  sys.M\sys.f_nl(x, x_dot)];
g_ext = @(t) [zeros(nDof, 1);
               sys.M\sys.f_ext(t, omega)];

%tspan= 0:dt:t_end;

tspan = [0 t_end];

odefun = @(t, y) L*y - g_nl(y(1:nDof), y(nDof+1:end)) + g_ext(t);
opt = odeset('RelTol',1e-7, 'AbsTol',1e-5);
[tt, y_sol] = ode45(odefun, tspan, y_0, opt);
%[tt, y_sol] = ode45(odefun, tspan, y_0);
x = y_sol(:, 1:nDof);