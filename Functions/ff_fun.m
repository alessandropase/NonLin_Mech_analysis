function h = ff_fun(y_0, sys, dt, t_end, omega)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% objective function for the fsolve in shooting algorithm
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
opt = odeset('RelTol',1e-12, 'AbsTol',1e-8);
[~, x] = ode45(odefun, tspan, y_0, opt);
%[tt, x] = ode45(odefun, tspan, y_0);

h = (x(end, :) - x(1, :)); 
