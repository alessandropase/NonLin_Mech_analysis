function [x_0_sol, xdot_0_sol, exitflag] = shooting(sys, T, dt, omega, y_guess)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function provide the ICs for which the motion is periodic
% sys = mechanical system - type: struct
% T   = time period
% dt  = timestep for integration
% y_guess = initial guess for the initial condition (nearby of the y_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_end = T;
nDof = length(sys.M);

fun = @(y) ff_fun(y, sys, dt, t_end, omega);

opt = optimoptions('fsolve','StepTolerance', 1e-10, 'functiontolerance', 1e-5, 'MaxFunctionEvaluations', 3e4, 'MaxIterations', 5e3, 'Display', 'off');
[y_sol, ~, exitflag] = fsolve(fun, y_guess, opt);

%y_sol = fsolve(fun, y_guess);
%fprintf('flag = %2.4d\n', exitflag);
x_0_sol = y_sol(1:nDof);
xdot_0_sol = y_sol(nDof+1:end);