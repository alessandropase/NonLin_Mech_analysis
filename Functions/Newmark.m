function [tt, x] = Newmark(M, C, K, S, f_nl, f_ext, x_0, x_dot_0, dt, t_end, gamma, beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newmark algorithm: time integration of the mechanical system using the
% system and his tangent S
% This algorithm is not longer used in the project to its slow
% computational time wrt ode45
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_dotdot_0 = M\(f_ext(0)-f_nl(x_0, x_dot_0));

f = @(x, x_dot) C*x_dot + K*x + f_nl(x, x_dot);

i = 1;

nDofs = length(M);
n = ceil(t_end/dt);
x = zeros(nDofs, n+1);
x(:, 1) = x_0;
x_dot = zeros(nDofs, n+1);
x_dot(:, 1) = x_dot_0;
x_dotdot = zeros(nDofs, n+1);
x_dotdot(:, 1) = x_dotdot_0;
epsilon = 1e-4;

for tt = dt:dt:t_end
    i = i+1;
    x_dot(:, i) = x_dot(:, i-1) + (1-gamma)*dt*x_dotdot(:, i-1);
    x(:, i) = x(:, i-1) + dt*x_dot(:, i-1) + (0.5-beta)*dt^2*x_dotdot(:, i-1);
    x_dotdot(:, i) = 0;
    res = M*x_dotdot(:, i) + f(x(:, i), x_dot(:, i)) - f_ext(tt);
    while norm(res) > epsilon*norm(f(x(:, i), x_dot(:, i)))
        dx = S(x(:, i))\(-res);
        x(:, i) = x(:, i) + dx;
        x_dot(:, i) = x_dot(:, i) + (gamma)/(beta*dt)*dx;       
        x_dotdot(:, i) = x_dotdot(:, i) + 1/(beta*dt^2)*dx;
        res = M*x_dotdot(:, i) + f(x(:, i), x_dot(:, i)) - f_ext(tt);
    end
end

tt = 0:dt:t_end;