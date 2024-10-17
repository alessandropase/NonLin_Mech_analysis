function basins = basins_attractions(sys, dt, t_end, omega, x_0_vect)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the amplitude of the motion at regime for a range of initial
% conditions
% sys      = mechanical system
% dt       = timestep integration
% t_end    = time integration window
% omega    = frequency of interest
% x_0_vect = initial condition span - equal for both x_0_1 and x_0_2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(x_0_vect);

basins = zeros(n);

for i = 1 : n
    for j = 1 : n
        sys.x_0 = [x_0_vect(i), x_0_vect(j)]';
        sys.x_d_0 = [0, 0]';
        y_0 = [sys.x_0; sys.x_d_0];
        [~, x] = time_integration(y_0, sys, dt, t_end, omega);
        x = x';
        A_1 = max(x(1, end-1000:end));
        basins(i, j) = A_1;
        clc;
        perc = ((i-1)*n+j)/(n^2)*100;
        fprintf('Percentage of steps perfomed: %.2f%%\n', perc)
    end
end