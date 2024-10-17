function [A_vect, omega] = continuation(sys, dt, omega_span)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This algorithm provides the continuation for the mechanical system
% sys        = mechanical system - type: struct
% dt         = timestep for numerical integration
% omega_span = continuation span - is only the first and points are
%                                  specified the algorithm will create a
%                                  vector with 101 points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% User choice
fprintf("Select the output of the sequential continuation algorithm: \n ")
fprintf("1. NLFRs \n ")
fprintf("2. NNMs \n ")

type = input('You are choosing: ');
if type == 2 
    sys.C = zeros(length(sys.M));
    sys.f_ext = @ (t, omega) zeros(length(sys.M), 1);
end

% Definition of the frequency range
omega_0 = omega_span(1);
omega_end = omega_span(end);
if length(omega_span)>2
    omega = omega_span;
else
    omega = linspace(omega_0, omega_end, 101);
end

% Initialization of the variables
A_vect_1 = zeros(length(sys.M), length(omega));
A_vect_2 = zeros(length(sys.M), length(omega));
j = 1;
y_0 = [sys.x_0; sys.x_d_0];
guard = 0;
omega_1 = omega;
% Sequential continuation
for w = omega
    clc;
    fprintf('Frequency number: %d / %d \n', j, length(omega));
    perc = (j+1)/length(omega)*100;
    fprintf('Percentage of the total steps: %.2f%%\n', perc);
    T = 2*pi/w;
    t_end = T;
    [x_0_sol, xdot_0_sol, exitflag] = shooting(sys, T, dt, w, y_0);
    if exitflag < 1
        omega_1 = omega(1:j-1);
        A_vect_1 = A_vect_1(:, 1:j-1);
        guard = 1;
        break
    end
    y_0 = [x_0_sol; xdot_0_sol];
    [~, x] = time_integration(y_0, sys, dt, t_end, w);
    for i = 1 : length(sys.M)
        A_vect_1(i, j) = max(abs(x(:, i)));
    end
    j = j + 1;    
end


if type == 2
    omega = omega_1;
    A_vect = A_vect_1;
elseif guard == 0
    fprintf('No turning point \n');
    omega = omega_1;
    A_vect = A_vect_1;
else
    omega = flip(omega);
    j = 1;
    fprintf('Turning point is reached: the algorithm will re-start from the last frequency \n')
    omega_2 = omega;
    for w = omega
        fprintf('Frequency number: %d \n', j);
        T = 2*pi/w;
        t_end = T;
        sys.x_0 = 1e-4*[1, -1]';
        [x_0_sol, xdot_0_sol, exitflag] = shooting(sys, T, dt, w, y_0);
        
        if exitflag ~= 1
            omega_2 = omega(1:j-1);
            A_vect_2 = A_vect_2(:, 1:j-1);
            break
        end
        y_0 = [x_0_sol; xdot_0_sol];
        [~, x] = time_integration(y_0, sys, dt, t_end, w);
        for i = 1 : length(sys.M)
            A_vect_2(i, j) = max(abs(x(:, i)));
        end
        j = j + 1;    
    end
    
    omega = [omega_1, flip(omega_2, 2)];
    A_vect = [A_vect_1, flip(A_vect_2, 2)];
    clc;
    fprintf("Reversed continuation performed")
end
end
