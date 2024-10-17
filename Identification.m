%% Intro
clc;
clear;
close all;
addpath("C:\Users\aless\OneDrive\Desktop\[AERO0035] - Ghysens, Michiels, Pase\EXPERIMENTAL RES");
addpath("C:\Users\aless\OneDrive\Desktop\[AERO0035] - Ghysens, Michiels, Pase\NI2D RES");

% to load properly the files, please change the line above with your path

%% Natural frequencies of the linear system

M = [1, 0; 0, 1]; % mass matrix of the system
% C = [3, -1; -1, 3]; % linear damping matrix of the system
K = 10^4*[2, -1; -1, 2]; % linear stiffness matrix of the system

% compute the linear natural frequencies and mode shapes
[modes, omega_squared] = eig(K, M, 'vector');
freq = sqrt(omega_squared)/(2*pi);
[freq,ind] = sort(freq);
modes = modes(:,ind);

fprintf('Natural frequencies of the linear system\n');
fprintf('----------------------------------------\n');
fprintf('omega_1 = %f Hz\n', freq(1));
fprintf('omega_2 = %f Hz\n\n', freq(2));

%% Response of DOF 2 to a sine sweep of amplitude 50 N

load('group2_test3-1.mat');

sweep_rate = ext_force.fex{1}{5};                                           % sweep rate [Hz/min]
freq_i = ext_force.fex{1}{3};                                               % starting frequency [Hz]
r = sweep_rate / 60;                                                        % sweep rate [Hz/s]
freq = freq_i + r * res.t;                                                  % frequency [Hz]

response_plot = figure();
set(response_plot, 'defaulttextinterpreter', 'latex');
plot(freq, res.x(2,:));
xlabel('$f$ [Hz]');
ylabel('$x_2$ [m]');
set(gca, 'fontsize', 18, 'fontname', 'Times');
grid on;
grid minor;
xlim([5 35]);

%% Comparison between sweep up and sweep down for a forcing amplitude of 50 N

% Sweep up
load('group2_test3-1.mat');
sweep_rate_up = ext_force.fex{1}{5};                                        % sweep rate, [Hz/min]
freq_i_up = ext_force.fex{1}{3};                                            % starting frequency, [Hz]
r_up = sweep_rate_up / 60;                                                  % sweep rate, [Hz/s]
freq_up = freq_i_up + r_up * res.t;                                         % frequency, [Hz]
disp_up = res.x(2,:);                                                       % displacement at DOF 2, [m]

% Sweep down
load('group2_test4-1.mat');
sweep_rate_down = ext_force.fex{1}{5};                                      % sweep rate, [Hz/min]
freq_i_down = ext_force.fex{1}{3};                                          % starting frequency, [Hz]
r_down = sweep_rate_down / 60;                                              % sweep rate, [Hz/s]
freq_2 = freq_i_down + r_down * res.t;                                      % frequency, [Hz]
disp_down = res.x(2,:);                                                     % displacement at DOF 2, [m]

comparison_sweep_up_down = figure();
set(comparison_sweep_up_down, 'defaulttextinterpreter', 'latex');
plot(freq_up, disp_up, freq_2, disp_down);
xlabel('$f$ [Hz]');
ylabel('$x_2$ [m]');
set(gca, 'fontsize', 18, 'fontname', 'Times');
grid on;
grid minor;
legend('Sine sweep up excitation', 'Sine sweep down excitation');
xlim([5 35]);

%% Comparison between forcing amplitudes of 50 N and 30 N

% Forcing amplitude of 50 N
load('group2_test3-1.mat');
sweep_rate_50 = ext_force.fex{1}{5};                                        % sweep rate, [Hz/min]
freq_i_50 = ext_force.fex{1}{3};                                            % starting frequency, [Hz]
r_50 = sweep_rate_50 / 60;                                                  % sweep rate, [Hz/s]
freq_50 = freq_i_50 + r_50 * res.t;                                         % frequency, [Hz]
disp_50 = res.x(2,:);                                                       % displacement at DOF 2, [m]

% Forcing amplitude of 30 N
load('group2_test4-2.mat');
sweep_rate_30 = ext_force.fex{1}{5};                                        % sweep rate, [Hz/min]
freq_i_30 = ext_force.fex{1}{3};                                            % starting frequency, [Hz]
r_30 = sweep_rate_30 / 60;                                                  % sweep rate, [Hz/s]
freq_30 = freq_i_30 + r_30 * res.t;                                         % frequency, [Hz]
disp_30 = res.x(2,:);                                                       % displacement at DOF 2, [m]

comparison_50_30 = figure();
set(comparison_50_30, 'defaulttextinterpreter', 'latex');
plot(freq_50, disp_50, freq_30, disp_30);
xlabel('$f$ [Hz]');
ylabel('$x_2$ [m]');
set(gca, 'fontsize', 18, 'fontname', 'Times');
grid on;
grid minor;
legend('Excitation of 50 N', 'Excitation of 30 N');
xlim([5 35]);

%% Acceleration surface method
load('group2_test3-1.mat');
sweep_rate = ext_force.fex{1}{5}; % sweep rate, [Hz/min]
freq_i = ext_force.fex{1}{3}; % starting frequency, [Hz]
r = sweep_rate / 60; % sweep rate, [Hz/s]
freq = freq_i + r * res.t; % frequency, [Hz]

% select only the second resonance mode starting at 24 Hz
index = find(freq >= 24, 1);
freq = freq(index:end);

rel_disp = res.x(2,index:end) - res.x(1,index:end); % relative displacement x_2-x_1, [m]
rel_vel = res.xd(2,index:end) - res.xd(1,index:end); % relative velocity v_2-v_1, [m/s]
acceleration = res.xdd(2,index:end); % acceleration at DOF 2, [m/s^2]

% Acceleration surface
acceleration_surface = figure();
set(acceleration_surface, 'defaulttextinterpreter', 'latex');
plot3(rel_disp, rel_vel, -acceleration);
%title('Acceleration surface');
xlabel('$x_2 - x_1$ [m]');
ylabel('$\dot{x}_2 - \dot{x}_1$ [m/s]');
zlabel('$ - \ddot{x}_2 \ [\mathrm{m}/\mathrm{s}^2]$');
set(gca, 'fontsize', 18, 'fontname', 'Times');
grid on;
grid minor;

% Computation of the points of the stiffness curve
acc_s = zeros(length(freq),1); % aceleration for zero velocity, [m/s^2]
disp = zeros(length(freq),1); % displacement for zero velocity, [m]
n = 0; % number of points in the stiffness curve
for i = 1:(length(freq)-1)
    if rel_vel(i)*rel_vel(i+1) < 0 % velocity is 0 between both points
        n = n + 1;
        % linear interpolation of the acceleration and displacement at zero velocity
        acc_s(n) = interp1([rel_vel(i), rel_vel(i+1)], [acceleration(i), acceleration(i+1)], 0);
        disp(n) = interp1([rel_vel(i), rel_vel(i+1)], [rel_disp(i), rel_disp(i+1)], 0);
    end
end
acc_s = acc_s(1:n);
disp = disp(1:n);

% Acceleration surface and points with zero velocity
ASM_stiffness = figure();
set(ASM_stiffness, 'defaulttextinterpreter', 'latex');
plot3(rel_disp, rel_vel, -acceleration);
xlabel('$x_2 - x_1$ [m]');
ylabel('$\dot{x}_2 - \dot{x}_1$ [m/s]');
zlabel('$ - \ddot{x}_2 \ [\mathrm{m}/\mathrm{s}^2]$');
set(gca, 'fontsize', 18, 'fontname', 'Times');
grid on;
grid minor;
hold on;
scatter3(disp, zeros(length(disp),1), -acc_s, 'x', 'MarkerEdgeColor', [0.6350 0.0780 0.1840]);

% Stiffness curve
stiffness_curve = figure();
set(stiffness_curve, 'defaulttextinterpreter', 'latex');
scatter(disp, -acc_s, 'x', 'MarkerEdgeColor', [0.6350 0.0780 0.1840]);
xlabel('$x_2-x_1$ [m]');
ylabel('$ - \ddot{x}_2 \ [\mathrm{m}/\mathrm{s}^2]$');
set(gca, 'fontsize', 18, 'fontname', 'Times');
grid on;
grid minor;

% Computation of the points of the damping curve
acc_d = zeros(length(freq),1); % aceleration for zero displacement, [m/s^2]
vel = zeros(length(freq),1); % velocity for zero displacement, [m/s]
n = 0; % number of points in the damping curve
for i = 1:(length(freq)-1)
    if rel_disp(i)*rel_disp(i+1) < 0 % displacement is 0 between both points
        n = n + 1;
        % linear interpolation of the acceleration and velocity at zero displacement
        acc_d(n) = interp1([rel_disp(i), rel_disp(i+1)], [acceleration(i), acceleration(i+1)], 0);
        vel(n) = interp1([rel_disp(i), rel_disp(i+1)], [rel_vel(i), rel_vel(i+1)], 0);
    end
end
acc_d = acc_d(1:n);
vel = vel(1:n);

% Acceleration surface and points with zero displacement
ASM_damping = figure();
set(ASM_damping, 'defaulttextinterpreter', 'latex');
plot3(rel_disp, rel_vel, -acceleration);
xlabel('$x_2 - x_1$ [m]');
ylabel('$\dot{x}_2 - \dot{x}_1$ [m/s]');
zlabel('$ - \ddot{x}_2 \ [\mathrm{m}/\mathrm{s}^2]$');
set(gca, 'fontsize', 18, 'fontname', 'Times');
grid on;
grid minor;
hold on;
scatter3(zeros(length(vel),1), vel, -acc_d, 'x', 'MarkerEdgeColor', [0.6350 0.0780 0.1840]);

% Damping curve
damping_curve = figure();
set(damping_curve, 'defaulttextinterpreter', 'latex');
scatter(vel, -acc_d, 'x', 'MarkerEdgeColor', [0.6350 0.0780 0.1840]);
%title('Damping curve');
xlabel('$\dot{x}_2-\dot{x}_1$ [m/s]');
ylabel('$ - \ddot{x}_2 \ [\mathrm{m}/\mathrm{s}^2]$');
set(gca, 'fontsize', 18, 'fontname', 'Times');
grid on;
grid minor;

%% Restoring force surface method

M = [1, 0; 0, 1]; % mass matrix of the system
C = [3, -1; -1, 3]; % linear damping matrix of the system
K = 10^4*[2, -1; -1, 2]; % linear stiffness matrix of the system

load('group2_test3-1.mat');
sweep_rate = ext_force.fex{1}{5}; % sweep rate, [Hz/min]
freq_i = ext_force.fex{1}{3}; % starting frequency, [Hz]
r = sweep_rate / 60; % sweep rate, [Hz/s]
freq = freq_i + r * res.t; % frequency, [Hz]

% select only the second resonance mode from 24 to 32 Hz
i_i = find(freq >= 24, 1);
i_f = find(freq >= 32, 1);
freq = freq(i_i:i_f);

disp_1 = res.x(1,i_i:i_f); % displacement at DOF 1, [m]
disp_2 = res.x(2,i_i:i_f); % displacement at DOF 2, [m]
vel_1 = res.xd(1,i_i:i_f); % velocity at DOF 1, [m/s]
vel_2 = res.xd(2,i_i:i_f); % velocity at DOF 2, [m/s]
acc_1 = res.xdd(1,i_i:i_f); % acceleration at DOF 1, [m/s^2]
acc_2 = res.xdd(2,i_i:i_f); % acceleration at DOF 2, [m/s^2]
force_1 = res.pex(1,i_i:i_f); % force at DOF 1, [N]
force_2 = zeros(1, length(force_1)); % force at DOF 2, [N]

rel_disp = disp_2 - disp_1; % relative displacement between DOFs 2 and 1, [m]
rel_vel = vel_2 - vel_1; % relative velocity between DOFs 2 and 1, [m/s]

% Assumed functional forms of the nonlinear force
%f_1 = [-rel_disp; rel_disp]; % linear stiffness force
%f_2 = [rel_disp.^2; -rel_disp.^2]; % quadratic stiffness force
f_3 = [-rel_disp.^3; rel_disp.^3]; % cubic stiffness force
f_4 = [rel_disp.^4; -rel_disp.^4]; % stiffness force of degree 4
%f_5 = [-rel_disp.^5; rel_disp.^5]; % stiffness force of degree 5
%f_6 = [rel_disp.^6; -rel_disp.^6]; % stiffness force of degree 6
f_7 = [-rel_disp.^7; rel_disp.^7]; % stiffness force of degree 7

% Computation of the right-hand side matrix of the system to solve
f_ext = [force_1; force_2];
q = [disp_1; disp_2];
q_dot = [vel_1; vel_2];
q_dotdot = [acc_1; acc_2];
right_term = f_ext - M * q_dotdot - C * q_dot - K * q;
right_term = right_term(:);

% Computation of the left-hand side matrix of the system to solve
%f_1 = f_1(:);
%f_2 = f_2(:);
f_3 = f_3(:);
f_4 = f_4(:);
%f_5 = f_5(:);
%f_6 = f_6(:);
f_7 = f_7(:);
f_matrix = [f_3, f_4, f_7];

% Computation of the coefficients k_i using the pseudoinverse
coeff = pinv(f_matrix) * right_term;

% Computation of the estimated nonlinear force
f_nl = f_matrix * coeff;

right_term_2 = right_term(2:2:end); % measured nonlinear force at DOF 2, [N]
f_nl_2 = f_nl(2:2:end); % estimated nonlinear force at DOF 2, [N]

% Measured nonlinear force at zero velocity
f_nl_disp = zeros(length(freq),1); % nonlinear force for zero relative velocity, [N]
disp = zeros(length(freq),1); % relative displacement for zero relative velocity, [m]
n = 0;
for i = 1:(length(freq)-1)
    if rel_vel(i)*rel_vel(i+1) < 0
        n = n + 1;
        % linear interpolation of the nonlinear force and displacement at zero velocity
        f_nl_disp(n) = interp1([rel_vel(i), rel_vel(i+1)], [right_term_2(i), right_term_2(i+1)], 0);
        disp(n) = interp1([rel_vel(i), rel_vel(i+1)], [rel_disp(i), rel_disp(i+1)], 0);
    end
end
f_nl_disp = f_nl_disp(1:n);
disp = disp(1:n);

% sort the elements in the ascending order of the relative displacement
[disp, indices_real] = sort(disp);
f_nl_disp = f_nl_disp(indices_real);
[rel_disp_sorted, indices] = sort(rel_disp);
f_nl_2 = f_nl_2(indices);

% Nonlinear force as a function of relative displacement for zero relative velocity
f_nl_disp_plot = figure();
set(f_nl_disp_plot, 'defaulttextinterpreter', 'latex');
plot(disp, f_nl_disp, 'linewidth', 2);
hold on;
plot(rel_disp_sorted, f_nl_2, 'linewidth', 1.2);
xlabel('$x_2 - x_1$ [m]');
ylabel('$f_{nl}(x_2-x_1,0)$ [N]');
set(gca, 'fontsize', 18, 'fontname', 'Times');
legend('Measured', 'Estimated', 'Location', 'Southeast');
legend('boxoff');
grid on;
grid minor;

% Measured nonlinear force at zero displacement
f_nl_vel = zeros(length(freq),1); % nonlinear force for zero relative displacement, [N]
vel = zeros(length(freq),1); % relative velocity for zero relative displacement, [m]
n = 0;
for i = 1:(length(freq)-1)
    if rel_disp(i)*rel_disp(i+1) < 0
        n = n + 1;
        % linear interpolation of the nonlinear force and velocity at zero displacement
        f_nl_vel(n) = interp1([rel_disp(i), rel_disp(i+1)], [right_term_2(i), right_term_2(i+1)], 0);
        vel(n) = interp1([rel_disp(i), rel_disp(i+1)], [rel_vel(i), rel_vel(i+1)], 0);
    end
end
f_nl_vel = f_nl_vel(1:n);
vel = vel(1:n);

% sort the elements in the ascending order of the relative velocity
[vel, indices] = sort(vel);
f_nl_vel = f_nl_vel(indices);

f_nl_vel_plot = figure();
set(f_nl_vel_plot, 'defaulttextinterpreter', 'latex');
scatter(vel, f_nl_vel, 'x');
hold on;
plot(vel, zeros(length(vel),1), 'linewidth', 1.2);
xlabel('$\dot{x}_2 - \dot{x}_1$ [m/s]');
ylabel('$f_{nl}(0,\dot{x}_2 - \dot{x}_1)$ [N]');
set(gca, 'fontsize', 18, 'fontname', 'Times');
legend('Measured', 'Estimated', 'Location', 'North');
legend('boxoff');
grid on;
grid minor;

%% Comparison between experimental results and simulations for a sine excitation at 10 Hz and 30 N

load('group2_test5-3.mat');

comparison_1 = figure();
set(comparison_1, 'defaulttextinterpreter', 'latex');
plot(res.t, res.x(1,:), 'r--', 'linewidth', 1.8);
xlabel('$t$ [s]');
ylabel('$x_1$ [m]');
set(gca, 'fontsize', 18, 'fontname', 'Times');
grid on;
grid minor;
hold on;
load('sim_5_3_disp_dof_1.mat');
plot(x, y, 'Color', [0 0.4470 0.7410], 'linewidth', 0.8);
legend('Experimental', 'Simulation');
xlim([0 1]);

%% Comparison between experimental results and simulations for a sine excitation at 30 Hz and 30 N

load('group2_test5-2.mat');

comparison_1 = figure();
set(comparison_1, 'defaulttextinterpreter', 'latex');
plot(res.t, res.x(1,:), 'r--', 'linewidth', 1.8);
xlabel('$t$ [s]');
ylabel('$x_1$ [m]');
set(gca, 'fontsize', 18, 'fontname', 'Times');
grid on;
grid minor;
hold on;
load('sim_5_2_disp_dof_1.mat');
plot(x, y, 'Color', [0 0.4470 0.7410], 'linewidth', 0.8);
legend('Experimental', 'Simulation');
xlim([0 1]);

%% Continuous Wavelet Transform (CWT)
load('EXPERIMENTAL RES\group2_test2-3.mat');
sweep_rate = 10;
freq_i = 5;
r = sweep_rate/60;
freq = freq_i + r * res.t;

% Now we are going to extrapolate the second mode
for i = 1 : length(freq)
    if freq(i) > 22  
        i_22 = i;
        break;
    end
end
for i = length(freq) : -1: 1
    if freq(i) < 32 
        i_32 = i;
        break;
    end
end
freq_cut = freq(i_22 : i_32);
displ1 = res.x(1, i_22:i_32);
displ2 = res.x(2, i_22:i_32);
veloc1 = res.xd(1, i_22:i_32);
veloc2 = res.xd(2, i_22:i_32);
accel1 = res.xdd(1, i_22:i_32);
accel2 = res.xdd(2, i_22:i_32);

index = 1 : 1 : length(displ1);
signal = accel1(index);
scales = freq_cut(1) : 1 : freq_cut(end);
wavelet = 'amor';
t_cut =  res.t(i_22 : i_32);
t = t_cut(index);
dt = t(2)-t(1);
fs = 1/dt;

% Perform CWT
[coefficients, f] = cwt(signal, wavelet, fs);

% Plot using contourf
close all;
figure('name', 'contourf');
h = contourf(t, f, abs(coefficients), 20, 'edgecolor','none');
colormap('jet');
colorbar;
xlabel('Time [s]', 'Interpreter', 'latex', 'Fontsize', 14);
ylabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize',14);
set(gca,'yscale','log')
axis tight;
clim([0, 800]);
ylim([1, 100]);

% Plot using pcolor
figure('name', 'Pcolor');
h = pcolor(t, f, abs(coefficients));
set(h, 'EdgeColor', 'none');
colormap('jet');
colorbar;
xlabel('Time [s]', 'Interpreter', 'latex', 'Fontsize', 14);
ylabel('Frequency [Hz]', 'Interpreter', 'latex', 'Fontsize', 14);
set(gca,'yscale','log')
axis tight;
clim([0, 800]);
ylim([1, 100]);
