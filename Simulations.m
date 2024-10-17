%% Intro
clc;
clear;
close all;
addpath('C:\Users\aless\OneDrive\Desktop\[AERO0035] - Ghysens, Michiels, Pase\MATLAB Functions')

% to load properly the files, please change the line above with your path

% all the simulation obtained with the continuation algorithm performed in
% during the project are stored in the folder 'CONTINUATION RES' 

%% System definition
load("M_C_K.mat");                                                                              % Linear matrix definitions
sys.M = M;
sys.K = K;
sys.C = C;

sys.f_nl = @(x, x_dot) [-2.1e+10*(x(2)-x(1))^7+5.3e+07*(x(2)-x(1))^4-2.5e+06*(x(2)-x(1))^3;     % Nonlinearities definitions
                        +2.1e+10*(x(2)-x(1))^7-5.3e+07*(x(2)-x(1))^4+2.5e+06*(x(2)-x(1))^3];

sys.f_ext = @(t, omega) [50*sin(omega*t); 0];                                                   % Force definitions

dt = 5e-5;
sys.x_0 = 1e-4*[1, 1]';                                                                         % ICs definitions
sys.x_d_0 = [0, 0]';
sys.x_dd_0 = [0, 0]';

%% Shooting algorithm test - periodic motion computation
sys.x_0 = [0.002, -0.001]';
sys.x_d_0 = [0.002, -0.001]';
y_guess = [sys.x_0; sys.x_d_0];
y_0 = y_guess;
f = 20;
omega = 2*pi*f;
T = 1/f;
t_end = T;
close all;
figure;
[tt, x] = time_integration(y_0, sys, dt, t_end, omega);
a = plot(tt, x(:, 1), 'Linewidth', 1.5);

[sys.x_0, sys.x_d_0] = shooting(sys, T, dt, omega, y_guess);
y_0 = [sys.x_0; sys.x_d_0];
t_end = T;
hold on;
[tt, x] = time_integration(y_0, sys, dt, t_end, omega);

b = plot(tt, x(:, 1),'Linewidth', 1.5);

grid on;
grid minor;
xlabel('Time [s]', 'FontSize', 14, 'Interpreter','latex')
ylabel('Displacement [m]', 'FontSize', 14, 'Interpreter','latex')
legend([a, b], 'Random I.C.', 'I.C. from shooting', 'Interpreter','latex', 'Fontsize', 12, 'Location', 'best')

plot([0, T], [sys.x_0(1), sys.x_0(1)], 'k--', 'LineWidth', 1.1)

%% Sequential continuation for linear system
sys.f_nl = @(x, x_dot) [0; 0];
freq_span = linspace(5, 45, 801); 
omega_span = 2*pi.*freq_span;
clc;
tic;
[A_vect, omega] = continuation(sys, dt, omega_span);
time_cont = toc;
fprintf("\nThe continuation took %f seconds \n", time_cont);

% Plot for linear system
close all;
[modes, omega_squared] = eig(K, M, 'vector');
freq = sqrt(omega_squared)/(2*pi);
[freq,ind] = sort(freq);
figure;
plot(omega/2/pi, A_vect(1, :),'LineWidth',1.35)
hold on;
plot(omega/2/pi, A_vect(2, :),'LineWidth',1.35)
xlim([5 35]);
grid on;
grid minor;
xx = [freq(1), freq(1)];
yy = [0, max(A_vect(1, :))];
plot(xx, yy, 'k--', 'LineWidth',1.15);
xx = [freq(2), freq(2)];
plot(xx, yy, 'k--', 'LineWidth',1.15);
xlabel("Frequency [Hz]", 'Interpreter','latex', 'FontSize',14)
ylabel("Amplitude [m]", 'Interpreter','latex','FontSize',14)
legend("DoF 1", "DoF 2", 'Interpreter', 'latex', 'FontSize', 12)


%% Continuation for FRC at 3 different forcing amplitude
freq_span = linspace(5, 40, 801); 
omega_span = 2*pi.*freq_span;
clc;
sys.f_ext = @(t, omega) [10*sin(omega*t); 0];
tic;
[A_vect_10, omega_10] = continuation(sys, dt, omega_span);
time_cont = toc;
fprintf("\nThe continuation took %f seconds \n", time_cont);

sys.f_ext = @(t, omega) [30*sin(omega*t); 0];
tic;
[A_vect_30, omega_30] = continuation(sys, dt, omega_span);
time_cont = toc;
fprintf("\nThe continuation took %f seconds \n", time_cont);

sys.f_ext = @(t, omega) [50*sin(omega*t); 0];
tic;
[A_vect_50, omega_50] = continuation(sys, dt, omega_span);
time_cont = toc;
fprintf("\nThe continuation took %f seconds \n", time_cont);


%% Backbone curve computation
sys.x_0 = [0.01, -0.01]';
freq_span = linspace(27.57, 30.64, 501); 
omega_span = 2*pi.*freq_span;
[A_vect, omega] = continuation(sys, dt, omega_span);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from now on until the basins of attractions parts are the plots for the
% reports - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot of the three FRCs
close all;
load("CONTINUATION RES\v2_MAT_3frcs.mat");
figure('name','DOF 1');
a  = plot(omega_10/2/pi, A_vect_10(1, :),'ro-', 'LineWidth', 1, 'MarkerSize',1.5 );
grid on;
hold on;
b = plot(omega_30(1:569)/2/pi, A_vect_30(1, 1:569),'bo-', 'LineWidth', 1, 'MarkerSize', 1.5);
plot(omega_30(570:end)/2/pi, A_vect_30(1, 570:end),'bo-', 'LineWidth', 1, 'MarkerSize', 1.5);

d = plot(omega_50(1:583)/2/pi, A_vect_50(1, 1:583), 'o-', 'Color', "#7E2F8E", 'LineWidth', 1, 'MarkerSize', 1.5);
plot(omega_50(584:end)/2/pi, A_vect_50(1, 584:end),'o-', 'Color', "#7E2F8E", 'LineWidth', 1, 'MarkerSize',1.5);

xx = [15.91, 15.91];
yy = [0, max(A_vect_50(1, :))];
plot(xx, yy, 'k--', 'LineWidth', 1.2);

load("NI2D RES\10N_NLFR_DOF1.mat")
%plot(x, y, 'k', 'LineWidth',1.2)

legend([a, b, d], '$|F| = 10$ N', '$|F| = 30$ N', '$|F| = 50$ N', 'Interpreter', 'latex', 'Fontsize', 12, 'Location', 'best')
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'Fontsize', 14)
ylabel('Amplitude [m]', 'Interpreter', 'latex', 'Fontsize', 14)
xlim([5 40])

% DOF 2
figure('name','DOF 2');
a  = plot(omega_10/2/pi, A_vect_10(2, :),'ro-', 'LineWidth', 1, 'MarkerSize',1.5 );
grid on;
grid minor;
hold on;
b = plot(omega_30(1:569)/2/pi, A_vect_30(2, 1:569),'bo-', 'LineWidth', 1, 'MarkerSize', 1.5);
plot(omega_30(570:end)/2/pi, A_vect_30(2, 570:end),'bo-', 'LineWidth', 1, 'MarkerSize', 1.5);

d = plot(omega_50(1:583)/2/pi, A_vect_50(2, 1:583), 'o-', 'Color', "#7E2F8E", 'LineWidth', 1, 'MarkerSize', 1.5);
plot(omega_50(584:end)/2/pi, A_vect_50(2, 584:end),'o-', 'Color', "#7E2F8E", 'LineWidth', 1, 'MarkerSize',1.5);

xx = [15.91, 15.91];
yy = [0, max(A_vect_50(1, :))];
plot(xx, yy, 'k--', 'LineWidth', 1.2);

legend([a, b, d], '$|F| = 10$ N', '$|F| = 30$ N', '$|F| = 50$ N', 'Interpreter', 'latex', 'Fontsize', 12, 'Location', 'best')
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'Fontsize', 14)
ylabel('Amplitude [m]', 'Interpreter', 'latex', 'Fontsize', 14)
xlim([5 40])

%% Validation with Sinesweep
close all;
load("EXPERIMENTAL RES\group2_test4-2.mat")
load("CONTINUATION RES\v2_MAT_3frcs.mat");
sweep_rate_up = ext_force.fex{1}{5};                                        
freq_i_up = ext_force.fex{1}{3};                                            
r_up = sweep_rate_up / 60;                                                  
freq_up = freq_i_up + r_up * res.t; 
a = plot(freq_up, res.x(1, :), 'Color', "#A2142F", 'Linewidth', 1.5);
hold on;
load("EXPERIMENTAL RES\group2_test4-3.mat")
sweep_rate_down = ext_force.fex{1}{5};                                      
freq_i_down = ext_force.fex{1}{3};                                          
r_down = sweep_rate_down / 60;                                              
freq_2 = freq_i_down + r_down * res.t;                                                                   
b = plot(freq_2, res.x(1, :), 'Color', "#EDB120",'Linewidth', 1.5);
xlim([10 35]);
c = plot(omega_30(1:569)/2/pi, A_vect_30(1, 1:569), 'b-', 'LineWidth', 2.5);
plot(omega_30(570:end)/2/pi, A_vect_30(1, 570:end),'b-', 'LineWidth', 2.5);

xlim([5 35]);
grid on;
grid minor;
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'Fontsize', 14)
ylabel('Amplitude [m]', 'Interpreter', 'latex', 'Fontsize', 14)
legend([a, b, c], 'Sinesweep up - time response', 'Sinesweep down - time response', 'NLFR', 'Interpreter', 'latex', 'Fontsize', 12, 'Location', 'best')


%% Comparison between our continuation and NI2D harmonic balance continuation
close all;
load("CONTINUATION RES\v2_MAT_3frcs.mat");
figure('name','DOF 1');
a  = plot(omega_10/2/pi, A_vect_10(1, :),'ro', 'MarkerSize',4.5 );
grid on;
grid minor;
hold on;

d = plot(omega_50(1:583)/2/pi, A_vect_50(1, 1:583), 'o', 'Color', "blu", 'MarkerSize', 4.5);
plot(omega_50(584:end)/2/pi, A_vect_50(1, 584:end),'o', 'Color', "blu",  'MarkerSize',4.5);

xx = [15.91, 15.91];
yy = [0, max(A_vect_50(1, :))];
plot(xx, yy, 'k--', 'LineWidth', 1.2);

load("NI2D RES\10N_NLFR_DOF1.mat")
k = plot(x, y, 'k', 'LineWidth',1.2);
load("NI2D RES\50N_NLFR_DOF1.mat")
plot(x, y, 'k', 'LineWidth',1.2);

legend([a, d, k], '$|F| = 10$ N', '$|F| = 50$ N', 'NI2D Harmonic Balance Continuation','Interpreter', 'latex', 'Fontsize', 12, 'Location', 'best')
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'Fontsize', 14)
ylabel('Amplitude [m]', 'Interpreter', 'latex', 'Fontsize', 14)
xlim([5 40])

% DOF 2
figure('name','DOF 2');
a  = plot(omega_10/2/pi, A_vect_10(2, :),'ro', 'LineWidth', 1, 'MarkerSize',4.5 );
grid on;
grid minor;
hold on;

d = plot(omega_50(1:583)/2/pi, A_vect_50(2, 1:583), 'o', 'Color', "blu", 'LineWidth', 1, 'MarkerSize', 4.5);
e = plot(omega_50(584:end)/2/pi, A_vect_50(2, 584:end),'o', 'Color', "blu", 'LineWidth', 1, 'MarkerSize',4.5);

xx = [15.91, 15.91];
yy = [0, max(A_vect_50(1, :))];
plot(xx, yy, 'k--', 'LineWidth', 1.2);

load("NI2D RES\10N_NLFR_DOF2.mat")
k = plot(x, y, 'k', 'LineWidth',1.2);
load("NI2D RES\50N_NLFR_DOF2.mat")
plot(x, y, 'k', 'LineWidth',1.2);

legend([a, d, k], '$|F| = 10$ N', '$|F| = 50$ N', 'NI2D Harmonic Balance Continuation','Interpreter', 'latex', 'Fontsize', 12, 'Location', 'best')
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'Fontsize', 14)
ylabel('Amplitude [m]', 'Interpreter', 'latex', 'Fontsize', 14)
xlim([5 40])

%% Comparison between our FEP and NI2D
load("CONTINUATION RES\Backbone_501pts.mat")
E = @(A, B) 1e4*A^2/2 + 1e4*A^2/2 + (1e4*(B-A)^2)/2 + (2.1e+10*(B-A)^8)/8 -(5.3e+07*(B-A)^5)/5 +(2.5e+06*(B-A)^4)/4;
V = zeros(size(omega));
for index = 1 : length(omega)
    A = A_vect(1, index);
    B = -A_vect(2, index);
    V(index) = E(A, B);
end
close all;
plot(log10(V), omega/2/pi, 'o', 'LineWidth',1.1);
hold on;
load("NI2D RES\FEP.mat")
plot(x, y, 'r-', 'LineWidth', 1.4);
xlim([-2 3]);
ylim([27.3 31]);
grid on;
legend('Sequential continuation', 'NI2D Harmonic continuation','Interpreter', 'latex', 'Fontsize', 12, 'Location', 'best')
xlabel('$\log_{10}(Energy)$ [J]', 'Interpreter', 'latex', 'Fontsize', 14)
ylabel('Frequency [Hz]', 'Interpreter', 'latex', 'Fontsize', 14)


%% Comparison between our backbone and NI2D backbone
close all;
load("CONTINUATION RES\Backbone_501pts.mat");
figure('name','DOF 1');
a = plot(omega/2/pi, A_vect(1, :),'o', 'MarkerSize', 3.5);
grid on;
hold on;
xx = [27.56, 27.56];
yy = [0, 0.04];
b = plot(xx, yy, 'k--', 'Linewidth', 1.2);
load("NI2D RES\Backbone DOF1.mat")
c = plot(x, y, '-', 'Color', "#7E2F8E", 'LineWidth', 1.2);
ylim([0 0.04])
xlim([27.2 31])
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'Fontsize', 14)
ylabel('Amplitude [m]', 'Interpreter', 'latex', 'Fontsize', 14)
legend([a c b],'Sequential continuation', 'NI2D continuation', 'Linear resonance frequency', 'Location', 'best')

figure('name','DOF 2');
a = plot(omega/2/pi, A_vect(2, :),'o', 'MarkerSize', 3.5);
grid on;
hold on;
xx = [27.56, 27.56];
yy = [0, 0.04];
b = plot(xx, yy, 'k--', 'Linewidth', 1.2);
load("NI2D RES\Backbone DOF2.mat")
c = plot(x, y, '-', 'Color', "#7E2F8E", 'LineWidth', 1.2);
ylim([0 0.04])
xlim([27.2 31])
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'Fontsize', 14)
ylabel('Amplitude [m]', 'Interpreter', 'latex', 'Fontsize', 14)
legend([a c b],'Sequential continuation', 'NI2D continuation', 'Linear resonance frequency', 'Location', 'best')


%% Superimposition of the FRFs with backbone
close all;
load("CONTINUATION RES\v2_MAT_3frcs.mat");
load("CONTINUATION RES\Backbone_501pts.mat");
figure('name','DOF 1');
a  = plot(omega_10/2/pi, A_vect_10(1, :),'ro-', 'LineWidth', 1, 'MarkerSize',1.5 );
grid on;
hold on;
bb = plot(omega/2/pi, A_vect(1, :), 'k-', 'LineWidth', 1.5);
b = plot(omega_30(1:569)/2/pi, A_vect_30(1, 1:569),'bo-', 'LineWidth', 1, 'MarkerSize', 1.5);
plot(omega_30(570:end)/2/pi, A_vect_30(1, 570:end),'bo-', 'LineWidth', 1, 'MarkerSize', 1.5);

d = plot(omega_50(1:583)/2/pi, A_vect_50(1, 1:583), 'o-', 'Color', "#7E2F8E", 'LineWidth', 1, 'MarkerSize', 1.5);
plot(omega_50(584:end)/2/pi, A_vect_50(1, 584:end),'o-', 'Color', "#7E2F8E", 'LineWidth', 1, 'MarkerSize',1.5);

xx = [15.91, 15.91];
yy = [0, max(A_vect_50(1, :))];
plot(xx, yy, 'k--', 'LineWidth', 1.2);

legend([a, b, d, bb], '$|F| = 10$ N', '$|F| = 30$ N', '$|F| = 50$ N', 'Backbone curve', 'Interpreter', 'latex', 'Fontsize', 12, 'Location', 'best')
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'Fontsize', 14)
ylabel('Amplitude [m]', 'Interpreter', 'latex', 'Fontsize', 14)
xlim([25 33])

% DOF 2
figure('name','DOF 2');
a  = plot(omega_10/2/pi, A_vect_10(2, :),'ro-', 'LineWidth', 1, 'MarkerSize',1.5 );
grid on;
hold on;
bb = plot(omega/2/pi, A_vect(2, :), 'k-', 'LineWidth', 1.5);
b = plot(omega_30(1:569)/2/pi, A_vect_30(2, 1:569),'bo-', 'LineWidth', 1, 'MarkerSize', 1.5);
plot(omega_30(570:end)/2/pi, A_vect_30(2, 570:end),'bo-', 'LineWidth', 1, 'MarkerSize', 1.5);

d = plot(omega_50(1:583)/2/pi, A_vect_50(2, 1:583), 'o-', 'Color', "#7E2F8E", 'LineWidth', 1, 'MarkerSize', 1.5);
plot(omega_50(584:end)/2/pi, A_vect_50(2, 584:end),'o-', 'Color', "#7E2F8E", 'LineWidth', 1, 'MarkerSize',1.5);

xx = [15.91, 15.91];
yy = [0, max(A_vect_50(2, :))];
plot(xx, yy, 'k--', 'LineWidth', 1.2);

legend([a, b, d, bb], '$|F| = 10$ N', '$|F| = 30$ N', '$|F| = 50$ N', 'Backbone curve', 'Interpreter', 'latex', 'Fontsize', 12, 'Location', 'best')
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'Fontsize', 14)
ylabel('Amplitude [m]', 'Interpreter', 'latex', 'Fontsize', 14)
xlim([25 33])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basins of attractions
% Test to check that the time window is sufficient
sys.x_0 = 1e-2*[1, -1]';
sys.x_d_0 = [0, 0]';
freq = 29.5;
omega = 2*pi*freq;
sys.f_ext = @(t, omega) [30*cos(omega*t); 0];
y_0 = [sys.x_0; sys.x_d_0];
t_end = 8;
[tt, x] = time_integration(y_0, sys, dt, t_end, omega);

figure;
plot(tt, x);
legend('Dof 1', 'Dof 2')
%% Performing the basins of attractions
clc;
freq = 28.9;
omega = 2*pi*freq;
sys.f_ext = @(t, omega) [30*cos(omega*t); 0];       % Forcing definition
t_end = 8;
x_0_vect = -0.15 : 0.005 : 0.15;                    % ICs span

basins = basins_attractions(sys, dt, t_end, omega, x_0_vect);
%% Plots of basins
close all;
clc;
load("BASINS RES\30N_29.5HZ_Basins_0.15.mat")
figure;
heatmap(basins)
figure;
h = pcolor(x_0_vect, x_0_vect, basins);
set(h, 'EdgeColor', 'none');
colormap("jet");
cbar = colorbar();
%clim([0.005 16.5e-3]);
clim([0.000 22e-3]);
cbar.Ticks = [basins(1, 1) basins(1, end)]; %Create 8 ticks from zero to 1
cbar.TickLabels = ['Attractor 1' , 'Attractor 2'] ;    %Replace the labels of these 8 ticks with the numbers 1 
set(cbar,'XTickLabel',{'Attractor 1','Attractor 2'})
cbar.TickLabelInterpreter = 'latex';
cbar.FontSize = 10;
shading flat;
xlabel('$x_{0,2}$', 'Interpreter', 'latex', 'Fontsize', 16)
ylabel('$x_{0,1}$', 'Interpreter', 'latex', 'Fontsize', 16)

%%
close all;
load("CONTINUATION RES\v2_MAT_3frcs.mat");
figure('name','DOF 1');
grid on;
grid minor;
hold on;
a = plot(omega_30(1:569)/2/pi, A_vect_30(1, 1:569),'bo-', 'LineWidth', 1, 'MarkerSize', 1.5);
plot(omega_30(570:end)/2/pi, A_vect_30(1, 570:end),'bo-', 'LineWidth', 1, 'MarkerSize', 1.5);


xx = [15, 15];
yy = [0, max(A_vect_30(1, :))];
plot(xx, yy, 'k--', 'LineWidth', 1.2);
xx = [29.5, 29.5];
yy = [0, 0.03];
b = plot(xx, yy, 'k--', 'LineWidth', 1.2);

xx = [28.9, 28.9];
yy = [0, 0.03];
c = plot(xx, yy, 'k--', 'LineWidth', 1.2);

%legend([a, b], 'FRCs', ', 'Interpreter', 'latex', 'Fontsize', 12, 'Location', 'best')
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'Fontsize', 14)
ylabel('Amplitude [m]', 'Interpreter', 'latex', 'Fontsize', 14)
xlim([10 35])





