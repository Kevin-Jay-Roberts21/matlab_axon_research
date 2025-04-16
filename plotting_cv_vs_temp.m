% This script is used to plot the conduction velocities of the HH, SC, and
% DC models given temperature data.
% Kevin Roberts
% March 2025

clear all
close all
clc

% Original data for HH, SC, and DC models
HH_cv1 = [6.3 8 10 12 14 16 18 20 22 24 26 28 30 31 32;
          13.6060718252499 14.4927536231884	15.5048076923077 16.6666666666667 17.7788220551378 18.9586357039187 20 21.1657801418440 22.1014492753623 22.9915433403805 23.6710963455149 23.8095238095238 23.5326688815061 22.9915433403805 21.6234967622572]
HH_cv = [6.3 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32;
        13.6060718252499 13.9377934272300 14.4409937888199 14.9819086386251	15.5048076923077 16.0650281618024 16.5983606557377 17.1683226183518	17.7788220551378 18.3501683501684 18.9586357039187 19.4193061840121 20.0000000000000 20.6207482993197 21.0549645390071 21.6234967622572	22.1014492753623 22.6010101010101 22.8594080338266 23.2558139534884	23.5326688815061 23.8095238095238 23.9547038327526 23.8095238095238	23.6710963455150 22.9915433403805 21.6234967622572];
SC_cv = [20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 51 53 54 55 56 57;
         0.631799234400533	0.679962721342032	0.725762527233116	0.781974637681160	0.818840579710145	0.882034632034632	0.918498168498169	0.964411027568922	1.01832399626517	1.06550802139037	1.11309523809524	1.14583333333333	1.19444444444444	1.23141186299081	1.25522138680033	1.24060150375940	1.21679197994988	1.16666666666667	1.13203463203463	1.03670634920635	0.893087898775071];
DC_cv = [20 21 22 23;
         0.1332    0.1314    0.1265    0.1126]; 


% Create a new figure with a custom size
figure('Position', [100, 100, 1200, 400]); % [left, bottom, width, height]

% HH cv plot (filtered for temperatures up to 25°C)
subplot(1, 3, 1); % 1 row, 3 columns, first subplot
plot(HH_cv1(1,:), HH_cv1(2,:), 'b-o', 'LineWidth', 2); % Plot temperature vs HH cv
xlabel('Temperature (°C)');
ylabel('Conduction Velocity (m/s)');
title('HH cv vs Temp.');
grid on;

% % HH cv plot (filtered for temperatures up to 25°C)
% subplot(1, 3, 1); % 1 row, 3 columns, first subplot
% plot(HH_cv(1,:), HH_cv(2,:), 'b-o', 'LineWidth', 2); % Plot temperature vs HH cv
% xlabel('Temperature (°C)');
% ylabel('Conduction Velocity (m/s)');
% title('HH cv vs Temp.');
% grid on;

% SC cv plot
subplot(1, 3, 2); % 1 row, 3 columns, second subplot
plot(SC_cv(1,:), SC_cv(2,:), 'r-o', 'LineWidth', 2); % Plot temperature vs SC cv
xlabel('Temperature (°C)');
ylabel('Conduction Velocity (m/s)');
title('SC cv vs Temp.');
grid on;

% DC cv plot
subplot(1, 3, 3); % 1 row, 3 columns, third subplot
plot(DC_cv(1,:), DC_cv(2,:), 'g-o', 'LineWidth', 2); % Plot temperature vs DC cv
xlabel('Temperature (°C)');
ylabel('Conduction Velocity (m/s)');
title('DC cv vs Temp.');
grid on;

% Adjust layout for better visibility
tight_layout();