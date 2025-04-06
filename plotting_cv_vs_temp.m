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
SC_cv = [20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 55;
         0.336788883981201	0.362930489214596	0.391565446591125	0.418661872150244	0.448992035838726	0.478236525634159	0.509018759018759	0.536669213139802	0.567540322580645	0.598864535343664	0.627413127413127	0.648148148148148	0.674531351001939	0.689915966386555	0.710784313725490	0.702450980392157	0.694047619047619	0.640650772229720	0.564350987370202];
DC_cv = [20 21 22 23;
         0.133183726139911	0.131413889193392	0.126295159866401	0.112502757586635]; 


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