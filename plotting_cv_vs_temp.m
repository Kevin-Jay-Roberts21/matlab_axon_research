% This script is used to plot the conduction velocities of the HH, SC, and
% DC models given temperature data.
% Kevin Roberts
% March 2025

clear all
close all
clc

% Original data for HH, SC, and DC models
HH_cv = [6.3 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32;
        13.6060718252499 13.9377934272300 14.4409937888199 14.9819086386251	15.5048076923077 16.0650281618024 16.5983606557377 17.1683226183518	17.7788220551378 18.3501683501684 18.9586357039187 19.4193061840121 20.0000000000000 20.6207482993197 21.0549645390071 21.6234967622572	22.1014492753623 22.6010101010101 22.8594080338266 23.2558139534884	23.5326688815061 23.8095238095238 23.9547038327526 23.8095238095238	23.6710963455150 22.9915433403805 21.6234967622572];
SC_cv = [20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50;
         0.404651922359089 0.439823096799841 0.476419413919414 0.513365360926337 0.553186650747627 0.588008158837956 0.626573617952929 0.664488017429194 0.697115384615384 0.727956989247312 0.764805097451275 0.793625914315569 0.805530676220332 0.810195551574862 0.787878787878788 0.690404040404040];
DC_cv = [20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50;
         0.402435610302352 0.437407637862643 0.476855992844365 0.513365360926336 0.552104141925179 0.589294727664847 0.624794745484401 0.664488017429194 0.700072843822844 0.727956989247312 0.755747126436782 0.793625914315570 0.805530676220331 0.805530676220332 0.778368147304647 0.694983619983621]; 

% Filter HH_cv data to only include values up to 25°C
HH_cv_filtered = HH_cv(:, HH_cv(1,:) <= 25);

% Create a new figure with a custom size
figure('Position', [100, 100, 1200, 400]); % [left, bottom, width, height]

% HH cv plot (filtered for temperatures up to 25°C)
subplot(1, 3, 1); % 1 row, 3 columns, first subplot
plot(HH_cv_filtered(1,:), HH_cv_filtered(2,:), 'b-o', 'LineWidth', 2); % Plot temperature vs HH cv
xlabel('Temperature (°C)');
ylabel('Conduction Velocity (m/s)');
title('HH cv vs Temp.');
grid on;

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