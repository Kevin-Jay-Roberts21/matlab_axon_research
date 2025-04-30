% Calculating R_i (resistivity)
% Kevin Roberts
% January 2025

clear all
close all
clc

% Constant parameters
R = 8.31454; % (J/(mole*K)) gas constant
T = 279.45; % (K) temperature in kelvin
F = 96485; % (C/mole) Faraday constant
D_Na = 1.33; % (um^2/ms) diffusion coefficient of sodium
D_K = 1.96; % (um^2/ms) diffusion coefficient of potassium
D_Cl = 2.0; % (um^2/ms) diffusion coefficient of chloride
z_Na = 1; % (dimless) valence of sodium
z_K = 1; % (dimless) valence of potassium
z_Cl = -1; % (dimless) valence of chloride

% Intracellular concentrations for squid and rat axons
squid_c_Na = 12; % (mM) intracellular concentration of sodium ions
squid_c_K = 155; % (mM) intracellular concentration of potassium ions
squid_c_Cl = 167; % (mM) intracellular concentration of chloride ions

rat_c_Na = 10; % (mM) intracellular concentration of sodium ions
rat_c_K = 140; % (mM) intracellular concentration of potassium ions
rat_c_Cl = 150; % (mM) intracellular concentration of chloride ions

% Compute R_i for squid and rat axons multiply by 10^8 to get kilo-ohm*cm
squid_R_i = 10^8*get_R_i(squid_c_Na, squid_c_K, squid_c_Cl, R, T, F, z_Na, D_Na, z_K, D_K, z_Cl, D_Cl)
rat_R_i = 10^8*get_R_i(rat_c_Na, rat_c_K, rat_c_Cl, R, T, F, z_Na, D_Na, z_K, D_K, z_Cl, D_Cl)

% Function to compute R_i
function R_i = get_R_i(c_Na, c_K, c_Cl, R, T, F, z_Na, D_Na, z_K, D_K, z_Cl, D_Cl)
    R_i = R * T / (F^2 * (z_Na^2 * D_Na * c_Na + z_K^2 * D_K * c_K + z_Cl^2 * D_Cl * c_Cl));
end