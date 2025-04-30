% computing C_m, R_m, C_my and R_my given a and a_my from Dr. Huang's using
% unit per length parameters given in Dr. Huang's paper.
% Kevin Roberts
% March 2025

clear all
close all
clc

a = 0.0001; % (cm) intracellular axon radius for myelin
d_m = 0.0000238; % (cm) thickness of the myelin
a_my = a + d_m; % (cm) effective axon radius when myelin is included

a_T = 0.0001; % (cm) intracellular axon radius used for SiGe Tube
d_T = 0.00001; % (cm) thickness of the SiGe Tube
a_my_T = a_T + d_T; % (cm) effective axon radius used for SiGe Tube

a_TP = 0.0001; % (cm) intracellular axon radius used for Tube+Paralyne
d_TP = 0.00003; % (cm) thickness of the Tube+Paralyne
a_my_TP = a_TP + d_T + d_TP; % (cm) effective axon radius used for Tube+Paralyne

% capacitance and resistance per unit length that was given for
% unmyelianted, myelinated, SiGe Tube and Tube+Paralyne
% Note: *_unmyelianted1 and *_unmyelianted2 are the upper and lower
% bounds given in the literature for c_m and r_m unmyelinated.
c_my_myelinated = 0.00000628; % (uF/cm)
c_my_T = 0.000123; % (uF/cm)
c_my_TP = 0.0000153; % (uF/cm)

r_my_myelinated = 1.59*10^5; % (kilo-ohms * cm)
r_my_T = 3.92*10^10; % (kilo-ohms * cm)
r_my_TP = 6*10^13; % (kilo-ohms * cm)

% the computed lists give values in a vector in the following order, where
% LB and UB are lower bound and upper bound:
% Unmyelinated LB, Unmyelinated UB, Myelinated, SiGe Tube, Tube+Paralyne

% computing myelin capacitances and resistances
M_computed_resistance = compute_resistance(a_my, r_my_myelinated)
M_computed_capacitance = compute_capacitance(a_my, c_my_myelinated)

% computing SiGe Tube capacitances and resistances
T_computed_resistance = compute_resistance(a_my_T, r_my_T)
T_computed_capacitance = compute_capacitance(a_my_T, c_my_T)

% computing Tube+Paralyne capacitances and resistances
TP_computed_resistance = compute_resistance(a_my_TP, r_my_TP)
TP_computed_capacitance = compute_capacitance(a_my_TP, c_my_TP)


function resistance = compute_resistance(a_my, r_my)
    resistance = 2*pi*a_my*r_my;
end 

function capacitance = compute_capacitance(a_my, c_my)
        capacitance = c_my/(2*pi*a_my);
end 






