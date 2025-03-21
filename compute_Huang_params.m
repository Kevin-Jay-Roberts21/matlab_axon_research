% computing C_m, R_m, C_my and R_my given a and a_my from Dr. Huang's using
% unit per length parameters given in Dr. Huang's paper.
% Kevin Roberts
% March 2025

clear all
close all
clc

a_T = 0.00015; % (cm) intracellular axon radius used for SiGe Tube
d_T = 0.00001; % (cm) thickness of the SiGe Tube
a_my_T = a_T + d_T; % (cm) effective axon radius used for SiGe Tube

a_TP = 0.0003; % (cm) intracellular axon radius used for Tube+Paralyne
d_TP = 0.00005; % (cm) thickness of the Tube+Paralyne
a_my_TP = a_TP + d_TP; % (cm) effective axon radius used for Tube+Paralyne

% capacitance and resistance per unit length that was given for
% unmyelianted, myelinated, SiGe Tube and Tube+Paralyne
% Note: *_unmyelianted1 and *_unmyelianted2 are the upper and lower
% bounds given in the literature for c_m and r_m unmyelinated.
c_m_unmyelinated1 = 0.00063; % (uF/cm)
c_m_unmyelinated2 = 0.00125; % (uF/cm)
c_my_myelinated = 0.00000628; % (uF/cm)
c_my_T = 0.000123; % (uF/cm)
c_my_TP = 0.0000153; % (uF/cm)

r_m_unmyelinated1 = 200; % (kilo-ohms * cm)
r_m_unmyelinated2 = 260; % (kilo-ohms * cm)
r_my_myelinated = 159000; % (kilo-ohms * cm)
r_my_T = 3.92*10^10; % (kilo-ohms * cm)
r_my_TP = 6*10^13; % (kilo-ohms * cm)

list_of_cms = [c_m_unmyelinated1, c_m_unmyelinated2, c_my_myelinated, c_my_T, c_my_TP];
list_of_rms = [r_m_unmyelinated1, r_m_unmyelinated2, r_my_myelinated, r_my_T, r_my_TP];

% the computed lists give values in a vector in the following order, where
% LB and UB are lower bound and upper bound:
% Unmyelinated LB, Unmyelinated UB, Myelinated, SiGe Tube, Tube+Paralyne

% computing SiGe Tube capacitances and resistances
T_computed_resistances = compute_resistances(a_T, a_my_T, list_of_rms)
T_computed_capacitances = compute_capacitances(a_T, a_my_T, list_of_cms)

% computing Tube+Paralyne capacitances and resistances
TP_computed_resistances = compute_resistances(a_TP, a_my_TP, list_of_rms)
TP_computed_capacitances = compute_capacitances(a_TP, a_my_TP, list_of_cms)


function resistances = compute_resistances(a, a_my, list_of_rms)
    
    resistances = zeros(1, length(list_of_rms));
    
    % computing R_m
    resistances(1) = 2*pi*a*list_of_rms(1);
    resistances(2) = 2*pi*a*list_of_rms(2);
    
    % computing R_my
    for i = 3:length(list_of_rms)
        resistances(i) = 2*pi*a_my*list_of_rms(i);
    end

end 

function capacitances = compute_capacitances(a, a_my, list_of_cms)
    
    capacitances = zeros(1, length(list_of_cms));
    
    % computing C_m
    capacitances(1) = list_of_cms(1)/(2*pi*a);
    capacitances(2) = list_of_cms(2)/(2*pi*a);
    
    % computing C_my
    for i = 3:length(list_of_cms)
        capacitances(i) = list_of_cms(i)/(2*pi*a_my);
    end

end 






