% computing C_m, R_m, C_my and R_my given a and a_my from Dr. Huang's using
% unit per length parameters given in Dr. Huang's paper.
% Kevin Roberts
% March 2025

clear all
close all
clc

a_T = 0.00015; % (cm) intracellular axon radius used for SiGe Tube
d_T = 0.00001; % (cm) thickness of the SiGe Tube
a_my_T = a_T + a_T; % (cm) effective axon radius used for SiGe Tube

a_TP = 0.0003; % (cm) intracellular axon radius used for Tube+Paralyne
d_TP = 0.00005; % (cm) thickness of the Tube+Paralyne
a_my_TP = a_TP + a_TP; % (cm) effective axon radius used for Tube+Paralyne

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

function R_my compute_R_my(a_my, )

end 





function compute_all_Cm_values(a)






