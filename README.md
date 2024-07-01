# matlab_axon_research

The code of this repository consists of 3 main files:

squid_giant_axon.m
A Backword Euler discretization of the Hodgkin Huxley Partial Differential Equation describing the voltage of an axon in time and space.

dimensionless_squid_giant_axon.m 
A nondimensionalized Hodgkin Huxley Partial Differential Equation that has been discretized using Backward Euler method.

myelinated_HH_model.m
A Backword Euler discretization of the Hodgkin Huxley Partial Differential Equation where the 4 variables (capacitance, radius, sodium conductance and potassium conductance) are modified to functions of space. Here we may consider axons of a much smaller length and radius, and include the functions of space in an attempt to understand ion channel distribution and mimic myelinated axons.