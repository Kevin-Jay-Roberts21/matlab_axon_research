# matlab_axon_research

This code was created to numerically solve the Hodgkin-Huxley (HH) (used for nonmyelianted axons), Single Cable (SC) and Double Cable (DC) (used for myelinated axons) models. Below we provide a detailed desciption of the contents of each folder, and how to use each of the codes for personal experimentation.

## main.m

In this file, one can may run an HH, SC or DC simulation. One must pick a mesh parameter set, a material parameter set, and a numerical scheme (some base mesh and material parameter sets are provided in \mesh_parameter_sets and \material_parameter_sets, as well as numerical schemes in \numerical_scheme_functions). After the sets and scheme are chosen, one can run a simulation and the resulting data will automatically be save (Vm, n, m, and h matrices for hh simulation, and Vm, Vmy, n, m and h for SC and DC simulation). The saved data can be later implemented in the \plotters to observe simulation results, or used in the \computation_scripts to compute cv and other quantitatives.

## \numerical_scheme_functions

Contains different schemes used on the HH, SC and DC models. We note that v1 for HH is implicit, but every version of the SC and DC models are semi-implicit since the models are coupled PDEs, where v1 is the "least implicit" and v3 is the "most implicit". More details to come on what "most" and "least" mean here.

## \mesh_parameter_sets

Different meshes can be experimented with here where we set hh_mesh_parameters and sc_dc_mesh_parameters as the "base" mesh parameters.

## \material_parameter_sets

A collection of material parameter sets used on the HH, SC and DC models. Of course, more parameters are included in the DC and SC than in the HH parameter set.

## \all_schemes

In each of the files here, we include the mesh, material parameters, and the schemes all in one file. These files are nice to use when changing one specific parameter. We include extra schemes here that are still under further experimentation and investigation such as stochastic_hh_code.m and myelinated_hh_code.m. Additionally, we include a nondimensionalized hh code dimensionless_hh_code.m.

## \computation_scripts

This folder includes files used to compute other interesting axon features like Vm equilibrium (computing_equilibrium.m) based on nernst potentials and conductances, the time it takes for an action potential to be generated (computing_stimulus_threshold_times.m), converting per unit length to material parameters (computing_unit_to_material.m), calculating the axon intracellular resistivity based on material parameters (computing_resistivity.m), and finally computing the conduction velocity given some Vm data (computing_conduction_velocity.m). Computing the conduction velocity and the time it takes for an action potential to generate require data from a simultion, whereas the rest do not.

## \plotters

Consists of of many functions used to plot time shots, space shots, and even animation given data from a simulation.



