PNMSolve 0.5
Georgy Borisochev, 2022
Geoenergy Research Group, Heriot-Watt University

Integrated software for generating running pore network models based on bulk rock pore characterisation data.
Main features include:
- generation of completely randomised PNM based on a variety of inputs provided by experimental data
- inclusion of real gas effects(Song, 2018) and gas adsorption

Instructions for use:
- Software allows to change input parameters and re-run the model without closing the program. Simply adjust and save the SETUP.csv and press Enter when prompted
- All input is provided in SETUP.csv, make sure it stays in the same folder as the .exe file
- easiest is to open .csv in Excel and ignore "Do you want to keep using that format?" prompt to make sure there are no delimiter issues when reading inputs
- there is little Error feedback at the moment. If software crashes unexpectedly, it is advised to revert SETUP.csv to its original state as distributed

Current input variables:

GENERATOR SETTINGS	UNIT	EXAMPLE		NOTES //settings used for generating the model, they both affect the geometry and resulting flow properties
PSD_diameter		nm	62.97865132	//pore diameter values from pore size distribution, in nanometers
PSD_specific_porosity	%	0.006948682	//respective specific porosity values. Currently only porosity is supported, bu tin the future selecting Incremental Pre Volume will be available too
Domain_size		nm	100		//side length of the model cubic domain, in nanometres. Software performs as O(n^3) based on size. Currently algorithm is developed to make generation of 1000nm models accessible,
						//standing at around 15 minute generation time. Generating larger models is not advised at the moment. 
Coordination			2.51		//Z = 2*n_throats/n_pores
Anisotropy			1		//anisotropy parameter between 0 (cannot be equal to 0!) and 1. Lower value means preferred troat connection in XY-plane, simulating mudrock bedding
SOLVER SETTINGS		UNIT			NOTES //settings for running the model, only affect the flow proprties
Pressure_in		Pa	1.00E+06	//upstream pressure, in Pascals
Pressure_drop		Pa	10000		//pressure drop across the model 
Temperature		K	400		//absolute model temperature, in Kelvin
ADVANSED SETUP		UNIT			NOTES //additional setup
Boundary		frac	0.05		//parameter for labelling boundary pores when network is generated. To be changed with care
Seed				100		//seed for network generation. Input any positive integer to create different model. Using same seed will result in same results, useful for tracjing effects of a single variable
Linear_iterations		2		//positive integer, number of linear iterations to be run. Properties are regenerated between iterations, higher value leads to better convergence 

Compilation notes and required libs:
- Using x64
- Boost, https://www.boost.org/users/history/version_1_78_0.html
- MTL4, installed as a part of boost package, https://github.com/simunova/mtl4
- Intel OneAPI MKL, additional dependencies: (mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib), https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#base-kit