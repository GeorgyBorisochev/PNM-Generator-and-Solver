#pragma once
#include "framework.h"
#include "utils.h"

class classPNM
{
public:


	//process variables

	//subclass for network constants and settings
	class classSettings
	{
	public:
		//constants
		int seed = 100;
		double tol; //tolerance for lables definition as a fraction of the cubic domain size

		//input variables
		double domain; //nanometers
		double coordination;  //2x nodes/pores
		double anisotropy; //(0, 1]
		double Pinit; //initial domain pressure for initialsation of properties
		double Tinit; //domain temperature
		double _Pin;
		double _dP;

		int n_linear_iterations; //number of iterations for the PARDISO solver

		vector<vector<std::string>> poreSizeDist_raw;
		vector<vector<double>> poreSizeDist;
		double expectedP;
	};
	classSettings settings;

	vector<vector<std::string>> poreSizeDist_raw;
	vector<vector<double>> poreSizeDist;
	double expectedP;

	//subclass for network and geometry parameters
	class classNetwork
	{
	public:
		vector<bool> left, right, back, front, top, bottom;  //label faces (true or false)
		vector<vector<double>> pores;			//{ X, Y, Z, diameter}
		vector<vector<double>> throats;			//{beginning, end, diameter}
		vector<double> diameter_pore, diameter_throat; //pore and throat diameters
		vector<double> area_pore, area_troat; //cross-section area
		vector<double> length_throat;
		vector<vector<double>> conduit_length_throat; //conduit length for each throat {pore1, throat, pore2}
		vector<vector<double>> shape_factor_throat; //shape factor for each throat {pore1, throat, pore2}
		vector<vector<vector<double>>> endpoints_throat; //end coordinates for each throat {pore1: {x, y, z}, pore2: {x, y, z}}

		//functions
		/*classPNM length_throat_piecewise(classPNM PNM, classUtils utils);
		classPNM endpoints_spherical(classPNM PNM, classUtils utils);
		classPNM cross_sectional_area_simple(classPNM PNM, classUtils utils);
		classPNM conduit_length_throat_def(classPNM PNM, classUtils utils);
		classPNM shape_factor_cone_stick(classPNM PNM, classUtils utils);

		classPNM generate_network_properties(classPNM PNM, classUtils utils);*/
	};
	classNetwork network;

	//subclasses for Phase properties calculation
	class classMethane
	{
	public:

		//constants
		double MW = 0.016; //molecular weight [mol/kg]
		double Pc = 4599200;
		double Tc = 190.564;
		double Vc = 1 / 162.66;
		double d_mol = 0.38e-9; //m

		//network phase properties
		vector<double> pressure_pore, pressure_throat, temperature_pore, temperature_throat;
		vector<double> density_throat, density_pore, viscosity_throat, viscosity_pore, knudsen_throat, knudsen_pore;

		//phase functions
		//classPNM density_Peng_Robinson_Methane(classPNM PNM, classUtils utils);
		//classPNM viscosity_natural_gas_Lee(classPNM PNM, classUtils utils);
		//classPNM knudsen_number(classPNM PNM, classUtils utils);

		//classPNM generate_phase_properties(classPNM PNM, classUtils utils);
	};
	class classCarbonDioxide
	{
	public:

		//constants
		double MW = 0.016; //molecular weight [mol/hg]
		double Pc = 4599200;
		double Tc = 190.564;
		double Vc = 1 / 162.66;

		//network phase properties
		vector<double> pressure_pore, pressure_throat, temperature_pore, temperature_throat;
		vector<double> density_throat, density_pore, viscosity_throat, viscosity_pore, knudsen_throat, knudsen_pore;

		//phase functions
		/*classPNM density_Peng_Robinson_Methane(classPNM PNM, classUtils utils);
		classPNM viscosity_natural_gas_Lee(classPNM PNM, classUtils utils);
		classPNM knudsen_number(classPNM PNM, classUtils utils);*/
	};
	//phase selection
	classMethane phase;
	//classCarbonDioxide phase;

	class classPhysics
	{
	public:
		//physics properties
		vector<double> hydraulic_conductance_throat;

		//physics functions
		/*classPNM hydraulic_conductance_Ma_2014(classPNM PNM, classUtils utils);
		classPNM hydraulic_conductance_Song_2018(classPNM PNM, classUtils utils);

		classPNM generate_physics_properties(classPNM PNM, classUtils utils);*/

	};
	classPhysics physics;

	class classSinglePhaseStokesFlow
	{
	public:
		//variables
		vector<bool> inlets, outlets;
		double Pin, dP;

		vector<double> conductance;
		vector<double> flow_rate_pore, flow_rate_throat;
		vector<double> effective_permeability;

		//ALL in MTL4 datatypes
		//Matrix A
		compressed2D<double> A;

		//RHS vector b
		dense_vector<double> b;

		//current and previous solution vectors X, Xold
		dense_vector<double> X, Xold;
		dense_vector<double> pressure_BC;

		//functions
		/*classPNM create_A_b(classPNM PNM, classUtils utils);
		dense_vector<double> populate_A_b(compressed2D<double>& A, dense_vector<double> b, classPNM PNM);
		double laplacian_A(compressed2D<double>& A, classPNM PNM);
		classPNM solve_A_b(classPNM PNM, classUtils utils);
		classPNM pass_result_to_phase(classPNM PNM, classUtils utils);
		classPNM calculate_flow_rate(classPNM PNM, classUtils utils);
		classPNM calculate_effective_permeability(classPNM PNM, classUtils utils);
		classPNM run_once(classPNM PNM, classUtils utils);
		classPNM run_linear_iterative(classPNM PNM, classUtils utils);*/

	};
	classSinglePhaseStokesFlow algorithm;


	//functions

	
	double actualPorosity(classPNM PNM);
	classPNM propertyCorrection(classPNM PNM, classUtils utils);
	vector<int> poresPsize(classPNM PNM);
	classPNM label_faces(classPNM PNM, classUtils utils);
	classPNM initial_conditions(classPNM PNM, classUtils utils);

	classPNM test(classPNM PNM, classUtils utils);

	classPNM csvReadNetwork_python(string loadPath, classPNM PNM, classUtils utils);

	class classMethods
	{
	public:
		class classMethodsGeneration
		{
		public:
			void poresGeneration2(classPNM::classNetwork& network, classPNM::classSettings& settings);
			void throatsGeneration(classPNM::classNetwork& network, classPNM::classSettings& settings);
			void throatsGeneration2(classPNM::classNetwork& network, classPNM::classSettings& settings);
			void connectedComponents(classPNM::classNetwork& network);
		};
		classMethodsGeneration methodsGeneration;

		class classMethodsNetwork
		{
		public:
			//functions
			void length_throat_piecewise(classPNM::classNetwork& network);
			void endpoints_spherical(classPNM::classNetwork& network);
			void cross_sectional_area_simple(classPNM::classNetwork& network);
			void conduit_length_throat_def(classPNM::classNetwork& network);
			void shape_factor_cone_stick(classPNM::classNetwork& network);

			void generate_network_properties(classPNM::classNetwork& network, classPNM::classMethods::classMethodsNetwork& methodsNetwork);
		};
		classMethodsNetwork methodsNetwork;

		class classMethodsMethane
		{
		public:
			//phase functions
			void density_Peng_Robinson_Methane(classPNM::classMethane& phase, classPNM::classNetwork& network);
			void viscosity_natural_gas_Lee(classPNM::classMethane& phase, classPNM::classNetwork& network);
			void knudsen_number(classPNM::classMethane& phase, classPNM::classNetwork& network);

			void generate_phase_properties(classPNM::classNetwork& network, classPNM::classMethane& phase, classPNM::classMethods::classMethodsMethane& methodsPhase);
		};
		//choose phase here too!!
		classMethodsMethane methodsPhase;

		class classMethodsPhysics
		{
		public:
			void hydraulic_conductance_Ma_2014(classPNM::classNetwork& network, classPNM::classMethane& phase, classPNM::classPhysics& physics);
			void hydraulic_conductance_Song_2018(classPNM::classNetwork& network, classPNM::classMethane& phase, classPNM::classPhysics& physics);

			void generate_physics_properties(classPNM::classNetwork& network, classPNM::classMethane& phase, classPNM::classPhysics& physics, classPNM::classMethods::classMethodsPhysics& methodsPhysics);
		};
		classMethodsPhysics methodsPhysics;

		class classMethodsSinglePhaseStokesFlow
		{
		public:
			//functions
			void create_A_b(classPNM::classNetwork& network,
								classPNM::classPhysics& physics,
								classPNM::classSinglePhaseStokesFlow& algorithm,
								classPNM::classMethods::classMethodsSinglePhaseStokesFlow& methodsAlgorithm);
			dense_vector<double> populate_A_b(compressed2D<double>& A,
												dense_vector<double>& b,
												classPNM::classNetwork& network,
												classPNM::classSinglePhaseStokesFlow& algorithm,
												classPNM::classMethods::classMethodsSinglePhaseStokesFlow& methodsAlgorithm);
			double laplacian_A(compressed2D<double>& laplacian_A,
								classPNM::classNetwork& network,
								classPNM::classSinglePhaseStokesFlow& algorithm);
			void solve_A_b(classPNM::classNetwork& network,
								classPNM::classSinglePhaseStokesFlow& algorithm);
			void solve_A_b_pypardiso(classPNM::classNetwork& network,
										classPNM::classSinglePhaseStokesFlow& algorithm);
			void pass_result_to_phase(classPNM::classNetwork& network,
											classPNM::classMethane& phase,
											classPNM::classSinglePhaseStokesFlow& algorithm);
			void calculate_flow_rate(classPNM::classNetwork& network,
										classPNM::classMethane& phase,
										classPNM::classSinglePhaseStokesFlow& algorithm);
			void calculate_effective_permeability(classPNM::classNetwork& network,
													classPNM::classSettings& settings,
													classPNM::classMethane& phase,
													classPNM::classSinglePhaseStokesFlow& algorithm);
			void run_once(classPNM::classNetwork& network,
							classPNM::classSettings& settings,
							classPNM::classMethane& phase,
							classPNM::classPhysics& physics,
							classPNM::classSinglePhaseStokesFlow& algorithm,
							classPNM::classMethods::classMethodsSinglePhaseStokesFlow& methodsAlgorithm);
			void run_linear_iterative(classPNM::classNetwork& network,
										classPNM::classSettings& settings,
										classPNM::classMethane& phase,
										classPNM::classPhysics& physics,
										classPNM::classSinglePhaseStokesFlow& algorithm,
										classPNM::classMethods& methods);
		};
		//choose algorythm here too!!
		classMethodsSinglePhaseStokesFlow methodsAlgorithm;
	};
	classMethods methods;

	class classIO
	{
	public:
		void csvWriteNetwork(string savePath, classPNM::classNetwork& network, classPNM::classSettings& settings);
		void csvReadNetwork(string loadPath, classPNM::classNetwork& network, classPNM::classSettings& settings, classUtils utils);
		void readSettings(classPNM::classSettings& settings);
		void readSetupFile(classPNM::classSettings& settings, classUtils utils);
	};
	classIO IO;
};