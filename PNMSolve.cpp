// PNMSolve.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "framework.h"
#include "utils.h"
#include "CSVRead.h"
#include "classPNM.h"

std::string currentDateTime() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

	return buf;
}

int main()
{

	classUtils utils;
	chrono::steady_clock::time_point begin;
	chrono::steady_clock::time_point end;
	string elapsedTime;

	cout << "PNMSolve 0.5\nSoftware for generation and simulation of pore network models.\n";
	cout << "Copyright Georgy Borisochev, 2022\nGeoenergy Research Group, Heriot-Watt University\n\n";
	cout << "This software version allows for continious running of models,\nUse SETUP.csv to change model parameters.\n\n";
	cout << "Press Enter to start the program:";
	getchar();

	int exit_flag = 0;
	while (exit_flag == 0)
	{
		
		classPNM PNM;

		PNM.IO.readSetupFile(PNM.settings, utils);



		//run generation algorithm

		bool generation = true;  //flag to generate or load model

		if (generation)
		{
			utils.writeLine("Generating Network...");
			PNM.methods.methodsGeneration.poresGeneration2(PNM.network, PNM.settings);
			PNM.methods.methodsGeneration.throatsGeneration2(PNM.network, PNM.settings);
			utils.writeLine("Network Generated");

			//remove unconnected components
			PNM.methods.methodsGeneration.connectedComponents(PNM.network);

			//adjust to a specific total porosity 
			PNM = PNM.propertyCorrection(PNM, utils);

			string name = currentDateTime().substr(0, 4) + currentDateTime().substr(5, 2) + currentDateTime().substr(8, 2) + "_" + to_string((int)(PNM.settings.domain * 1e9)) + "_nm_stochastic";
			//string savePath = "C:/Users/gb168/Desktop/Opalinus_OpenPNM/Opalinus_521_multi/cpp_networks/";
			string savePath = "";
			savePath += name;
			PNM.IO.csvWriteNetwork(savePath, PNM.network, PNM.settings);
			utils.writeLine("Network Saved to " + savePath);
			utils.emptyLine();
		}
		else
		{
			string name = "20220124_100_nm_stochastic";
			string Path = "C:/Users/gb168/Desktop/Opalinus_OpenPNM/Opalinus_521_multi/cpp_networks/";
			string loadPath = Path + name;
			//string loadPath = "C:/Users/gb168/Desktop/Opalinus_OpenPNM/Opalinus_521_multi/test_small";
			PNM.IO.csvReadNetwork(loadPath, PNM.network, PNM.settings, utils);
			utils.writeLine("Network Uploaded");
			utils.emptyLine();
		}
		
		cout << setprecision(4) << std::fixed;
		//finalise the network properties
		PNM.methods.methodsNetwork.generate_network_properties(PNM.network, PNM.methods.methodsNetwork);

		//label boundary faces
		PNM = PNM.label_faces(PNM, utils);

		//calculation of initial phase properties 
		PNM.algorithm.Pin = PNM.settings._Pin;
		PNM.algorithm.dP = PNM.settings._dP;
		PNM.settings.Pinit = PNM.algorithm.Pin - PNM.algorithm.dP / 2;
		PNM = PNM.initial_conditions(PNM, utils);

		PNM.methods.methodsPhase.generate_phase_properties(PNM.network, PNM.phase, PNM.methods.methodsPhase);

		//calculate initial physics properties
		PNM.methods.methodsPhysics.generate_physics_properties(PNM.network, PNM.phase, PNM.physics, PNM.methods.methodsPhysics);

		//TEST NETWORK
		//PNM = PNM.test(PNM, utils);

		//initalise algorithm
		PNM.algorithm.inlets = PNM.network.left;
		PNM.algorithm.outlets = PNM.network.right;
		/*PNM.algorithm.inlets = PNM.network.bottom;
		PNM.algorithm.outlets = PNM.network.top;*/

		begin = chrono::steady_clock::now();

		PNM.methods.methodsAlgorithm.run_linear_iterative(PNM.network, PNM.settings, PNM.phase, PNM.physics, PNM.algorithm, PNM.methods);
		//PNM.methods.methodsAlgorithm.run_once(PNM.network, PNM.settings, PNM.phase, PNM.physics, PNM.algorithm, PNM.methods.methodsAlgorithm);

		end = chrono::steady_clock::now();
		elapsedTime = "Elapsed time: " + to_string(chrono::duration_cast<std::chrono::seconds>(end - begin).count()) + " [s]";
		cout << elapsedTime << endl;

		cout << "Finish!" << endl;
		cout << "Press Enter to repeat, any other key to exit:\n";
		string in;
		getline(std::cin, in);
		if (!in.empty())
		{
			exit_flag++;
		}
	}
}