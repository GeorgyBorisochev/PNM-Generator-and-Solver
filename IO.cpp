#include "classPNM.h"

//
bool isNumber(const std::string& s)
{
	char* end = nullptr;
	double val = strtod(s.c_str(), &end);
	return end != s.c_str() && *end == '\0' && val != HUGE_VAL;
}

//FUNCTIONS FOR SAVING AND READING NETWORKS AND GENERATOR INPUTS


void classPNM::classIO::csvWriteNetwork(string savePath, classPNM::classNetwork& network, classPNM::classSettings& settings)
{
	ofstream results;
	results.open(savePath + ".csv");

	//Header
	results << "Line entries: 1 - xPore; 2 - yPore; 3 - zPore; 4 - diameterPore; 5 - startThroat; 6 - endThroat; 7 - diameterThroat; 8 - other settings\n";

	//pores X
	results << network.pores[0][0] ;
	for (int i = 1; i < network.pores.size(); i++)
	{
		results << "," << network.pores[i][0] ;
	}
	results << "\n";
	//pores Y
	results << network.pores[0][1];
	for (int i = 1; i < network.pores.size(); i++)
	{
		results << "," << network.pores[i][1];
	}
	results << "\n";
	//pores Z
	results << network.pores[0][2];
	for (int i = 1; i < network.pores.size(); i++)
	{
		results << "," << network.pores[i][2];
	}
	results << "\n";
	//pores Diameter
	results << network.pores[0][3];
	for (int i = 1; i < network.pores.size(); i++)
	{
		results << "," << network.pores[i][3];
	}
	results << "\n";
	
	//throats Start
	results << (int)network.throats[0][0];
	for (int i = 1; i < network.throats.size(); i++)
	{
		results << "," << (int)network.throats[i][0] ;
	}
	results << "\n";
	//throats End
	results << (int)network.throats[0][1];
	for (int i = 1; i < network.throats.size(); i++)
	{
		results << "," << (int)network.throats[i][1];
	}
	results << "\n";
	//throats Diameter
	results << network.throats[0][2];
	for (int i = 1; i < network.throats.size(); i++)
	{
		results << "," << network.throats[i][2];
	}
	results << "\n";
	//domain size
	results << settings.domain;

	results.close();
}

void classPNM::classIO::csvReadNetwork(string loadPath, classPNM::classNetwork& network, classPNM::classSettings& settings, classUtils utils)
{
	network.pores.clear();
	network.throats.clear();
	network.diameter_pore.clear();
	network.diameter_throat.clear();

	loadPath = loadPath + ".csv";

	std::vector<std::vector<std::string> > raw = utils.csvRead(loadPath);

	//converting string into double PSD
	vector<double> temp;

	//pore data

	for (int i = 0; i < raw[1].size(); i++)
	{
		temp.clear();
		temp.push_back(stod(raw[1][i])); //x
		temp.push_back(stod(raw[2][i])); //y
		temp.push_back(stod(raw[3][i])); //z
		temp.push_back(stod(raw[4][i])); //D
		network.pores.push_back(temp);
		network.diameter_pore.push_back(stod(raw[4][i]));
	}
	temp.clear();
	for (int i = 0; i < raw[5].size(); i++)
	{
		temp.clear();
		temp.push_back(stoi(raw[5][i])); //start
		temp.push_back(stoi(raw[6][i])); //end
		temp.push_back(stod(raw[7][i])); //D
		network.throats.push_back(temp);
		network.diameter_throat.push_back(stod(raw[7][i]));
	}
	settings.domain = stod(raw[8][0]);
}

void classPNM::classIO::readSettings(classPNM::classSettings& settings)
{
	string in;
	double temp;

	std::cout << "\nProvide Network domain size in [nm] between 10 and 1000:\n";
	getline(std::cin, in);
	while ((!isNumber(in) || (stod(in) < 10) || stod(in) > 1000) && !in.empty())
	{
		if (!isNumber(in))
		{
			std::cout << "Not a number, try again:\n";
			getline(std::cin, in);
		}
		else
		{
			std::cout << "Outside of range, try again:\n";
			getline(std::cin, in);
		}
	}
	if (!in.empty())
	{
		settings.domain = stod(in)*1e-9;
	}
	else
	{
		std::cout << "Using default domain size of 100nm!";
		settings.domain = 100e-9;
	}

	std::cout << "\nProvide Network coordination number between 2 and 6:\n";
	getline(std::cin, in);
	while ((!isNumber(in) || (stod(in) < 2) || (stod(in) > 6)) && !in.empty())
	{
		if (!isNumber(in))
		{
			std::cout << "Not a number, try again:\n";
			getline(std::cin, in);
		}
		else
		{
			std::cout << "Outside of range, try again:\n";
			getline(std::cin, in);
		}
	}
	if (!in.empty())
	{
		settings.coordination = stod(in);
	}
	else
	{
		std::cout << "Using default domain size of 2.51!\n";
		settings.coordination = 2.51;
	}

}

void classPNM::classIO::readSetupFile(classPNM::classSettings& settings, classUtils utils)
{
	//string loadPath = "C:/Users/gb168/Desktop/PhD research/Network Generation application/PNMGen/x64/Release/SETUP.csv";
	string loadPath = "SETUP.csv";

	std::vector<std::vector<std::string> > raw = utils.csvRead(loadPath);
	vector<string> temp_vec;
	string temp;
	for (int i = 0; i < raw.size(); i++)
	{

		//generator settings 
		if (raw[i][0] == "PSD_diameter")
		{
			settings.poreSizeDist_raw.clear();
			temp_vec = raw[i]; //diameter
			temp_vec.erase(temp_vec.begin(), temp_vec.begin() + 2);
			settings.poreSizeDist_raw.push_back(temp_vec);
			temp_vec = raw[i+1]; //specific porosity
			temp_vec.erase(temp_vec.begin(), temp_vec.begin() + 2);
			settings.poreSizeDist_raw.push_back(temp_vec);
			settings.poreSizeDist_raw = utils.transpose(settings.poreSizeDist_raw);
		}
		else if(raw[i][0] == "Domain_size")
		{
			temp = raw[i][2];
			settings.domain = stod(temp)*1e-9;
		}
		else if (raw[i][0] == "Coordination")
		{
			temp = raw[i][2];
			settings.coordination = stod(temp);
		}
		else if (raw[i][0] == "Anisotropy")
		{
			temp = raw[i][2];
			settings.anisotropy = stod(temp);
		}
		

		//solver settings
		else if (raw[i][0] == "Pressure_in")
		{
			temp = raw[i][2];
			settings._Pin = stod(temp);
		}
		else if (raw[i][0] == "Pressure_drop")
		{
			temp = raw[i][2];
			settings._dP = stod(temp);
		}
		else if (raw[i][0] == "Temperature")
		{
			temp = raw[i][2];
			settings.Tinit = stod(temp);
		}

		//advanced settings
		else if (raw[i][0] == "Boundary")
		{
			temp = raw[i][2];
			settings.tol = stod(temp);
		}
		else if (raw[i][0] == "Seed")
		{
			temp = raw[i][2];
			settings.seed = stoi(temp);
		}
		else if (raw[i][0] == "Linear_iterations")
		{
			temp = raw[i][2];
			settings.n_linear_iterations = stoi(temp);
		}
	}
}

