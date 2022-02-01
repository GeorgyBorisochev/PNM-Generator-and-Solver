#include "classPNM.h"

// Forward declarations of functions included in this code module:
bool sortColumn0(const vector<double>& v1, const vector<double>& v2) { return v1[0] < v2[0]; }
bool sortColumn2(const vector<double>& v1, const vector<double>& v2) { return v1[2] < v2[2]; }

vector<int> poresPsize_(vector<vector<double>>& poreSizeDist, double coordination, double domain)
{
	double pi = 3.1415926;
	vector<int> result;
	double volume_per_element;
	double throatLengthFactor = 2.1; //ratio of throat length to pore diameter 

	for (int i = 0; i < poreSizeDist.size(); i++)
	{
		volume_per_element = 4 / 3 * pi*pow(poreSizeDist[i][0] / 2, 3) + throatLengthFactor * coordination / 2 * pi*pow(poreSizeDist[i][0] / 2, 2)*(poreSizeDist[i][0]);
		result.push_back(round(poreSizeDist[i][1] * pow(domain, 3) / volume_per_element));

	}
	return result;

}

//function for determining generated pore diameter, given the current PSD entry, upper and lower bounds
double diameter_pore_gen(int current_psd_index, classPNM::classSettings& settings)
{
	//PSD is provided from biggst to smallest pore size
	double d_mean, d_max, d_min, diameter_pore;
	d_mean = settings.poreSizeDist[current_psd_index][0];
	int upper_limit = settings.poreSizeDist.size() - 1; //for bounds def

	//define maximum bound
	if ((current_psd_index > 0) && (current_psd_index < upper_limit))
	{
		d_max = settings.poreSizeDist[current_psd_index - 1][0];
		d_min = settings.poreSizeDist[current_psd_index + 1][0];
	}
	else if (current_psd_index <= 0)
	{
		d_min = settings.poreSizeDist[current_psd_index + 1][0];
		d_max = 2 * d_mean - d_min;
	}
	else if (current_psd_index >= upper_limit)
	{
		d_max = settings.poreSizeDist[current_psd_index - 1][0];
		d_min = 2 * d_mean - d_max;
	}
	d_max = 0.5*(d_max + d_mean); //halfway point between PSD entries
	d_min = 0.5*(d_min + d_mean);

	double ratio = (d_mean - d_min) / (d_max - d_min); //if deviations to left and right are uneven, it is important to generate on "smaller" deviation side to make sure over large enough sample, average result is d_mean
	
	double sample = ((double)rand() / RAND_MAX);//random saple in [0, 1] to define on which side of d_mean we are generating

	//ratio is calculated from d_min, but selection of the side is done in inverse (froe example, if ratio = 0.1, that means 0.9 times out of 1 we need to generate on the d_min side to even out the result)
	if (sample < ratio) //generating on d_max side
	{
		double x = ((double)rand() / RAND_MAX);
		diameter_pore = d_mean + (d_max - d_mean) * x;
	}
	else if (sample > ratio) //generating on d_min side
	{
		double x = ((double)rand() / RAND_MAX);
		diameter_pore = d_min + (d_mean - d_min) * x;
	}
	else { diameter_pore = d_mean; }

	//diameter_pore = d_mean;
	return diameter_pore;
}

///
///ClassGeneration functions
///

void classPNM::classMethods::classMethodsGeneration::poresGeneration2(classPNM::classNetwork& network, classPNM::classSettings& settings)
{
	chrono::steady_clock::time_point begin = chrono::steady_clock::now();

	settings.poreSizeDist.clear();
	network.pores.clear();
	network.diameter_pore.clear();

	//converting string into double PSD
	vector<double> temp;
	for (int i = 0; i < settings.poreSizeDist_raw.size(); i++)
	{
		temp.clear();
		for (int j = 0; j < settings.poreSizeDist_raw[i].size(); j++)
		{
			if (j == 0) { temp.push_back(stod(settings.poreSizeDist_raw[i][j])*1e-9); }
			else { temp.push_back(stod(settings.poreSizeDist_raw[i][j])); }
		}
		settings.poreSizeDist.push_back(temp);
	}

	//calculate number of pores per each poe size with given domain
	vector<int> nPores = poresPsize_(settings.poreSizeDist, settings.coordination, settings.domain);

	//variable definition
	int total_pn, current_psd_index, nCells, xCell, yCell, zCell, cellIndex, nCellsSmall, xCellSmall, yCellSmall, zCellSmall, cellIndexSmall;
	double expected_porosity, cellSize, pores_per_cell, cellSizeSmall;
	int nCellsMicro, xCellMicro, yCellMicro, zCellMicro, cellIndexMicro;
	double cellSizeMicro;

	total_pn = 0;
	expected_porosity = 0;
	//calculate total pore number
	for (int i = 0; i < nPores.size(); i++)
	{
		if (!(nPores[i] < 1))
		{
			total_pn += nPores[i];
			expected_porosity += settings.poreSizeDist[i][1];
		}
	}
	settings.expectedP = expected_porosity;
	cout << "Expected network porosity is: " << expected_porosity * 100.0 << " %\n";

	//define which pore size we are starting generation from
	current_psd_index = 0;
	while (nPores[current_psd_index] < 1)
	{
		current_psd_index++;
	}

	//define how many pores we want to generate as a "large" grid:
	double fraction = 0.0;
	double cut_off = current_psd_index; //everything up to but NOT including THIS entry into PSD to be generated as "large" grid
	double cut_off_micro = current_psd_index; ////everything up to but NOT including THIS entry into PSD to be generated as "small" grid
	int large_pores = 0; //number of large pores to be generated
	int small_pores = 0;
	while (fraction < 0.001)
	{
		fraction += (double)nPores[cut_off] / total_pn;
		large_pores += nPores[cut_off];
		cut_off++;

	};
	//cout << "Cut-off: " << cut_off << endl;

	//divide donain into "large" cells, just larger than biggest pore radius
	cellSize = settings.poreSizeDist[current_psd_index][0] / 2;
	nCells = (int)floor(settings.domain / cellSize);
	cellSize = settings.domain / nCells;

	//divide domain into "small" cells, based on the cut-off pore size
	cellSizeSmall = settings.poreSizeDist[cut_off][0] / 2;
	nCellsSmall = (int)floor(settings.domain / cellSizeSmall);
	cellSizeSmall = settings.domain / nCellsSmall;

	//variable for storing each "large" cell pore coords+diameter
	vector<vector<vector<double>>> cellCoords(pow(nCells, 3));

	//variable for storing each "small" cell pore coords+diameter
	vector<vector<vector<double>>> cellCoordsSmall(pow(nCellsSmall, 3));


	int n_generated = 0;

	//give seed to the random generation algorithm
	srand(settings.seed);

	//generate the first pore
	//variable for storing randomly generated pore coords plus its diameter before pushing into the Pores array
	vector<double> temp_pore = { ((double)rand() / RAND_MAX) * settings.domain, ((double)rand() / RAND_MAX) * settings.domain, ((double)rand() / RAND_MAX) * settings.domain };
	double temp_diameter = diameter_pore_gen(current_psd_index, settings); //generate pore size
	temp_pore.push_back(temp_diameter); //assign value produced by the genertation function
	network.pores.push_back(temp_pore);
	network.diameter_pore.push_back(temp_diameter);


	n_generated++;
	if (n_generated >= nPores[current_psd_index])
	{
		current_psd_index++;
		n_generated = 0;
	}

	//assign pore to a cell
	xCell = (int)floor(temp_pore[0] / cellSize);
	yCell = (int)floor(temp_pore[1] / cellSize);
	zCell = (int)floor(temp_pore[2] / cellSize);
	cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
	cellCoords[cellIndex].push_back(temp_pore);

	cout << "Pore generation progress:\n";
	cout << setprecision(2) << std::fixed;

	//loop for generating the rest of the LARGE pores
	int flag;
	double progress;
	for (int i = 1; i < large_pores; i++)
	{
		temp_diameter = diameter_pore_gen(current_psd_index, settings); //generate pore size

		//V1: assign pores to cells
		flag = 1; //flag is 1 (overlap)to enter the loop, but reset to zero(no overlap) straight after and raised in the overlap loop
		while (flag > 0)
		{
			flag = 0; //flag reset
			//take a random coordinate and calculate its location in the cell grid
			temp_pore.clear();
			double failsafe = 0.999999999;
			temp_pore = { ((double)rand() / RAND_MAX) * settings.domain * failsafe, ((double)rand() / RAND_MAX) * settings.domain * failsafe, ((double)rand() / RAND_MAX) * settings.domain * failsafe };
			temp_pore.push_back(temp_diameter);
			

			xCell = (int)floor(temp_pore[0] / cellSize);
			yCell = (int)floor(temp_pore[1] / cellSize);
			zCell = (int)floor(temp_pore[2] / cellSize);
			//cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
			//cellCoords[cellIndex].push_back(temp_pore);

			//run search in 25 cells around the cell of interest
			int cellIndexSearch;
			double distance2;
			for (int iCell = max(0, xCell - 2); iCell < min(nCells, xCell + 3); iCell++)
			{
				for (int jCell = max(0, yCell - 2); jCell < min(nCells, yCell + 3); jCell++)
				{
					for (int kCell = max(0, zCell - 2); kCell < min(nCells, zCell + 3); kCell++)
					{
						cellIndexSearch = pow(nCells, 2)*iCell + nCells * jCell + kCell;
						for (int pore_in_cell = 0; pore_in_cell < cellCoords[cellIndexSearch].size(); pore_in_cell++)
						{
							distance2 = pow((temp_pore[0] - cellCoords[cellIndexSearch][pore_in_cell][0]), 2) + pow((temp_pore[1] - cellCoords[cellIndexSearch][pore_in_cell][1]), 2) + pow((temp_pore[2] - cellCoords[cellIndexSearch][pore_in_cell][2]), 2);
							if (distance2 <= pow(cellCoords[cellIndexSearch][pore_in_cell][3] / 2 + temp_pore[3] / 2, 2))
							{
								flag = 1;
							}
							if (flag > 0) { break; }
						}
						if (flag > 0) { break; }
					}
					if (flag > 0) { break; }
				}
				if (flag > 0) { break; }
			}
		}


		cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
		cellCoords[cellIndex].push_back(temp_pore);
		network.pores.push_back(temp_pore);
		//pore HAS TO BE this diameter, assign now
		network.diameter_pore.push_back(temp_diameter);

		n_generated++;
		if (n_generated >= nPores[current_psd_index])
		{
			current_psd_index++;
			n_generated = 0;
		}

		if (i % 10000 == 0)
		{
			progress = 100 * (i + 1) / (double)total_pn;
			cout << "\r" << progress << "%" << flush;
		}

	}

	//loop for generating SMALL pores
	for (int i = large_pores; i < total_pn; i++)
	{
		temp_diameter = diameter_pore_gen(current_psd_index, settings); //generate pore size

		//V1: assign pores to cells
		flag = 1; //flag is 1 (overlap)to enter the loop, but reset to zero(no overlap) straight after and raised in the overlap loop
		while (flag > 0)
		{
			flag = 0; //flag reset
			//take a random coordinate and calculate its location in the cell grid
			temp_pore.clear();
			double failsafe = 0.999999999;
			temp_pore = { ((double)rand() / RAND_MAX) * settings.domain * failsafe, ((double)rand() / RAND_MAX) * settings.domain * failsafe, ((double)rand() / RAND_MAX) * settings.domain * failsafe };
			temp_pore.push_back(temp_diameter);
			

			//TEST LARGE GRID FIRST
			xCell = (int)floor(temp_pore[0] / cellSize);
			yCell = (int)floor(temp_pore[1] / cellSize);
			zCell = (int)floor(temp_pore[2] / cellSize);
			//cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
			//cellCoords[cellIndex].push_back(temp_pore);

			//run search in 25 cells around the cell of interest
			int cellIndexSearch;
			double distance2;
			for (int iCell = max(0, xCell - 2); iCell < min(nCells, xCell + 3); iCell++)
			{
				for (int jCell = max(0, yCell - 2); jCell < min(nCells, yCell + 3); jCell++)
				{
					for (int kCell = max(0, zCell - 2); kCell < min(nCells, zCell + 3); kCell++)
					{
						cellIndexSearch = pow(nCells, 2)*iCell + nCells * jCell + kCell;
						for (int pore_in_cell = 0; pore_in_cell < cellCoords[cellIndexSearch].size(); pore_in_cell++)
						{
							distance2 = pow((temp_pore[0] - cellCoords[cellIndexSearch][pore_in_cell][0]), 2) + pow((temp_pore[1] - cellCoords[cellIndexSearch][pore_in_cell][1]), 2) + pow((temp_pore[2] - cellCoords[cellIndexSearch][pore_in_cell][2]), 2);
							if (distance2 <= pow(cellCoords[cellIndexSearch][pore_in_cell][3] / 2 + temp_pore[3] / 2, 2))
							{
								flag = 1;
							}
							if (flag > 0) { break; }
						}
						if (flag > 0) { break; }
					}
					if (flag > 0) { break; }
				}
				if (flag > 0) { break; }
			}
			//if fist check exits with flag > 0, it will exit instantly in the next loop, loss of efficiency is minimal

			//TEST SMALL GRID
			xCellSmall = (int)floor(temp_pore[0] / cellSizeSmall);
			yCellSmall = (int)floor(temp_pore[1] / cellSizeSmall);
			zCellSmall = (int)floor(temp_pore[2] / cellSizeSmall);
			//cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
			//cellCoords[cellIndex].push_back(temp_pore);

			//run search in 25 cells around the cell of interest
			/*int cellIndexSearch;
			double distance2;*/
			for (int iCell = max(0, xCellSmall - 2); iCell < min(nCellsSmall, xCellSmall + 3); iCell++)
			{
				for (int jCell = max(0, yCellSmall - 2); jCell < min(nCellsSmall, yCellSmall + 3); jCell++)
				{
					for (int kCell = max(0, zCellSmall - 2); kCell < min(nCellsSmall, zCellSmall + 3); kCell++)
					{
						cellIndexSearch = pow(nCellsSmall, 2)*iCell + nCellsSmall * jCell + kCell;
						for (int pore_in_cell = 0; pore_in_cell < cellCoordsSmall[cellIndexSearch].size(); pore_in_cell++)
						{
							distance2 = pow((temp_pore[0] - cellCoordsSmall[cellIndexSearch][pore_in_cell][0]), 2) + pow((temp_pore[1] - cellCoordsSmall[cellIndexSearch][pore_in_cell][1]), 2) + pow((temp_pore[2] - cellCoordsSmall[cellIndexSearch][pore_in_cell][2]), 2);
							if (distance2 <= pow(cellCoordsSmall[cellIndexSearch][pore_in_cell][3] / 2 + temp_pore[3] / 2, 2))
							{
								flag = 1;
							}
							if (flag > 0) { break; }
						}
						if (flag > 0) { break; }
					}
					if (flag > 0) { break; }
				}
				if (flag > 0) { break; }
			}
		}

		//push back generated pore into the small grid storage
		cellIndexSmall = pow(nCellsSmall, 2)*xCellSmall + nCellsSmall * yCellSmall + zCellSmall;
		cellCoordsSmall[cellIndexSmall].push_back(temp_pore);
		network.pores.push_back(temp_pore);
		//pore HAS TO BE this diameter, assign now
		network.diameter_pore.push_back(temp_diameter);

		n_generated++;
		if (n_generated >= nPores[current_psd_index])
		{
			current_psd_index++;
			n_generated = 0;
		}

		if (i % 10000 == 0)
		{
			progress = 100 * (i + 1) / (double)total_pn;
			cout << "\r" << progress << "%" << flush;
		}

	}
	progress = 100;
	cout << "\r" << progress << "%" << flush;

	cout << "\nNumber of pores generated: " << network.pores.size() << endl;

	chrono::steady_clock::time_point end = chrono::steady_clock::now();
	string elapsedTime = "Elapsed time: " + to_string(chrono::duration_cast<std::chrono::seconds>(end - begin).count()) + " [s]";
	cout << elapsedTime << endl;
}

void classPNM::classMethods::classMethodsGeneration::throatsGeneration(classPNM::classNetwork& network, classPNM::classSettings& settings)
{
	chrono::steady_clock::time_point begin = chrono::steady_clock::now();

	//vector of throat connectivities {start pore index; end pore index; throat diameter}
	network.throats.clear();
	network.diameter_throat.clear();

	double running_z = 0;
	int n_nearest = 10; //randomly connect to Z/2 out of nearest N pores, to provide some variability
	srand(settings.seed);

	//converting string into double PSD
	vector<double> temp;
	settings.poreSizeDist.clear();
	for (int i = 0; i < settings.poreSizeDist_raw.size(); i++)
	{
		temp.clear();
		for (int j = 0; j < settings.poreSizeDist_raw[i].size(); j++)
		{
			if (j == 0) { temp.push_back(stod(settings.poreSizeDist_raw[i][j])*1e-9); }
			else { temp.push_back(stod(settings.poreSizeDist_raw[i][j])); }
		}
		settings.poreSizeDist.push_back(temp);
	}

	//calculate number of pores per each poe size with given domain
	vector<int> nPores = poresPsize_(settings.poreSizeDist, settings.coordination, settings.domain);

	//variable definition
	int total_pn, current_psd_index;
	double expected_porosity;

	total_pn = 0;
	//calculate total pore number
	for (int i = 0; i < nPores.size(); i++)
	{
		if (!(nPores[i] < 1))
		{
			total_pn += nPores[i];
		}
	}

	//define which pore size we are starting generation from
	current_psd_index = 0;
	while (nPores[current_psd_index] < 1)
	{
		current_psd_index++;
	}

	//define how many pores we want to generate as a "large" grid:
	double fraction = 0.0;
	int cut_off = current_psd_index; //everything up to but NOT including THIS entry into PSD to be generated as "large" grid
	int large_pores = 0; //number of large pores to be generated
	while (fraction < 0.001)
	{
		fraction += (double)nPores[cut_off] / total_pn;
		large_pores += nPores[cut_off];
		cut_off++;

	};
	cout << "Cut-off: " << cut_off << endl;


	//assign large pores into cells:
	//divide donain into cells, just larger than biggest pore radius
	double cellSize = network.pores[0][3] / 2;
	int nCells = (int)floor(settings.domain / cellSize);
	cellSize = settings.domain / nCells;

	vector<vector<int>> cellCoords(pow(nCells, 3));
	int xCell, yCell, zCell, cellIndex;
	for (int i = 0; i < large_pores; i++)
	{
		xCell = (int)floor(network.pores[i][0] / cellSize);
		yCell = (int)floor(network.pores[i][1] / cellSize);
		zCell = (int)floor(network.pores[i][2] / cellSize);
		cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
		cellCoords[cellIndex].push_back(i);
	}
	int assignment_check = 0; //check if all pors have been assigned to cells correctly
	for (int i = 0; i < cellCoords.size(); i++)
	{
		assignment_check += cellCoords[i].size();
	}

	//assign SMALL pores into cells:
	//divide donain into cells:
	//divide domain into "small" cells, based on the cut-off pore size
	double cellSizeSmall = settings.poreSizeDist[cut_off][0] / 2;
	int nCellsSmall = (int)floor(settings.domain / cellSizeSmall);
	cellSizeSmall = settings.domain / nCellsSmall;

	vector<vector<int>> cellCoordsSmall(pow(nCellsSmall, 3));

	int xCellSmall, yCellSmall, zCellSmall, cellIndexSmall;
	for (int i = large_pores; i < total_pn; i++)
	{
		xCellSmall = (int)floor(network.pores[i][0] / cellSizeSmall);
		yCellSmall = (int)floor(network.pores[i][1] / cellSizeSmall);
		zCellSmall = (int)floor(network.pores[i][2] / cellSizeSmall);
		cellIndexSmall = pow(nCellsSmall, 2)*xCellSmall + nCellsSmall * yCellSmall + zCellSmall;
		cellCoordsSmall[cellIndexSmall].push_back(i);
	}
	for (int i = 0; i < cellCoordsSmall.size(); i++)
	{
		assignment_check += cellCoordsSmall[i].size();
	}
	cout << assignment_check << endl;

	//vector for storing existing connections for checking at generation
	vector<vector<int>> existing_conns(total_pn);
	//vector<int> temp_throats; //temporary storage of connectivity to be pushed back at the end of the loop

	vector <vector<double>> distance_array; //array for storing sorted distance array for closest pores, {distance; pore index}

	cout << "Throat generation progress:\n";
	cout << setprecision(2) << std::fixed;
	double progress;


	//main loop for creacting Z/2 throats for each pore
	for (int i = 0; i < total_pn; i++)
	{
		distance_array.clear();
		int large_grid_search_radius = 2;
		int small_grid_search_radius = 4;
		//search for closest pores within nearest 5 cells of the LARGE grid

		xCell = (int)floor(network.pores[i][0] / cellSize);
		yCell = (int)floor(network.pores[i][1] / cellSize);
		zCell = (int)floor(network.pores[i][2] / cellSize);
		cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
		int cellIndexSearch;
		double distance2;
		for (int iCell = max(0, xCell - large_grid_search_radius); iCell < min(nCells, xCell + large_grid_search_radius + 1); iCell++)
		{
			for (int jCell = max(0, yCell - large_grid_search_radius); jCell < min(nCells, yCell + large_grid_search_radius + 1); jCell++)
			{
				for (int kCell = max(0, zCell - large_grid_search_radius); kCell < min(nCells, zCell + large_grid_search_radius + 1); kCell++)
				{
					cellIndexSearch = pow(nCells, 2)*iCell + nCells * jCell + kCell;
					for (int pore_in_cell = 0; pore_in_cell < cellCoords[cellIndexSearch].size(); pore_in_cell++)
					{
						//construct the distance array for the pores within the searched cells 
						distance_array.push_back({ (pow((network.pores[i][0] - network.pores[cellCoords[cellIndexSearch][pore_in_cell]][0]), 2) + pow((network.pores[i][1] - network.pores[cellCoords[cellIndexSearch][pore_in_cell]][1]), 2) + pow((network.pores[i][2] - network.pores[cellCoords[cellIndexSearch][pore_in_cell]][2]), 2)), (double)cellCoords[cellIndexSearch][pore_in_cell] });
					}
				}
			}
		}

		//search for closest pores within nearest x cells of the SMALL grid
		xCellSmall = (int)floor(network.pores[i][0] / cellSizeSmall);
		yCellSmall = (int)floor(network.pores[i][1] / cellSizeSmall);
		zCellSmall = (int)floor(network.pores[i][2] / cellSizeSmall);
		cellIndexSmall = pow(nCellsSmall, 2)*xCellSmall + nCellsSmall * yCellSmall + zCellSmall;
		//int cellIndexSearchSmall;
		//double distance2;
		for (int iCell = max(0, xCellSmall - small_grid_search_radius); iCell < min(nCellsSmall, xCellSmall + small_grid_search_radius + 1); iCell++)
		{
			for (int jCell = max(0, yCellSmall - small_grid_search_radius); jCell < min(nCellsSmall, yCellSmall + small_grid_search_radius + 1); jCell++)
			{
				for (int kCell = max(0, zCellSmall - small_grid_search_radius); kCell < min(nCellsSmall, zCellSmall + small_grid_search_radius + 1); kCell++)
				{
					cellIndexSearch = pow(nCellsSmall, 2)*iCell + nCellsSmall * jCell + kCell;
					for (int pore_in_cell = 0; pore_in_cell < cellCoordsSmall[cellIndexSearch].size(); pore_in_cell++)
					{
						//construct the distance array for the pores within the searched cells 
						distance_array.push_back({ (pow((network.pores[i][0] - network.pores[cellCoordsSmall[cellIndexSearch][pore_in_cell]][0]), 2) + pow((network.pores[i][1] - network.pores[cellCoordsSmall[cellIndexSearch][pore_in_cell]][1]), 2) + pow((network.pores[i][2] - network.pores[cellCoordsSmall[cellIndexSearch][pore_in_cell]][2]), 2)), (double)cellCoordsSmall[cellIndexSearch][pore_in_cell] });
					}
				}
			}
		}
		//sorting the distance array 
		if (distance_array.size() < n_nearest)
		{
			cout << "Not enough neighbours for pore " << i << ", cannot proceed generation!" << endl;
			exit(1);
		}

		sort(distance_array.begin(), distance_array.end(), sortColumn0);

		int connection, flag; //which closest pore do we try to connect to?
		if (running_z < settings.coordination)
		{
			for (int j = 0; j < ceil(settings.coordination / 2); j++)
			{
				//check for existing connection on the pore of interest
				flag = 0;
				connection = rand() % n_nearest + 1; //desired pore to connect to from the list of closest pores

				while (flag >= 0)
				{
					if (existing_conns[i].empty()) { flag = -1; } //if there are no existing connections on pore i, skip
					else
					{
						for (int k = 0; k < existing_conns[i].size(); k++)
						{
							if (existing_conns[i][k] == distance_array[connection][1])
							{
								connection = rand() % n_nearest + 1;
								flag = 0;
								break;
							}
							flag = -1;
						}
					}
				}


				//update the connection array adn the diameter array
				if (i < distance_array[connection][1])
				{
					network.throats.push_back({ (double)i, distance_array[connection][1], network.pores[i][3] });
				}
				else
				{
					network.throats.push_back({ distance_array[connection][1], (double)i, network.pores[i][3] });
				}
				network.diameter_throat.push_back(network.pores[i][3]);

				//update the existing connection array with the new connection
				existing_conns[i].push_back((int)distance_array[connection][1]);
				existing_conns[(int)distance_array[connection][1]].push_back(i);
			}
		}
		else
		{
			for (int j = 0; j < floor(settings.coordination / 2); j++)
			{
				//check for existing connection on the pore of interest
				flag = 0;
				connection = rand() % n_nearest + 1; //desired pore to connect to from the list of closest pores

				while (flag >= 0)
				{
					if (existing_conns[i].empty()) { flag = -1; } //if there are no existing connections on pore i, skip
					else
					{
						for (int k = 0; k < existing_conns[i].size(); k++)
						{
							if (existing_conns[i][k] == distance_array[connection][1])
							{
								connection = rand() % n_nearest + 1;
								flag = 0;
								break;
							}
							flag = -1;
						}
					}
				}


				//update the connection array adn the diameter array
				if (i < distance_array[connection][1])
				{
					network.throats.push_back({ (double)i, distance_array[connection][1], network.pores[i][3] });
				}
				else
				{
					network.throats.push_back({ distance_array[connection][1], (double)i, network.pores[i][3] });
				}
				network.diameter_throat.push_back(network.pores[i][3]);

				//update the existing connection array with the new connection
				existing_conns[i].push_back((int)distance_array[connection][1]);
				existing_conns[(int)distance_array[connection][1]].push_back(i);
			}
		}
		running_z = 2 * (double)network.throats.size() / (double)(i + 1);

		if (i % 1000 == 0)
		{
			progress = 100 * (i + 1) / (double)total_pn;
			cout << "\r" << progress << "%" << flush;
		}
	}
	progress = 100;
	cout << "\r" << progress << "%" << flush;
	cout << "\nNumber of throats generated: " << network.throats.size() << "\n";


	chrono::steady_clock::time_point end = chrono::steady_clock::now();
	string elapsedTime = "Elapsed time: " + to_string(chrono::duration_cast<std::chrono::seconds>(end - begin).count()) + " [s]";
	cout << elapsedTime << endl;
}

void classPNM::classMethods::classMethodsGeneration::throatsGeneration2(classPNM::classNetwork& network, classPNM::classSettings& settings)
{
	chrono::steady_clock::time_point begin = chrono::steady_clock::now();

	//vector of throat connectivities {start pore index; end pore index; throat diameter}
	network.throats.clear();
	network.diameter_throat.clear();

	double running_z = 0;
	int n_nearest = 10; //randomly connect to Z/2 out of nearest N pores, to provide some variability
	srand(settings.seed);

	//converting string into double PSD
	vector<double> temp;
	settings.poreSizeDist.clear();
	for (int i = 0; i < settings.poreSizeDist_raw.size(); i++)
	{
		temp.clear();
		for (int j = 0; j < settings.poreSizeDist_raw[i].size(); j++)
		{
			if (j == 0) { temp.push_back(stod(settings.poreSizeDist_raw[i][j])*1e-9); }
			else { temp.push_back(stod(settings.poreSizeDist_raw[i][j])); }
		}
		settings.poreSizeDist.push_back(temp);
	}

	//calculate number of pores per each poe size with given domain
	vector<int> nPores = poresPsize_(settings.poreSizeDist, settings.coordination, settings.domain);

	//variable definition
	int total_pn, current_psd_index;
	double expected_porosity;

	total_pn = 0;
	//calculate total pore number
	for (int i = 0; i < nPores.size(); i++)
	{
		if (!(nPores[i] < 1))
		{
			total_pn += nPores[i];
		}
	}

	//define which pore size we are starting generation from
	current_psd_index = 0;
	while (nPores[current_psd_index] < 1)
	{
		current_psd_index++;
	}

	//define how many pores we want to generate as a "large" grid:
	double fraction = 0.0;
	int cut_off = current_psd_index; //everything up to but NOT including THIS entry into PSD to be generated as "large" grid
	int large_pores = 0; //number of large pores to be generated
	while (fraction < 0.001)
	{
		fraction += (double)nPores[cut_off] / total_pn;
		large_pores += nPores[cut_off];
		cut_off++;

	};
	//cout << "Cut-off: " << cut_off << endl;


	//assign large pores into cells:
	//divide donain into cells, just larger than biggest pore radius
	double cellSize = network.pores[0][3] / 2;
	int nCells = (int)floor(settings.domain / cellSize);
	cellSize = settings.domain / nCells;

	vector<vector<int>> cellCoords(pow(nCells, 3));
	int xCell, yCell, zCell, cellIndex;
	for (int i = 0; i < large_pores; i++)
	{
		xCell = (int)floor(network.pores[i][0] / cellSize);
		yCell = (int)floor(network.pores[i][1] / cellSize);
		zCell = (int)floor(network.pores[i][2] / cellSize);
		cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
		cellCoords[cellIndex].push_back(i);
	}
	int assignment_check = 0; //check if all pors have been assigned to cells correctly
	for (int i = 0; i < cellCoords.size(); i++)
	{
		assignment_check += cellCoords[i].size();
	}

	//assign SMALL pores into cells:
	//divide donain into cells:
	//divide domain into "small" cells, based on the cut-off pore size
	double cellSizeSmall = settings.poreSizeDist[cut_off][0] / 2;
	int nCellsSmall = (int)floor(settings.domain / cellSizeSmall);
	cellSizeSmall = settings.domain / nCellsSmall;

	vector<vector<int>> cellCoordsSmall(pow(nCellsSmall, 3));

	int xCellSmall, yCellSmall, zCellSmall, cellIndexSmall;
	for (int i = large_pores; i < total_pn; i++)
	{
		xCellSmall = (int)floor(network.pores[i][0] / cellSizeSmall);
		yCellSmall = (int)floor(network.pores[i][1] / cellSizeSmall);
		zCellSmall = (int)floor(network.pores[i][2] / cellSizeSmall);
		cellIndexSmall = pow(nCellsSmall, 2)*xCellSmall + nCellsSmall * yCellSmall + zCellSmall;
		cellCoordsSmall[cellIndexSmall].push_back(i);
	}
	for (int i = 0; i < cellCoordsSmall.size(); i++)
	{
		assignment_check += cellCoordsSmall[i].size();
	}
	//cout << assignment_check << endl;

	//vector for storing existing connections for checking at generation
	vector<vector<int>> existing_conns(total_pn);
	//vector<int> temp_throats; //temporary storage of connectivity to be pushed back at the end of the loop

	vector <vector<double>> distance_array; //array for storing sorted distance array for closest pores, {distance; pore index}
	vector <vector<double>> distance_array_temp; //array for storing sorted distance array for closest pores, {distance; pore index}

	cout << "Throat generation progress:\n";
	cout << setprecision(2) << std::fixed;
	double progress;

	int existing_connections = 0;

	//variables for distance and direction calc
	double distance2;
	double cosine2;
	//variables declared outside of the loop 
	int large_grid_search_radius = 2;
	int small_grid_search_radius = 4;
	int cellIndexSearch;
	double failsafe = 0.999999999;
	double anisotropy_pow = 1 / settings.anisotropy;
	//main loop for creacting Z/2 throats for each pore
	for (int i = 0; i < total_pn; i++)
	{
		distance_array.clear();
		n_nearest = 10;
		
		//search for closest pores within nearest 5 cells of the LARGE grid

		xCell = (int)floor(network.pores[i][0] / cellSize);
		yCell = (int)floor(network.pores[i][1] / cellSize);
		zCell = (int)floor(network.pores[i][2] / cellSize);
		cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
		
		for (int iCell = max(0, xCell - large_grid_search_radius); iCell < min(nCells, xCell + large_grid_search_radius + 1); iCell++)
		{
			for (int jCell = max(0, yCell - large_grid_search_radius); jCell < min(nCells, yCell + large_grid_search_radius + 1); jCell++)
			{
				for (int kCell = max(0, zCell - large_grid_search_radius); kCell < min(nCells, zCell + large_grid_search_radius + 1); kCell++)
				{
					cellIndexSearch = pow(nCells, 2)*iCell + nCells * jCell + kCell;
					for (int pore_in_cell = 0; pore_in_cell < cellCoords[cellIndexSearch].size(); pore_in_cell++)
					{
						//construct the distance array for the pores within the searched cells 
						//calculate distance and diraction with retation to XY plane 
						distance2 = (pow((network.pores[i][0] - network.pores[cellCoords[cellIndexSearch][pore_in_cell]][0]), 2) + pow((network.pores[i][1] - network.pores[cellCoords[cellIndexSearch][pore_in_cell]][1]), 2) + pow((network.pores[i][2] - network.pores[cellCoords[cellIndexSearch][pore_in_cell]][2]), 2));
						cosine2 = pow((network.pores[i][2] - network.pores[cellCoords[cellIndexSearch][pore_in_cell]][2]), 2) / (pow((network.pores[i][0] - network.pores[cellCoords[cellIndexSearch][pore_in_cell]][0]), 2) + pow((network.pores[i][1] - network.pores[cellCoords[cellIndexSearch][pore_in_cell]][1]), 2));
						distance_array.push_back({ distance2, (double)cellCoords[cellIndexSearch][pore_in_cell], cosine2 });
					}
				}
			}
		}

		//search for closest pores within nearest x cells of the SMALL grid
		xCellSmall = (int)floor(network.pores[i][0] / cellSizeSmall);
		yCellSmall = (int)floor(network.pores[i][1] / cellSizeSmall);
		zCellSmall = (int)floor(network.pores[i][2] / cellSizeSmall);
		cellIndexSmall = pow(nCellsSmall, 2)*xCellSmall + nCellsSmall * yCellSmall + zCellSmall;
		//int cellIndexSearchSmall;
		for (int iCell = max(0, xCellSmall - small_grid_search_radius); iCell < min(nCellsSmall, xCellSmall + small_grid_search_radius + 1); iCell++)
		{
			for (int jCell = max(0, yCellSmall - small_grid_search_radius); jCell < min(nCellsSmall, yCellSmall + small_grid_search_radius + 1); jCell++)
			{
				for (int kCell = max(0, zCellSmall - small_grid_search_radius); kCell < min(nCellsSmall, zCellSmall + small_grid_search_radius + 1); kCell++)
				{
					cellIndexSearch = pow(nCellsSmall, 2)*iCell + nCellsSmall * jCell + kCell;
					for (int pore_in_cell = 0; pore_in_cell < cellCoordsSmall[cellIndexSearch].size(); pore_in_cell++)
					{
						//construct the distance array for the pores within the searched cells 
						distance2 = (pow((network.pores[i][0] - network.pores[cellCoordsSmall[cellIndexSearch][pore_in_cell]][0]), 2) + pow((network.pores[i][1] - network.pores[cellCoordsSmall[cellIndexSearch][pore_in_cell]][1]), 2) + pow((network.pores[i][2] - network.pores[cellCoordsSmall[cellIndexSearch][pore_in_cell]][2]), 2));
						cosine2 = pow((network.pores[i][2] - network.pores[cellCoordsSmall[cellIndexSearch][pore_in_cell]][2]), 2) / (pow((network.pores[i][0] - network.pores[cellCoordsSmall[cellIndexSearch][pore_in_cell]][0]), 2) + pow((network.pores[i][1] - network.pores[cellCoordsSmall[cellIndexSearch][pore_in_cell]][1]), 2));
						distance_array.push_back({ distance2 , (double)cellCoordsSmall[cellIndexSearch][pore_in_cell], cosine2 });
					}
				}
			}
		}
		//sorting the distance array 
		if (distance_array.size() <= n_nearest)
		{
			//cout << "Not enough neighbours for pore " << i << ", cannot proceed generation!" << endl;
			//exit(1);
			n_nearest = distance_array.size() - 1;
		}
		//sort by distance
		sort(distance_array.begin(), distance_array.end(), sortColumn0);
		//cut to n_nearest + 1 length (to account for the fact that this pore found itself as a closest neigbour)
		distance_array.resize(n_nearest + 1);
		//sort by direction
		sort(distance_array.begin(), distance_array.end(), sortColumn2);

		int connection, flag; //which closest pore do we try to connect to?
		if (running_z < settings.coordination)
		{
			for (int j = 0; j < ceil(settings.coordination / 2); j++)
			{
				//check for existing connection on the pore of interest
				flag = 0;
				//connection = rand() % n_nearest + 1; //desired pore to connect to from the list of closest pores, not including itself
				connection = (int)(pow((double)rand() / RAND_MAX, anisotropy_pow) * failsafe * n_nearest) + 1; //produce value [0, 1), change to desired range, truncate to int and shift 1 up to avoid connection with itself

				while (flag >= 0)
				{
					if (existing_conns[i].empty()) { flag = -1; } //if there are no existing connections on pore i, skip
					else
					{
						for (int k = 0; k < existing_conns[i].size(); k++)
						{
							if (existing_conns[i][k] == distance_array[connection][1])
							{
								connection = (int)(pow((double)rand() / RAND_MAX, anisotropy_pow) * failsafe * n_nearest) + 1; //produce value [0, 1), change to desired range, truncate to int and shift 1 up to avoid connection with itself
								flag = 0;
								existing_connections++;
								break;
							}
							flag = -1;
						}
					}
				}


				//update the connection array adn the diameter array
				if (i < distance_array[connection][1])
				{
					network.throats.push_back({ (double)i, distance_array[connection][1], network.pores[i][3] });
				}
				else
				{
					network.throats.push_back({ distance_array[connection][1], (double)i, network.pores[i][3] });
				}
				network.diameter_throat.push_back(network.pores[i][3]);

				//update the existing connection array with the new connection
				existing_conns[i].push_back((int)distance_array[connection][1]);
				existing_conns[(int)distance_array[connection][1]].push_back(i);
			}
		}
		else
		{
			for (int j = 0; j < floor(settings.coordination / 2); j++)
			{
				//check for existing connection on the pore of interest
				flag = 0;
				connection = (int)(pow((double)rand() / RAND_MAX, anisotropy_pow) * failsafe * n_nearest) + 1; //produce value [0, 1), change to desired range, truncate to int and shift 1 up to avoid connection with itself

				while (flag >= 0)
				{
					if (existing_conns[i].empty()) { flag = -1; } //if there are no existing connections on pore i, skip
					else
					{
						for (int k = 0; k < existing_conns[i].size(); k++)
						{
							if (existing_conns[i][k] == distance_array[connection][1])
							{
								connection = (int)(pow((double)rand() / RAND_MAX, anisotropy_pow) * failsafe * n_nearest) + 1; //produce value [0, 1), change to desired range, truncate to int and shift 1 up to avoid connection with itself
								flag = 0;
								existing_connections++;
								break;
							}
							flag = -1;
						}
					}
				}


				//update the connection array adn the diameter array
				if (i < distance_array[connection][1])
				{
					network.throats.push_back({ (double)i, distance_array[connection][1], network.pores[i][3] });
				}
				else
				{
					network.throats.push_back({ distance_array[connection][1], (double)i, network.pores[i][3] });
				}
				network.diameter_throat.push_back(network.pores[i][3]);

				//update the existing connection array with the new connection
				existing_conns[i].push_back((int)distance_array[connection][1]);
				existing_conns[(int)distance_array[connection][1]].push_back(i);
			}
		}
		running_z = 2 * (double)network.throats.size() / (double)(i + 1);

		if (i % 1000 == 0)
		{
			progress = 100 * (i + 1) / (double)total_pn;
			cout << "\r" << progress << "%" << flush;
		}
	}
	progress = 100;
	cout << "\r" << progress << "%" << flush;
	cout << "\nNumber of throats generated: " << network.throats.size() << "\n";
	//cout << "Existing connections encountered: " << existing_connections << endl;


	chrono::steady_clock::time_point end = chrono::steady_clock::now();
	string elapsedTime = "Elapsed time: " + to_string(chrono::duration_cast<std::chrono::seconds>(end - begin).count()) + " [s]";
	cout << elapsedTime << endl;
}

void classPNM::classMethods::classMethodsGeneration::connectedComponents(classPNM::classNetwork& network)
{
	using namespace boost;
	
		int Nt = network.throats.size();
		int Np = network.pores.size();
		typedef adjacency_list <vecS, vecS, undirectedS> Graph;

		//using example code from boost library for definiton of graph clusters
		Graph G;
		for (int i = 0; i < Nt; i++)
		{
			add_edge(network.throats[i][0], network.throats[i][1], G);
		}

		cout << "Searching for disconnected clusters..." << endl;
		chrono::steady_clock::time_point begin = chrono::steady_clock::now();

		std::vector<int> component(num_vertices(G));
		int num = connected_components(G, &component[0]);

		chrono::steady_clock::time_point end = chrono::steady_clock::now();
		string elapsedTime = "Elapsed time: " + to_string(chrono::duration_cast<std::chrono::seconds>(end - begin).count()) + " [s]";
		//cout << elapsedTime << endl;

		//cout << "Rest of the health code:" << endl;
		begin = chrono::steady_clock::now();

		vector<int> clusters(num, 0);

		//define how many pores are in each cluster 
		for (int i = 0; i < component.size(); i++)
		{
			clusters[component[i]]++;
		}

		//which cluster is NOT removed? (he biggest one)
		int max = std::max_element(clusters.begin(), clusters.end()) - clusters.begin();

		//remove all pore entries apart from the BIGGEST cluster
		vector<int> old_indexes; //exentully stores a vector that contains the old index of a pore (but at the new position): 
									//for example - if network.pores[5] is to be removed, old_indexes[5] = 6, because 6th pore was shifted one back
		vector<int> new_indexes(Np); //constructed from the old_indices, shows what new index value each pore has. continiuing example above, new_indexes[6] = 5, because what previously was 6th pore is now 5th etc. 
										//Elements for deleted pores left undefined since they are never used 
		int j;
		vector<vector<double>> temp_pores; //temporary containers for faster deletion
		vector<double> temp_diameter_pore;
		for (int i = 0; i < Np; i++)
		{
			if (component[i] == max)
			{
				old_indexes.push_back(i);//pass old index to the storage
				temp_pores.push_back(network.pores[i]);
				temp_diameter_pore.push_back(network.diameter_pore[i]);
			}
		}
		network.pores = temp_pores;
		network.diameter_pore = temp_diameter_pore;
		temp_pores.clear();
		temp_diameter_pore.clear();

		//fill new_indexes
		for (int i = 0; i < old_indexes.size(); i++)
		{
			new_indexes[old_indexes[i]] = i;
		}

		//remove all throat entries for which Start or End pores are NOT in the BIGGEST cluster
		//if pores are in the biggest cluster, change indexes from old to new using new_indexes vector.
		vector<vector<double>> temp_throats; //temporary containers for faster deletion
		vector<double> temp_diameter_throat;
		vector<double> temp;
		for (int j = 0; j < Nt; j++)
		{
			temp.clear();
			if ((component[network.throats[j][0]] == max) && (component[network.throats[j][1]] == max)) //only if both start and end are a part of the BIGGEST cluster
			{
				temp_diameter_throat.push_back(network.diameter_throat[j]);

				temp.push_back(new_indexes[network.throats[j][0]]);
				temp.push_back(new_indexes[network.throats[j][1]]);
				temp.push_back(network.diameter_throat[j]);

				temp_throats.push_back(temp);
			}
		}
		network.throats = temp_throats;
		network.diameter_throat = temp_diameter_throat;
		temp_throats.clear();
		temp_diameter_throat.clear();

		end = chrono::steady_clock::now();
		elapsedTime = "Elapsed time: " + to_string(chrono::duration_cast<std::chrono::seconds>(end - begin).count()) + " [s]";
		//cout << elapsedTime << endl;

		cout << "Final number of pores: " << network.pores.size() << endl;
		cout << "Final number of throats: " << network.throats.size() << endl;

}