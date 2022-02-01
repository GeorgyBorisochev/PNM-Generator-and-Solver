#include "classPNM.h"



///
///ClassPNM functions
///

//classPNM classPNM::poresGeneration(classPNM PNM, classUtils utils)
//{
//
//	PNM.poreSizeDist.clear();
//	PNM.network.pores.clear();
//	PNM.network.diameter_pore.clear();
//
//	//variable for holding the array of pore coordinates
//	vector<vector<double>> pores;
//
//	//converting string into double PSD
//	vector<double> temp;
//	for (int i = 0; i < PNM.poreSizeDist_raw.size(); i++)
//	{
//		temp.clear();
//		for (int j = 0; j < PNM.poreSizeDist_raw[i].size(); j++)
//		{
//			if (j == 0) { temp.push_back(stod(PNM.poreSizeDist_raw[i][j])*1e-9); }
//			else { temp.push_back(stod(PNM.poreSizeDist_raw[i][j])); }
//		}
//		PNM.poreSizeDist.push_back(temp);
//	}
//
//	//calculate number of pores per each poe size with given domain
//	vector<int> nPores = PNM.poresPsize(PNM);
//
//	//variable definition
//	int total_pn, current_psd_index, nCells, xCell, yCell, zCell, cellIndex;
//
//	double expected_porosity, cellSize, pores_per_cell;
//
//	total_pn = 0;
//	expected_porosity = 0;
//	//calculate total pore number
//	for (int i = 0; i < nPores.size(); i++)
//	{
//		if (!(nPores[i] < 1))
//		{
//			total_pn += nPores[i];
//			expected_porosity += PNM.poreSizeDist[i][1];
//		}
//	}
//	PNM.expectedP = expected_porosity;
//	cout << "Expected network porosity is: " << expected_porosity*100.0 << " %\n";
//
//	//define which pore size we are starting generation from
//	current_psd_index = 0;
//	while (nPores[current_psd_index] < 1)
//	{
//		current_psd_index++;
//	}
//
//	//divide donain into cells, just larger than biggest pore radius
//	cellSize = PNM.poreSizeDist[current_psd_index][0] / 2;
//	nCells = (int)floor(PNM.settings.domain / cellSize);
//	cellSize = PNM.settings.domain / nCells;
//
//	//variable for storing each cell pore coords+diameter
//	vector<vector<vector<double>>> cellCoords(pow(nCells, 3));
//
//	int n_generated = 0;
//
//	//give seed to the random generation algorithm
//	srand(PNM.settings.seed);
//
//	//variable for storing randomly generated pore coords plus its diameter before pushing into the Pores array
//	vector<double> temp_pore = { ((double)rand() / RAND_MAX) * PNM.settings.domain, ((double)rand() / RAND_MAX) * PNM.settings.domain, ((double)rand() / RAND_MAX) * PNM.settings.domain };
//	temp_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//	pores.push_back(temp_pore);
//	PNM.network.diameter_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//
//
//	n_generated++;
//	if (n_generated >= nPores[current_psd_index])
//	{
//		current_psd_index++;
//		n_generated = 0;
//	}
//
//	//assign pore to a cell
//	xCell = (int)floor(temp_pore[0] / cellSize);
//	yCell = (int)floor(temp_pore[1] / cellSize);
//	zCell = (int)floor(temp_pore[2] / cellSize);
//	cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//	cellCoords[cellIndex].push_back(temp_pore);
//
//	//loop for generating the rest of the pores
//	int flag;
//	for (int i = 1; i < total_pn; i++)
//	{
//		//V1: assign pores to cells
//		flag = 1; //flag is 1 (overlap)to enter the loop, but reset to zero(no overlap) straight after and raised in the overlap loop
//		while (flag > 0)
//		{
//			flag = 0; //flag reset
//			//take a random coordinate and calculate its location in the cell grid
//			temp_pore.clear();
//			double failsafe = 0.999999999;
//			temp_pore = { ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe, ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe, ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe };
//			temp_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//			//pore HAS TO BE this diameter, assign now
//			PNM.network.diameter_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//
//			xCell = (int)floor(temp_pore[0] / cellSize);
//			yCell = (int)floor(temp_pore[1] / cellSize);
//			zCell = (int)floor(temp_pore[2] / cellSize);
//			//cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//			//cellCoords[cellIndex].push_back(temp_pore);
//
//			//run search in 25 cells around the cell of interest
//			int cellIndexSearch;
//			double distance2;
//			for (int iCell = max(0, xCell - 2); iCell < min(nCells, xCell + 3); iCell++)
//			{
//				for (int jCell = max(0, yCell - 2); jCell < min(nCells, yCell + 3); jCell++)
//				{
//					for (int kCell = max(0, zCell - 2); kCell < min(nCells, zCell + 3); kCell++)
//					{
//						cellIndexSearch = pow(nCells, 2)*iCell + nCells * jCell + kCell;
//						for (int pore_in_cell = 0; pore_in_cell < cellCoords[cellIndexSearch].size(); pore_in_cell++)
//						{
//							distance2 = pow((temp_pore[0] - cellCoords[cellIndexSearch][pore_in_cell][0]), 2) + pow((temp_pore[1] - cellCoords[cellIndexSearch][pore_in_cell][1]), 2) + pow((temp_pore[2] - cellCoords[cellIndexSearch][pore_in_cell][2]), 2);
//							if (distance2 <= pow(cellCoords[cellIndexSearch][pore_in_cell][3] / 2 + temp_pore[3] / 2, 2))
//							{
//								flag = 1;
//							}
//							if (flag > 0) { break; }
//						}
//						if (flag > 0) { break; }
//					}
//					if (flag > 0) { break; }
//				}
//				if (flag > 0) { break; }
//			}
//		}
//
//
//		cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//		cellCoords[cellIndex].push_back(temp_pore);
//		pores.push_back(temp_pore);
//
//		n_generated++;
//		if (n_generated >= nPores[current_psd_index])
//		{
//			current_psd_index++;
//			n_generated = 0;
//		}
//
//		if (i % 1000 == 0)
//		{
//			//cout << i << "\n";
//		}
//
//	}
//
//	cout << "Number of pores generated: " << pores.size() << endl;
//	PNM.network.pores = pores;
//
//	return PNM;
//}
//
//classPNM classPNM::poresGeneration2(classPNM PNM, classUtils utils)
//{
//	chrono::steady_clock::time_point begin = chrono::steady_clock::now();
//
//	PNM.poreSizeDist.clear();
//	PNM.network.pores.clear();
//	PNM.network.diameter_pore.clear();
//
//	//variable for holding the array of pore coordinates
//	vector<vector<double>> pores;
//
//	//converting string into double PSD
//	vector<double> temp;
//	for (int i = 0; i < PNM.poreSizeDist_raw.size(); i++)
//	{
//		temp.clear();
//		for (int j = 0; j < PNM.poreSizeDist_raw[i].size(); j++)
//		{
//			if (j == 0) { temp.push_back(stod(PNM.poreSizeDist_raw[i][j])*1e-9); }
//			else { temp.push_back(stod(PNM.poreSizeDist_raw[i][j])); }
//		}
//		PNM.poreSizeDist.push_back(temp);
//	}
//
//	//calculate number of pores per each poe size with given domain
//	vector<int> nPores = PNM.poresPsize(PNM);
//
//	//variable definition
//	int total_pn, current_psd_index, nCells, xCell, yCell, zCell, cellIndex, nCellsSmall, xCellSmall, yCellSmall, zCellSmall, cellIndexSmall;
//	double expected_porosity, cellSize, pores_per_cell, cellSizeSmall;
//	int nCellsMicro, xCellMicro, yCellMicro, zCellMicro, cellIndexMicro;
//	double cellSizeMicro;
//
//	total_pn = 0;
//	expected_porosity = 0;
//	//calculate total pore number
//	for (int i = 0; i < nPores.size(); i++)
//	{
//		if (!(nPores[i] < 1))
//		{
//			total_pn += nPores[i];
//			expected_porosity += PNM.poreSizeDist[i][1];
//		}
//	}
//	PNM.expectedP = expected_porosity;
//	cout << "Expected network porosity is: " << expected_porosity * 100.0 << " %\n";
//
//	//define which pore size we are starting generation from
//	current_psd_index = 0;
//	while (nPores[current_psd_index] < 1)
//	{
//		current_psd_index++;
//	}
//
//	//define how many pores we want to generate as a "large" grid:
//	double fraction = 0.0;
//	double cut_off = current_psd_index; //everything up to but NOT including THIS entry into PSD to be generated as "large" grid
//	double cut_off_micro = current_psd_index; ////everything up to but NOT including THIS entry into PSD to be generated as "small" grid
//	int large_pores = 0; //number of large pores to be generated
//	int small_pores = 0;
//	while (fraction < 0.001)
//	{
//		fraction += (double)nPores[cut_off] / total_pn;
//		large_pores += nPores[cut_off];
//		cut_off++;
//
//	};
//	cout << "Cut-off: " << cut_off << endl;
//
//	//divide donain into "large" cells, just larger than biggest pore radius
//	cellSize = PNM.poreSizeDist[current_psd_index][0] / 2;
//	nCells = (int)floor(PNM.settings.domain / cellSize);
//	cellSize = PNM.settings.domain / nCells;
//
//	//divide domain into "small" cells, based on the cut-off pore size
//	cellSizeSmall = PNM.poreSizeDist[cut_off][0] / 2;
//	nCellsSmall = (int)floor(PNM.settings.domain / cellSizeSmall);
//	cellSizeSmall = PNM.settings.domain / nCellsSmall;
//
//	//variable for storing each "large" cell pore coords+diameter
//	vector<vector<vector<double>>> cellCoords(pow(nCells, 3));
//	
//	//variable for storing each "small" cell pore coords+diameter
//	vector<vector<vector<double>>> cellCoordsSmall(pow(nCellsSmall, 3));
//
//
//	int n_generated = 0;
//
//	//give seed to the random generation algorithm
//	srand(PNM.settings.seed);
//	
//	//generate the first pore
//	//variable for storing randomly generated pore coords plus its diameter before pushing into the Pores array
//	vector<double> temp_pore = { ((double)rand() / RAND_MAX) * PNM.settings.domain, ((double)rand() / RAND_MAX) * PNM.settings.domain, ((double)rand() / RAND_MAX) * PNM.settings.domain };
//	temp_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//	pores.push_back(temp_pore);
//	PNM.network.diameter_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//
//
//	n_generated++;
//	if (n_generated >= nPores[current_psd_index])
//	{
//		current_psd_index++;
//		n_generated = 0;
//	}
//
//	//assign pore to a cell
//	xCell = (int)floor(temp_pore[0] / cellSize);
//	yCell = (int)floor(temp_pore[1] / cellSize);
//	zCell = (int)floor(temp_pore[2] / cellSize);
//	cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//	cellCoords[cellIndex].push_back(temp_pore);
//
//	cout << "Pore generation progress:\n";
//	cout << setprecision(2) << std::fixed;
//
//	//loop for generating the rest of the LARGE pores
//	int flag;
//	double progress;
//	for (int i = 1; i < large_pores; i++)
//	{
//		//V1: assign pores to cells
//		flag = 1; //flag is 1 (overlap)to enter the loop, but reset to zero(no overlap) straight after and raised in the overlap loop
//		while (flag > 0)
//		{
//			flag = 0; //flag reset
//			//take a random coordinate and calculate its location in the cell grid
//			temp_pore.clear();
//			double failsafe = 0.999999999;
//			temp_pore = { ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe, ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe, ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe };
//			temp_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//			//pore HAS TO BE this diameter, assign now
//			PNM.network.diameter_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//
//			xCell = (int)floor(temp_pore[0] / cellSize);
//			yCell = (int)floor(temp_pore[1] / cellSize);
//			zCell = (int)floor(temp_pore[2] / cellSize);
//			//cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//			//cellCoords[cellIndex].push_back(temp_pore);
//
//			//run search in 25 cells around the cell of interest
//			int cellIndexSearch;
//			double distance2;
//			for (int iCell = max(0, xCell - 2); iCell < min(nCells, xCell + 3); iCell++)
//			{
//				for (int jCell = max(0, yCell - 2); jCell < min(nCells, yCell + 3); jCell++)
//				{
//					for (int kCell = max(0, zCell - 2); kCell < min(nCells, zCell + 3); kCell++)
//					{
//						cellIndexSearch = pow(nCells, 2)*iCell + nCells * jCell + kCell;
//						for (int pore_in_cell = 0; pore_in_cell < cellCoords[cellIndexSearch].size(); pore_in_cell++)
//						{
//							distance2 = pow((temp_pore[0] - cellCoords[cellIndexSearch][pore_in_cell][0]), 2) + pow((temp_pore[1] - cellCoords[cellIndexSearch][pore_in_cell][1]), 2) + pow((temp_pore[2] - cellCoords[cellIndexSearch][pore_in_cell][2]), 2);
//							if (distance2 <= pow(cellCoords[cellIndexSearch][pore_in_cell][3] / 2 + temp_pore[3] / 2, 2))
//							{
//								flag = 1;
//							}
//							if (flag > 0) { break; }
//						}
//						if (flag > 0) { break; }
//					}
//					if (flag > 0) { break; }
//				}
//				if (flag > 0) { break; }
//			}
//		}
//
//
//		cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//		cellCoords[cellIndex].push_back(temp_pore);
//		pores.push_back(temp_pore);
//
//		n_generated++;
//		if (n_generated >= nPores[current_psd_index])
//		{
//			current_psd_index++;
//			n_generated = 0;
//		}
//
//		if (i % 10000 == 0)
//		{
//			progress = 100 * (i + 1) / (double)total_pn;
//			cout << "\r" << progress << "%" << flush;
//		}
//
//	}
//
//	//loop for generating SMALL pores
//	for (int i = large_pores; i < total_pn; i++)
//	{
//		//V1: assign pores to cells
//		flag = 1; //flag is 1 (overlap)to enter the loop, but reset to zero(no overlap) straight after and raised in the overlap loop
//		while (flag > 0)
//		{
//			flag = 0; //flag reset
//			//take a random coordinate and calculate its location in the cell grid
//			temp_pore.clear();
//			double failsafe = 0.999999999;
//			temp_pore = { ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe, ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe, ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe };
//			temp_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//			//pore HAS TO BE this diameter, assign now
//			PNM.network.diameter_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//
//			//TEST LARGE GRID FIRST
//			xCell = (int)floor(temp_pore[0] / cellSize);
//			yCell = (int)floor(temp_pore[1] / cellSize);
//			zCell = (int)floor(temp_pore[2] / cellSize);
//			//cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//			//cellCoords[cellIndex].push_back(temp_pore);
//
//			//run search in 25 cells around the cell of interest
//			int cellIndexSearch;
//			double distance2;
//			for (int iCell = max(0, xCell - 2); iCell < min(nCells, xCell + 3); iCell++)
//			{
//				for (int jCell = max(0, yCell - 2); jCell < min(nCells, yCell + 3); jCell++)
//				{
//					for (int kCell = max(0, zCell - 2); kCell < min(nCells, zCell + 3); kCell++)
//					{
//						cellIndexSearch = pow(nCells, 2)*iCell + nCells * jCell + kCell;
//						for (int pore_in_cell = 0; pore_in_cell < cellCoords[cellIndexSearch].size(); pore_in_cell++)
//						{
//							distance2 = pow((temp_pore[0] - cellCoords[cellIndexSearch][pore_in_cell][0]), 2) + pow((temp_pore[1] - cellCoords[cellIndexSearch][pore_in_cell][1]), 2) + pow((temp_pore[2] - cellCoords[cellIndexSearch][pore_in_cell][2]), 2);
//							if (distance2 <= pow(cellCoords[cellIndexSearch][pore_in_cell][3] / 2 + temp_pore[3] / 2, 2))
//							{
//								flag = 1;
//							}
//							if (flag > 0) { break; }
//						}
//						if (flag > 0) { break; }
//					}
//					if (flag > 0) { break; }
//				}
//				if (flag > 0) { break; }
//			}
//			//if fist check exits with flag > 0, it will exit instantly in the next loop, loss of efficiency is minimal
//
//			//TEST SMALL GRID
//			xCellSmall = (int)floor(temp_pore[0] / cellSizeSmall);
//			yCellSmall = (int)floor(temp_pore[1] / cellSizeSmall);
//			zCellSmall = (int)floor(temp_pore[2] / cellSizeSmall);
//			//cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//			//cellCoords[cellIndex].push_back(temp_pore);
//
//			//run search in 25 cells around the cell of interest
//			/*int cellIndexSearch;
//			double distance2;*/
//			for (int iCell = max(0, xCellSmall - 2); iCell < min(nCellsSmall, xCellSmall + 3); iCell++)
//			{
//				for (int jCell = max(0, yCellSmall - 2); jCell < min(nCellsSmall, yCellSmall + 3); jCell++)
//				{
//					for (int kCell = max(0, zCellSmall - 2); kCell < min(nCellsSmall, zCellSmall + 3); kCell++)
//					{
//						cellIndexSearch = pow(nCellsSmall, 2)*iCell + nCellsSmall * jCell + kCell;
//						for (int pore_in_cell = 0; pore_in_cell < cellCoordsSmall[cellIndexSearch].size(); pore_in_cell++)
//						{
//							distance2 = pow((temp_pore[0] - cellCoordsSmall[cellIndexSearch][pore_in_cell][0]), 2) + pow((temp_pore[1] - cellCoordsSmall[cellIndexSearch][pore_in_cell][1]), 2) + pow((temp_pore[2] - cellCoordsSmall[cellIndexSearch][pore_in_cell][2]), 2);
//							if (distance2 <= pow(cellCoordsSmall[cellIndexSearch][pore_in_cell][3] / 2 + temp_pore[3] / 2, 2))
//							{
//								flag = 1;
//							}
//							if (flag > 0) { break; }
//						}
//						if (flag > 0) { break; }
//					}
//					if (flag > 0) { break; }
//				}
//				if (flag > 0) { break; }
//			}
//		}
//
//		//push back generated pore into the small grid storage
//		cellIndexSmall = pow(nCellsSmall, 2)*xCellSmall + nCellsSmall * yCellSmall + zCellSmall;
//		cellCoordsSmall[cellIndexSmall].push_back(temp_pore);
//		pores.push_back(temp_pore);
//
//		n_generated++;
//		if (n_generated >= nPores[current_psd_index])
//		{
//			current_psd_index++;
//			n_generated = 0;
//		}
//
//		if (i % 10000 == 0)
//		{
//			progress = 100 * (i + 1) / (double)total_pn;
//			cout << "\r" << progress << "%" << flush;
//		}
//
//	}
//	progress = 100;
//	cout << "\r" << progress << "%" << flush;
//
//	cout << "\nNumber of pores generated: " << pores.size() << endl;
//	PNM.network.pores = pores;
//
//	return PNM;
//
//	chrono::steady_clock::time_point end = chrono::steady_clock::now();
//	string elapsedTime = "Elapsed time: " + to_string(chrono::duration_cast<std::chrono::seconds>(end - begin).count()) + " [s]";
//	cout << elapsedTime << endl;
//}
//
//classPNM classPNM::poresGeneration3(classPNM PNM, classUtils utils)
//{
//
//	PNM.poreSizeDist.clear();
//	PNM.network.pores.clear();
//	PNM.network.diameter_pore.clear();
//
//	//variable for holding the array of pore coordinates
//	vector<vector<double>> pores;
//
//	//converting string into double PSD
//	vector<double> temp;
//	for (int i = 0; i < PNM.poreSizeDist_raw.size(); i++)
//	{
//		temp.clear();
//		for (int j = 0; j < PNM.poreSizeDist_raw[i].size(); j++)
//		{
//			if (j == 0) { temp.push_back(stod(PNM.poreSizeDist_raw[i][j])*1e-9); }
//			else { temp.push_back(stod(PNM.poreSizeDist_raw[i][j])); }
//		}
//		PNM.poreSizeDist.push_back(temp);
//	}
//
//	//calculate number of pores per each poe size with given domain
//	vector<int> nPores = PNM.poresPsize(PNM);
//
//	//variable definition
//	int total_pn, current_psd_index, nCells, xCell, yCell, zCell, cellIndex;
//	int nCellsSmall, xCellSmall, yCellSmall, zCellSmall, cellIndexSmall;
//	double expected_porosity, cellSize, pores_per_cell;
//	double cellSizeSmall;
//	int nCellsMicro, xCellMicro, yCellMicro, zCellMicro, cellIndexMicro;
//	double cellSizeMicro;
//
//	total_pn = 0;
//	expected_porosity = 0;
//	//calculate total pore number
//	for (int i = 0; i < nPores.size(); i++)
//	{
//		if (!(nPores[i] < 1))
//		{
//			total_pn += nPores[i];
//			expected_porosity += PNM.poreSizeDist[i][1];
//		}
//	}
//	PNM.expectedP = expected_porosity;
//	cout << "Expected network porosity is: " << expected_porosity * 100.0 << " %\n";
//
//	//define which pore size we are starting generation from
//	current_psd_index = 0;
//	while (nPores[current_psd_index] < 1)
//	{
//		current_psd_index++;
//	}
//
//	//define how many pores we want to generate as a "large" grid:
//	double fraction = 0.0;
//	double cut_off = current_psd_index; //everything up to but NOT including THIS entry into PSD to be generated as "large" grid
//	double cut_off_micro = current_psd_index; ////everything up to but NOT including THIS entry into PSD to be generated as "small" grid
//	int large_pores = 0; //number of large pores to be generated
//	int small_pores = 0;
//	while (fraction < 0.0001)
//	{
//		fraction += (double)nPores[cut_off] / total_pn;
//		large_pores += nPores[cut_off];
//		cut_off++;
//
//	};
//	cout << "Cut-off: " << cut_off << endl;
//
//	fraction = 0.0;
//	while (fraction < 0.001)
//	{
//		fraction += (double)nPores[cut_off_micro] / total_pn;
//		small_pores += nPores[cut_off_micro];
//		cut_off_micro++;
//
//	};
//	cout << "Cut-off_micro: " << cut_off_micro << endl;
//
//	//divide donain into "large" cells, just larger than biggest pore radius
//	cellSize = PNM.poreSizeDist[current_psd_index][0] / 2;
//	nCells = (int)floor(PNM.settings.domain / cellSize);
//	cellSize = PNM.settings.domain / nCells;
//
//	//divide domain into "small" cells, based on the cut-off pore size
//	cellSizeSmall = PNM.poreSizeDist[cut_off][0] / 2;
//	nCellsSmall = (int)floor(PNM.settings.domain / cellSizeSmall);
//	cellSizeSmall = PNM.settings.domain / nCellsSmall;
//	
//	//divide domain into "micro" cells, based on the cut-off_micro pore size
//	cellSizeMicro = PNM.poreSizeDist[cut_off_micro][0] / 2;
//	nCellsMicro = (int)floor(PNM.settings.domain / cellSizeMicro);
//	cellSizeMicro = PNM.settings.domain / nCellsSmall;
//
//	//variable for storing each "large" cell pore coords+diameter
//	vector<vector<vector<double>>> cellCoords(pow(nCells, 3));
//
//	//variable for storing each "small" cell pore coords+diameter
//	vector<vector<vector<double>>> cellCoordsSmall(pow(nCellsSmall, 3));
//	
//	//variable for storing each "micro" cell pore coords+diameter
//	vector<vector<vector<double>>> cellCoordsMicro(pow(nCellsMicro, 3));
//
//
//	int n_generated = 0;
//
//	//give seed to the random generation algorithm
//	srand(PNM.settings.seed);
//
//	//generate the first pore
//	//variable for storing randomly generated pore coords plus its diameter before pushing into the Pores array
//	vector<double> temp_pore = { ((double)rand() / RAND_MAX) * PNM.settings.domain, ((double)rand() / RAND_MAX) * PNM.settings.domain, ((double)rand() / RAND_MAX) * PNM.settings.domain };
//	temp_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//	pores.push_back(temp_pore);
//	PNM.network.diameter_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//
//
//	n_generated++;
//	if (n_generated >= nPores[current_psd_index])
//	{
//		current_psd_index++;
//		n_generated = 0;
//	}
//
//	//assign pore to a cell
//	xCell = (int)floor(temp_pore[0] / cellSize);
//	yCell = (int)floor(temp_pore[1] / cellSize);
//	zCell = (int)floor(temp_pore[2] / cellSize);
//	cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//	cellCoords[cellIndex].push_back(temp_pore);
//
//	//loop for generating the rest of the LARGE pores
//	int flag;
//	for (int i = 1; i < large_pores; i++)
//	{
//		//V1: assign pores to cells
//		flag = 1; //flag is 1 (overlap)to enter the loop, but reset to zero(no overlap) straight after and raised in the overlap loop
//		while (flag > 0)
//		{
//			flag = 0; //flag reset
//			//take a random coordinate and calculate its location in the cell grid
//			temp_pore.clear();
//			double failsafe = 0.999999999;
//			temp_pore = { ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe, ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe, ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe };
//			temp_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//			//pore HAS TO BE this diameter, assign now
//			PNM.network.diameter_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//
//			xCell = (int)floor(temp_pore[0] / cellSize);
//			yCell = (int)floor(temp_pore[1] / cellSize);
//			zCell = (int)floor(temp_pore[2] / cellSize);
//			//cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//			//cellCoords[cellIndex].push_back(temp_pore);
//
//			//run search in 25 cells around the cell of interest
//			int cellIndexSearch;
//			double distance2;
//			for (int iCell = max(0, xCell - 2); iCell < min(nCells, xCell + 3); iCell++)
//			{
//				for (int jCell = max(0, yCell - 2); jCell < min(nCells, yCell + 3); jCell++)
//				{
//					for (int kCell = max(0, zCell - 2); kCell < min(nCells, zCell + 3); kCell++)
//					{
//						cellIndexSearch = pow(nCells, 2)*iCell + nCells * jCell + kCell;
//						for (int pore_in_cell = 0; pore_in_cell < cellCoords[cellIndexSearch].size(); pore_in_cell++)
//						{
//							distance2 = pow((temp_pore[0] - cellCoords[cellIndexSearch][pore_in_cell][0]), 2) + pow((temp_pore[1] - cellCoords[cellIndexSearch][pore_in_cell][1]), 2) + pow((temp_pore[2] - cellCoords[cellIndexSearch][pore_in_cell][2]), 2);
//							if (distance2 <= pow(cellCoords[cellIndexSearch][pore_in_cell][3] / 2 + temp_pore[3] / 2, 2))
//							{
//								flag = 1;
//							}
//							if (flag > 0) { break; }
//						}
//						if (flag > 0) { break; }
//					}
//					if (flag > 0) { break; }
//				}
//				if (flag > 0) { break; }
//			}
//		}
//
//
//		cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//		cellCoords[cellIndex].push_back(temp_pore);
//		pores.push_back(temp_pore);
//
//		n_generated++;
//		if (n_generated >= nPores[current_psd_index])
//		{
//			current_psd_index++;
//			n_generated = 0;
//		}
//
//		if (i % 1000 == 0)
//		{
//			//cout << i << "\n";
//		}
//
//	}
//
//	//loop for generating SMALL pores
//	for (int i = large_pores; i < small_pores; i++)
//	{
//		//V1: assign pores to cells
//		flag = 1; //flag is 1 (overlap)to enter the loop, but reset to zero(no overlap) straight after and raised in the overlap loop
//		while (flag > 0)
//		{
//			flag = 0; //flag reset
//			//take a random coordinate and calculate its location in the cell grid
//			temp_pore.clear();
//			double failsafe = 0.999999999;
//			temp_pore = { ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe, ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe, ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe };
//			temp_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//			//pore HAS TO BE this diameter, assign now
//			PNM.network.diameter_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//
//			//TEST LARGE GRID FIRST
//			xCell = (int)floor(temp_pore[0] / cellSize);
//			yCell = (int)floor(temp_pore[1] / cellSize);
//			zCell = (int)floor(temp_pore[2] / cellSize);
//			//cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//			//cellCoords[cellIndex].push_back(temp_pore);
//
//			//run search in 25 cells around the cell of interest
//			int cellIndexSearch;
//			double distance2;
//			for (int iCell = max(0, xCell - 2); iCell < min(nCells, xCell + 3); iCell++)
//			{
//				for (int jCell = max(0, yCell - 2); jCell < min(nCells, yCell + 3); jCell++)
//				{
//					for (int kCell = max(0, zCell - 2); kCell < min(nCells, zCell + 3); kCell++)
//					{
//						cellIndexSearch = pow(nCells, 2)*iCell + nCells * jCell + kCell;
//						for (int pore_in_cell = 0; pore_in_cell < cellCoords[cellIndexSearch].size(); pore_in_cell++)
//						{
//							distance2 = pow((temp_pore[0] - cellCoords[cellIndexSearch][pore_in_cell][0]), 2) + pow((temp_pore[1] - cellCoords[cellIndexSearch][pore_in_cell][1]), 2) + pow((temp_pore[2] - cellCoords[cellIndexSearch][pore_in_cell][2]), 2);
//							if (distance2 <= pow(cellCoords[cellIndexSearch][pore_in_cell][3] / 2 + temp_pore[3] / 2, 2))
//							{
//								flag = 1;
//							}
//							if (flag > 0) { break; }
//						}
//						if (flag > 0) { break; }
//					}
//					if (flag > 0) { break; }
//				}
//				if (flag > 0) { break; }
//			}
//			//if fist check exits with flag > 0, it will exit instantly in the next loop, loss of efficiency is minimal
//
//			//TEST SMALL GRID
//			xCellSmall = (int)floor(temp_pore[0] / cellSizeSmall);
//			yCellSmall = (int)floor(temp_pore[1] / cellSizeSmall);
//			zCellSmall = (int)floor(temp_pore[2] / cellSizeSmall);
//			//cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//			//cellCoords[cellIndex].push_back(temp_pore);
//
//			//run search in 25 cells around the cell of interest
//			/*int cellIndexSearch;
//			double distance2;*/
//			for (int iCell = max(0, xCellSmall - 2); iCell < min(nCellsSmall, xCellSmall + 3); iCell++)
//			{
//				for (int jCell = max(0, yCellSmall - 2); jCell < min(nCellsSmall, yCellSmall + 3); jCell++)
//				{
//					for (int kCell = max(0, zCellSmall - 2); kCell < min(nCellsSmall, zCellSmall + 3); kCell++)
//					{
//						cellIndexSearch = pow(nCellsSmall, 2)*iCell + nCellsSmall * jCell + kCell;
//						for (int pore_in_cell = 0; pore_in_cell < cellCoordsSmall[cellIndexSearch].size(); pore_in_cell++)
//						{
//							distance2 = pow((temp_pore[0] - cellCoordsSmall[cellIndexSearch][pore_in_cell][0]), 2) + pow((temp_pore[1] - cellCoordsSmall[cellIndexSearch][pore_in_cell][1]), 2) + pow((temp_pore[2] - cellCoordsSmall[cellIndexSearch][pore_in_cell][2]), 2);
//							if (distance2 <= pow(cellCoordsSmall[cellIndexSearch][pore_in_cell][3] / 2 + temp_pore[3] / 2, 2))
//							{
//								flag = 1;
//							}
//							if (flag > 0) { break; }
//						}
//						if (flag > 0) { break; }
//					}
//					if (flag > 0) { break; }
//				}
//				if (flag > 0) { break; }
//			}
//		}
//
//		//push back generated pore into the small grid storage
//		cellIndexSmall = pow(nCellsSmall, 2)*xCellSmall + nCellsSmall * yCellSmall + zCellSmall;
//		cellCoordsSmall[cellIndexSmall].push_back(temp_pore);
//		pores.push_back(temp_pore);
//
//		n_generated++;
//		if (n_generated >= nPores[current_psd_index])
//		{
//			current_psd_index++;
//			n_generated = 0;
//		}
//
//		if (i % 1000 == 0)
//		{
//			//cout << i << "\n";
//		}
//
//	}
//
//	//loop for generating MICRO pores
//	for (int i = small_pores; i < total_pn; i++)
//	{
//		//V1: assign pores to cells
//		flag = 1; //flag is 1 (overlap)to enter the loop, but reset to zero(no overlap) straight after and raised in the overlap loop
//		while (flag > 0)
//		{
//			flag = 0; //flag reset
//			//take a random coordinate and calculate its location in the cell grid
//			temp_pore.clear();
//			double failsafe = 0.999999999;
//			temp_pore = { ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe, ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe, ((double)rand() / RAND_MAX) * PNM.settings.domain * failsafe };
//			temp_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//			//pore HAS TO BE this diameter, assign now
//			PNM.network.diameter_pore.push_back(PNM.poreSizeDist[current_psd_index][0]);
//
//			//TEST LARGE GRID FIRST
//			xCell = (int)floor(temp_pore[0] / cellSize);
//			yCell = (int)floor(temp_pore[1] / cellSize);
//			zCell = (int)floor(temp_pore[2] / cellSize);
//			//cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//			//cellCoords[cellIndex].push_back(temp_pore);
//
//			//run search in 25 cells around the cell of interest
//			int cellIndexSearch;
//			double distance2;
//			for (int iCell = max(0, xCell - 2); iCell < min(nCells, xCell + 3); iCell++)
//			{
//				for (int jCell = max(0, yCell - 2); jCell < min(nCells, yCell + 3); jCell++)
//				{
//					for (int kCell = max(0, zCell - 2); kCell < min(nCells, zCell + 3); kCell++)
//					{
//						cellIndexSearch = pow(nCells, 2)*iCell + nCells * jCell + kCell;
//						for (int pore_in_cell = 0; pore_in_cell < cellCoords[cellIndexSearch].size(); pore_in_cell++)
//						{
//							distance2 = pow((temp_pore[0] - cellCoords[cellIndexSearch][pore_in_cell][0]), 2) + pow((temp_pore[1] - cellCoords[cellIndexSearch][pore_in_cell][1]), 2) + pow((temp_pore[2] - cellCoords[cellIndexSearch][pore_in_cell][2]), 2);
//							if (distance2 <= pow(cellCoords[cellIndexSearch][pore_in_cell][3] / 2 + temp_pore[3] / 2, 2))
//							{
//								flag = 1;
//							}
//							if (flag > 0) { break; }
//						}
//						if (flag > 0) { break; }
//					}
//					if (flag > 0) { break; }
//				}
//				if (flag > 0) { break; }
//			}
//			//if fist check exits with flag > 0, it will exit instantly in the next loop, loss of efficiency is minimal
//
//			//TEST SMALL GRID
//			xCellSmall = (int)floor(temp_pore[0] / cellSizeSmall);
//			yCellSmall = (int)floor(temp_pore[1] / cellSizeSmall);
//			zCellSmall = (int)floor(temp_pore[2] / cellSizeSmall);
//			//cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//			//cellCoords[cellIndex].push_back(temp_pore);
//
//			//run search in 25 cells around the cell of interest
//			/*int cellIndexSearch;
//			double distance2;*/
//			for (int iCell = max(0, xCellSmall - 2); iCell < min(nCellsSmall, xCellSmall + 3); iCell++)
//			{
//				for (int jCell = max(0, yCellSmall - 2); jCell < min(nCellsSmall, yCellSmall + 3); jCell++)
//				{
//					for (int kCell = max(0, zCellSmall - 2); kCell < min(nCellsSmall, zCellSmall + 3); kCell++)
//					{
//						cellIndexSearch = pow(nCellsSmall, 2)*iCell + nCellsSmall * jCell + kCell;
//						for (int pore_in_cell = 0; pore_in_cell < cellCoordsSmall[cellIndexSearch].size(); pore_in_cell++)
//						{
//							distance2 = pow((temp_pore[0] - cellCoordsSmall[cellIndexSearch][pore_in_cell][0]), 2) + pow((temp_pore[1] - cellCoordsSmall[cellIndexSearch][pore_in_cell][1]), 2) + pow((temp_pore[2] - cellCoordsSmall[cellIndexSearch][pore_in_cell][2]), 2);
//							if (distance2 <= pow(cellCoordsSmall[cellIndexSearch][pore_in_cell][3] / 2 + temp_pore[3] / 2, 2))
//							{
//								flag = 1;
//							}
//							if (flag > 0) { break; }
//						}
//						if (flag > 0) { break; }
//					}
//					if (flag > 0) { break; }
//				}
//				if (flag > 0) { break; }
//			}
//
//			//TEST MICRO GRID
//			xCellMicro = (int)floor(temp_pore[0] / cellSizeMicro);
//			yCellMicro = (int)floor(temp_pore[1] / cellSizeMicro);
//			zCellMicro = (int)floor(temp_pore[2] / cellSizeMicro);
//			//cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//			//cellCoords[cellIndex].push_back(temp_pore);
//
//			//run search in 25 cells around the cell of interest
//			/*int cellIndexSearch;
//			double distance2;*/
//			for (int iCell = max(0, xCellMicro - 2); iCell < min(nCellsMicro, xCellMicro + 3); iCell++)
//			{
//				for (int jCell = max(0, yCellMicro - 2); jCell < min(nCellsMicro, yCellMicro + 3); jCell++)
//				{
//					for (int kCell = max(0, zCellMicro - 2); kCell < min(nCellsMicro, zCellMicro + 3); kCell++)
//					{
//						cellIndexSearch = pow(nCellsMicro, 2)*iCell + nCellsMicro * jCell + kCell;
//						for (int pore_in_cell = 0; pore_in_cell < cellCoordsMicro[cellIndexSearch].size(); pore_in_cell++)
//						{
//							distance2 = pow((temp_pore[0] - cellCoordsMicro[cellIndexSearch][pore_in_cell][0]), 2) + pow((temp_pore[1] - cellCoordsMicro[cellIndexSearch][pore_in_cell][1]), 2) + pow((temp_pore[2] - cellCoordsMicro[cellIndexSearch][pore_in_cell][2]), 2);
//							if (distance2 <= pow(cellCoordsMicro[cellIndexSearch][pore_in_cell][3] / 2 + temp_pore[3] / 2, 2))
//							{
//								flag = 1;
//							}
//							if (flag > 0) { break; }
//						}
//						if (flag > 0) { break; }
//					}
//					if (flag > 0) { break; }
//				}
//				if (flag > 0) { break; }
//			}
//		}
//
//		//push back generated pore into the small grid storage
//		cellIndexMicro = pow(nCellsMicro, 2)*xCellMicro + nCellsMicro * yCellMicro + zCellMicro;
//		cellCoordsMicro[cellIndexMicro].push_back(temp_pore);
//		pores.push_back(temp_pore);
//
//		n_generated++;
//		if (n_generated >= nPores[current_psd_index])
//		{
//			current_psd_index++;
//			n_generated = 0;
//		}
//
//		if (i % 1000 == 0)
//		{
//			//cout << i << "\n";
//		}
//
//	}
//
//	cout << "Number of pores generated: " << pores.size() << endl;
//	PNM.network.pores = pores;
//
//	return PNM;
//}
//
//classPNM classPNM::throatsGeneration(classPNM PNM, classUtils utils)
//{
//	chrono::steady_clock::time_point begin = chrono::steady_clock::now();
//
//	PNM.network.throats.clear();
//	PNM.network.diameter_throat.clear();
//
//	//vector of throat connectivities {start pore index; end pore index; throat diameter}
//	vector<vector<double>> throats;
//
//	double running_z = 0;
//	int n_nearest = 10; //randomly connect to Z/2 out of nearest N pores, to provide some variability
//	srand(PNM.settings.seed);
//
//	//converting string into double PSD
//	vector<double> temp;
//	PNM.poreSizeDist.clear();
//	for (int i = 0; i < PNM.poreSizeDist_raw.size(); i++)
//	{
//		temp.clear();
//		for (int j = 0; j < PNM.poreSizeDist_raw[i].size(); j++)
//		{
//			if (j == 0) { temp.push_back(stod(PNM.poreSizeDist_raw[i][j])*1e-9); }
//			else { temp.push_back(stod(PNM.poreSizeDist_raw[i][j])); }
//		}
//		PNM.poreSizeDist.push_back(temp);
//	}
//
//	//calculate number of pores per each poe size with given domain
//	vector<int> nPores = PNM.poresPsize(PNM);
//
//	//variable definition
//	int total_pn, current_psd_index;
//	double expected_porosity;
//
//	total_pn = 0;
//	//calculate total pore number
//	for (int i = 0; i < nPores.size(); i++)
//	{
//		if (!(nPores[i] < 1))
//		{
//			total_pn += nPores[i];
//		}
//	}
//
//	//define which pore size we are starting generation from
//	current_psd_index = 0;
//	while (nPores[current_psd_index] < 1)
//	{
//		current_psd_index++;
//	}
//
//	//define how many pores we want to generate as a "large" grid:
//	double fraction = 0.0;
//	int cut_off = current_psd_index; //everything up to but NOT including THIS entry into PSD to be generated as "large" grid
//	int large_pores = 0; //number of large pores to be generated
//	while (fraction < 0.001)
//	{
//		fraction += (double)nPores[cut_off] / total_pn;
//		large_pores += nPores[cut_off];
//		cut_off++;
//
//	};
//	cout << "Cut-off: " << cut_off << endl;
//
//
//	//assign large pores into cells:
//	//divide donain into cells, just larger than biggest pore radius
//	double cellSize = PNM.network.pores[0][3] / 2;
//	int nCells = (int)floor(PNM.settings.domain / cellSize);
//	cellSize = PNM.settings.domain / nCells;
//
//	vector<vector<int>> cellCoords(pow(nCells, 3));
//	int xCell, yCell, zCell, cellIndex;
//	for (int i = 0; i < large_pores; i++)
//	{
//		xCell = (int)floor(PNM.network.pores[i][0] / cellSize);
//		yCell = (int)floor(PNM.network.pores[i][1] / cellSize);
//		zCell = (int)floor(PNM.network.pores[i][2] / cellSize);
//		cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//		cellCoords[cellIndex].push_back(i);
//	}
//	int assignment_check = 0; //check if all pors have been assigned to cells correctly
//	for (int i = 0; i < cellCoords.size(); i++)
//	{
//		assignment_check += cellCoords[i].size();
//	}
//
//	//assign SMALL pores into cells:
//	//divide donain into cells:
//	//divide domain into "small" cells, based on the cut-off pore size
//	double cellSizeSmall = PNM.poreSizeDist[cut_off][0] / 2;
//	int nCellsSmall = (int)floor(PNM.settings.domain / cellSizeSmall);
//	cellSizeSmall = PNM.settings.domain / nCellsSmall;
//
//	vector<vector<int>> cellCoordsSmall(pow(nCellsSmall, 3));
//
//	int xCellSmall, yCellSmall, zCellSmall, cellIndexSmall;
//	for (int i = large_pores; i < total_pn; i++)
//	{
//		xCellSmall = (int)floor(PNM.network.pores[i][0] / cellSizeSmall);
//		yCellSmall = (int)floor(PNM.network.pores[i][1] / cellSizeSmall);
//		zCellSmall = (int)floor(PNM.network.pores[i][2] / cellSizeSmall);
//		cellIndexSmall = pow(nCellsSmall, 2)*xCellSmall + nCellsSmall * yCellSmall + zCellSmall;
//		cellCoordsSmall[cellIndexSmall].push_back(i);
//	}
//	for (int i = 0; i < cellCoordsSmall.size(); i++)
//	{
//		assignment_check += cellCoordsSmall[i].size();
//	}
//	cout << assignment_check << endl;
//
//	//vector for storing existing connections for checking at generation
//	vector<vector<int>> existing_conns(total_pn);
//	//vector<int> temp_throats; //temporary storage of connectivity to be pushed back at the end of the loop
//
//	vector <vector<double>> distance_array; //array for storing sorted distance array for closest pores, {distance; pore index}
//
//	cout << "Throat generation progress:\n";
//	cout << setprecision(2) << std::fixed;
//	double progress;
//	
//
////main loop for creacting Z/2 throats for each pore
//	for (int i = 0; i < total_pn; i++)
//	{
//		distance_array.clear();
//		int large_grid_search_radius = 2;
//		int small_grid_search_radius = 4;
//		//search for closest pores within nearest 5 cells of the LARGE grid
//		
//		xCell = (int)floor(PNM.network.pores[i][0] / cellSize);
//		yCell = (int)floor(PNM.network.pores[i][1] / cellSize);
//		zCell = (int)floor(PNM.network.pores[i][2] / cellSize);
//		cellIndex = pow(nCells, 2)*xCell + nCells * yCell + zCell;
//		int cellIndexSearch;
//		double distance2;
//		for (int iCell = max(0, xCell - large_grid_search_radius); iCell < min(nCells, xCell + large_grid_search_radius + 1); iCell++)
//		{
//			for (int jCell = max(0, yCell - large_grid_search_radius); jCell < min(nCells, yCell + large_grid_search_radius + 1); jCell++)
//			{
//				for (int kCell = max(0, zCell - large_grid_search_radius); kCell < min(nCells, zCell + large_grid_search_radius + 1); kCell++)
//				{
//					cellIndexSearch = pow(nCells, 2)*iCell + nCells * jCell + kCell;
//					for (int pore_in_cell = 0; pore_in_cell < cellCoords[cellIndexSearch].size(); pore_in_cell++)
//					{
//						//construct the distance array for the pores within the searched cells 
//						distance_array.push_back({ (pow((PNM.network.pores[i][0] - PNM.network.pores[cellCoords[cellIndexSearch][pore_in_cell]][0]), 2) + pow((PNM.network.pores[i][1] - PNM.network.pores[cellCoords[cellIndexSearch][pore_in_cell]][1]), 2) + pow((PNM.network.pores[i][2] - PNM.network.pores[cellCoords[cellIndexSearch][pore_in_cell]][2]), 2)), (double)cellCoords[cellIndexSearch][pore_in_cell] });
//					}
//				}
//			}
//		}
//
//		//search for closest pores within nearest x cells of the SMALL grid
//		xCellSmall = (int)floor(PNM.network.pores[i][0] / cellSizeSmall);
//		yCellSmall = (int)floor(PNM.network.pores[i][1] / cellSizeSmall);
//		zCellSmall = (int)floor(PNM.network.pores[i][2] / cellSizeSmall);
//		cellIndexSmall = pow(nCellsSmall, 2)*xCellSmall + nCellsSmall * yCellSmall + zCellSmall;
//		//int cellIndexSearchSmall;
//		//double distance2;
//		for (int iCell = max(0, xCellSmall - small_grid_search_radius); iCell < min(nCellsSmall, xCellSmall + small_grid_search_radius + 1); iCell++)
//		{
//			for (int jCell = max(0, yCellSmall - small_grid_search_radius); jCell < min(nCellsSmall, yCellSmall + small_grid_search_radius + 1); jCell++)
//			{ 
//				for (int kCell = max(0, zCellSmall - small_grid_search_radius); kCell < min(nCellsSmall, zCellSmall + small_grid_search_radius + 1); kCell++)
//				{
//					cellIndexSearch = pow(nCellsSmall, 2)*iCell + nCellsSmall * jCell + kCell;
//					for (int pore_in_cell = 0; pore_in_cell < cellCoordsSmall[cellIndexSearch].size(); pore_in_cell++)
//					{
//						//construct the distance array for the pores within the searched cells 
//						distance_array.push_back({ (pow((PNM.network.pores[i][0] - PNM.network.pores[cellCoordsSmall[cellIndexSearch][pore_in_cell]][0]), 2) + pow((PNM.network.pores[i][1] - PNM.network.pores[cellCoordsSmall[cellIndexSearch][pore_in_cell]][1]), 2) + pow((PNM.network.pores[i][2] - PNM.network.pores[cellCoordsSmall[cellIndexSearch][pore_in_cell]][2]), 2)), (double)cellCoordsSmall[cellIndexSearch][pore_in_cell] });
//					}
//				}
//			}
//		}
//		//sorting the distance array 
//		if (distance_array.size() < n_nearest)
//		{
//			cout << "Not enough neighbours for pore " << i << ", cannot proceed generation!" << endl;
//			exit(1);
//		}
//
//		sort(distance_array.begin(), distance_array.end(), sortColumn0);
//
//		int connection, flag; //which closest pore do we try to connect to?
//		if (running_z < PNM.settings.coordination)
//		{
//			for (int j = 0; j < ceil(PNM.settings.coordination / 2); j++)
//			{
//				//check for existing connection on the pore of interest
//				flag = 0;
//				connection = rand() % n_nearest + 1; //desired pore to connect to from the list of closest pores
//
//				while (flag >= 0)
//				{
//					if (existing_conns[i].empty()) { flag = -1; } //if there are no existing connections on pore i, skip
//					else
//					{
//						for (int k = 0; k < existing_conns[i].size(); k++)
//						{
//							if (existing_conns[i][k] == distance_array[connection][1])
//							{
//								connection = rand() % n_nearest + 1;
//								flag = 0;
//								break;
//							}
//							flag = -1;
//						}
//					}
//				}
//
//
//				//update the connection array adn the diameter array
//				if (i < distance_array[connection][1])
//				{
//					throats.push_back({ (double)i, distance_array[connection][1], PNM.network.pores[i][3] });
//				}
//				else
//				{
//					throats.push_back({ distance_array[connection][1], (double)i, PNM.network.pores[i][3] });
//				}
//				PNM.network.diameter_throat.push_back(PNM.network.pores[i][3]);
//
//				//update the existing connection array with the new connection
//				existing_conns[i].push_back((int)distance_array[connection][1]);
//				existing_conns[(int)distance_array[connection][1]].push_back(i);
//			}
//		}
//		else
//		{
//			for (int j = 0; j < floor(PNM.settings.coordination / 2); j++)
//			{
//				//check for existing connection on the pore of interest
//				flag = 0;
//				connection = rand() % n_nearest + 1; //desired pore to connect to from the list of closest pores
//
//				while (flag >= 0)
//				{
//					if (existing_conns[i].empty()) { flag = -1; } //if there are no existing connections on pore i, skip
//					else
//					{
//						for (int k = 0; k < existing_conns[i].size(); k++)
//						{
//							if (existing_conns[i][k] == distance_array[connection][1])
//							{
//								connection = rand() % n_nearest + 1;
//								flag = 0;
//								break;
//							}
//							flag = -1;
//						}
//					}
//				}
//
//
//				//update the connection array adn the diameter array
//				if (i < distance_array[connection][1])
//				{
//					throats.push_back({ (double)i, distance_array[connection][1], PNM.network.pores[i][3] });
//				}
//				else
//				{
//					throats.push_back({ distance_array[connection][1], (double)i, PNM.network.pores[i][3] });
//				}
//				PNM.network.diameter_throat.push_back(PNM.network.pores[i][3]);
//
//				//update the existing connection array with the new connection
//				existing_conns[i].push_back((int)distance_array[connection][1]);
//				existing_conns[(int)distance_array[connection][1]].push_back(i);
//			}
//		}
//		running_z = 2 * (double)throats.size() / (double)(i + 1);
//
//		if (i % 1000 == 0)
//		{
//			progress = 100 * (i + 1) / (double)total_pn;
//			cout << "\r" << progress << "%" << flush;
//		}
//	}
//	progress = 100;
//	cout << "\r" << progress << "%" << flush;
//	cout << "\nNumber of throats generated: " << throats.size() << "\n";
//	PNM.network.throats = throats;
//
//	return PNM;
//
//	chrono::steady_clock::time_point end = chrono::steady_clock::now();
//	string elapsedTime = "Elapsed time: " + to_string(chrono::duration_cast<std::chrono::seconds>(end - begin).count()) + " [s]";
//	cout << elapsedTime << endl;
//}


classPNM classPNM::label_faces(classPNM PNM, classUtils utils)
{
	for (int i = 0; i < PNM.network.pores.size(); i++)
	{
		//X-axis: left to right
		if (PNM.network.pores[i][0] <= PNM.settings.tol*PNM.settings.domain) { PNM.network.left.push_back(true); }
		else { PNM.network.left.push_back(false); }
		if (PNM.network.pores[i][0] >= (1 - PNM.settings.tol)*PNM.settings.domain) { PNM.network.right.push_back(true); }
		else { PNM.network.right.push_back(false); }

		//Y-axis: back to front
		if (PNM.network.pores[i][1] <= PNM.settings.tol*PNM.settings.domain) { PNM.network.back.push_back(true); }
		else { PNM.network.back.push_back(false); }
		if (PNM.network.pores[i][1] >= (1 - PNM.settings.tol)*PNM.settings.domain) { PNM.network.front.push_back(true); }
		else { PNM.network.front.push_back(false); }

		//Z-axis: bottom to top
		if (PNM.network.pores[i][2] <= PNM.settings.tol*PNM.settings.domain) { PNM.network.bottom.push_back(true); }
		else { PNM.network.bottom.push_back(false); }
		if (PNM.network.pores[i][2] >= (1 - PNM.settings.tol)*PNM.settings.domain) { PNM.network.top.push_back(true); }
		else { PNM.network.top.push_back(false); }
	}
	utils.writeLine("Faces Labelled");

	return PNM;
}

vector<int> classPNM::poresPsize(classPNM PNM)
{
	double pi = 3.1415926;
	vector<int> result;
	double volume_per_element;
	double throatLengthFactor = 2.1; //ratio of throat length to pore diameter 

	for (int i = 0; i < PNM.poreSizeDist.size(); i++)
	{
		volume_per_element = 4 / 3 * pi*pow(PNM.poreSizeDist[i][0] / 2, 3) + throatLengthFactor * PNM.settings.coordination / 2 * pi*pow(PNM.poreSizeDist[i][0] / 2, 2)*(PNM.poreSizeDist[i][0]); 
		result.push_back(round(PNM.poreSizeDist[i][1] * pow(PNM.settings.domain, 3) / volume_per_element));

	}
	return result;

}

double classPNM::actualPorosity(classPNM PNM)
{
	double actualP;
	double totalV = pow(PNM.settings.domain, 3);
	double poreV = 0;
	double pi = 3.1415926;
	int Np = PNM.network.pores.size();
	int Nt = PNM.network.throats.size();
	//calculate volume for each element
	for (int i = 0; i < Np; i++)
	{
		poreV += 4 / 3 * pi * pow(PNM.network.diameter_pore[i] / 2, 3);
	}
	for (int i = 0; i < Nt; i++)
	{
		poreV += pi * pow(PNM.network.diameter_throat[i] / 2, 2) * PNM.network.conduit_length_throat[i][1];
	}
	actualP = poreV / totalV;

	return actualP;
}

classPNM classPNM::propertyCorrection(classPNM PNM, classUtils utils)
{
	double actualP;
	PNM.methods.methodsNetwork.generate_network_properties(PNM.network, PNM.methods.methodsNetwork);
	actualP = PNM.actualPorosity(PNM);
	cout << "Actual network porosity is: " << actualP*100.0 << " %\n";
	return PNM;
}


classPNM classPNM::csvReadNetwork_python(string loadPath, classPNM PNM, classUtils utils)
{
	PNM.network.pores.clear();
	PNM.network.throats.clear();
	PNM.network.diameter_pore.clear();
	PNM.network.diameter_throat.clear();
	string path = loadPath + ".csv";
	//string loadPathThroats = loadPath + "_throats.csv";

	//std::vector<std::vector<std::string> > pores_raw = utils.csvRead(loadPathPores);
	std::vector<std::vector<std::string> > raw = utils.csvRead(path);

	//converting string into double PSD
	vector<double> temp;
	for (int i = 1; i < raw.size(); i++)
	{
		temp.clear();
		if (!raw[i][10].empty())
		{
			temp.push_back(stod(raw[i][8]));
			temp.push_back(stod(raw[i][9]));
			temp.push_back(stod(raw[i][10]));
			temp.push_back(stod(raw[i][7]));

			PNM.network.pores.push_back(temp);
			PNM.network.diameter_pore.push_back(stod(raw[i][7]));
		}
		temp.clear();

		temp.push_back(stod(raw[i][3]));
		temp.push_back(stod(raw[i][4]));
		temp.push_back(stod(raw[i][2]));
		PNM.network.throats.push_back(temp);
		PNM.network.diameter_throat.push_back(stod(raw[i][2]));
	}
	

	return PNM;
}


classPNM classPNM::initial_conditions(classPNM PNM, classUtils utils)
{
	//calculate for pores
	for (int i = 0; i < PNM.network.pores.size(); i++)
	{
		PNM.phase.pressure_pore.push_back(PNM.settings.Pinit);
		PNM.phase.temperature_pore.push_back(PNM.settings.Tinit);
	}
	//calculate for throats
	for (int i = 0; i < PNM.network.throats.size(); i++)
	{
		PNM.phase.pressure_throat.push_back(PNM.settings.Pinit);
		PNM.phase.temperature_throat.push_back(PNM.settings.Tinit);
	}
	return PNM;
}

classPNM classPNM::test(classPNM PNM, classUtils utils)
{
	int Nt = PNM.network.throats.size();
	PNM.physics.hydraulic_conductance_throat = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200 };
	PNM.network.left = { true, true, true , true , false, false , false , false , false , false , false , false };
	PNM.network.right = { false, false , false , false , false , false , false , false, true, true, true , true };
	return PNM;
}

//

