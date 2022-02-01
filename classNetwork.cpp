#include "classPNM.h"

//
//NETWORK
//

//classPNM classPNM::classNetwork::cross_sectional_area_simple(classPNM PNM, classUtils utils)
//{
//	PNM.network.area_pore.clear();
//	PNM.network.area_troat.clear();
//	//calculate for pores
//	for (int i = 0; i < PNM.network.pores.size(); i++)
//	{
//		PNM.network.area_pore.push_back(3.1415926 / 4 * pow(PNM.network.diameter_pore[i], 2));
//	}
//	//calculate for throats
//	for (int i = 0; i < PNM.network.throats.size(); i++)
//	{
//		PNM.network.area_troat.push_back(3.1415926 / 4 * pow(PNM.network.diameter_throat[i], 2));
//	}
//
//	return PNM;
//}
//
//classPNM classPNM::classNetwork::endpoints_spherical(classPNM PNM, classUtils utils)
//{
//	//Calculate the coordinates of throat endpoints, assuming spherical pores.
//		//This model accounts for the overlapping lens between pores and throats.
//
//	PNM.network.endpoints_throat.clear();
//
//	double L, L1, L2, D1, D2, Dt;
//	int p1, p2; //pore indexes for each throat
//	vector<double> temp_coord, unit_vec_P1T, unit_vec_P2T;
//	vector<vector<double>> temp_endpoints;
//
//	//calculate for throats
//	for (int i = 0; i < PNM.network.throats.size(); i++)
//	{
//		temp_coord.clear();
//		temp_endpoints.clear();
//		unit_vec_P1T.clear();
//		unit_vec_P2T.clear();
//
//		p1 = PNM.network.throats[i][0];
//		p2 = PNM.network.throats[i][1];
//		Dt = PNM.network.diameter_throat[i];
//		D1 = PNM.network.diameter_pore[p1];
//		D2 = PNM.network.diameter_pore[p2];
//
//		L = sqrt(pow(PNM.network.pores[p1][0] - PNM.network.pores[p2][0], 2) + pow(PNM.network.pores[p1][1] - PNM.network.pores[p2][1], 2) + pow(PNM.network.pores[p1][2] - PNM.network.pores[p2][2], 2)) + 1e-15;
//
//		//L1
//		if (Dt > D1) {
//			L1 = 0.5*D1;
//		}
//		else {
//			L1 = sqrt(pow(D1, 2) - pow(Dt, 2)) / 2;
//		}
//		//L2
//		if (Dt > D2) {
//			L2 = 0.5*D2;
//		}
//		else {
//			L2 = sqrt(pow(D2, 2) - pow(Dt, 2)) / 2;
//		}
//		//unit vectors for each pore
//		unit_vec_P1T.push_back((PNM.network.pores[p2][0] - PNM.network.pores[p1][0]) / L);
//		unit_vec_P1T.push_back((PNM.network.pores[p2][1] - PNM.network.pores[p1][1]) / L);
//		unit_vec_P1T.push_back((PNM.network.pores[p2][2] - PNM.network.pores[p1][2]) / L);
//
//		unit_vec_P2T.push_back(-(PNM.network.pores[p2][0] - PNM.network.pores[p1][0]) / L);
//		unit_vec_P2T.push_back(-(PNM.network.pores[p2][1] - PNM.network.pores[p1][1]) / L);
//		unit_vec_P2T.push_back(-(PNM.network.pores[p2][2] - PNM.network.pores[p1][2]) / L);
//
//		//create each endpoint and push into temp vector
//		temp_coord.push_back(PNM.network.pores[p1][0] + L1 * unit_vec_P1T[0]);
//		temp_coord.push_back(PNM.network.pores[p1][1] + L1 * unit_vec_P1T[1]);
//		temp_coord.push_back(PNM.network.pores[p1][2] + L1 * unit_vec_P1T[2]);
//		temp_endpoints.push_back(temp_coord);
//		temp_coord.clear();
//		temp_coord.push_back(PNM.network.pores[p2][0] + L2 * unit_vec_P2T[0]);
//		temp_coord.push_back(PNM.network.pores[p2][1] + L2 * unit_vec_P2T[1]);
//		temp_coord.push_back(PNM.network.pores[p2][2] + L2 * unit_vec_P2T[2]);
//		temp_endpoints.push_back(temp_coord);
//
//		PNM.network.endpoints_throat.push_back(temp_endpoints);
//
//	}
//
//	return PNM;
//}
//
//classPNM classPNM::classNetwork::length_throat_piecewise(classPNM PNM, classUtils utils)
//{
//	PNM.network.length_throat.clear();
//
//	vector<double> EP1, EP2; //temp pore cooredinates and throat endpoints
//	double Lt;
//
//	//calculate for throats
//	for (int i = 0; i < PNM.network.throats.size(); i++)
//	{
//		EP1.clear(), EP2.clear();
//		EP1 = PNM.network.endpoints_throat[i][0];
//		EP2 = PNM.network.endpoints_throat[i][1];
//
//		Lt = sqrt(pow(EP1[0] - EP2[0], 2) + pow(EP1[1] - EP2[1], 2) + pow(EP1[2] - EP2[2], 2));
//		PNM.network.length_throat.push_back(Lt);
//	}
//	return PNM;
//}
//
//classPNM classPNM::classNetwork::conduit_length_throat_def(classPNM PNM, classUtils utils)
//{
//	PNM.network.conduit_length_throat.clear();
//	vector<double> C1, C2, EP1, EP2; //temp pore cooredinates and throat endpoints
//	int p1, p2; //pore indexes for each throat
//	double Lt, L1, L2;
//	vector<double> temp; //temp storage for conduit lengths {pore1, hroat, pore2}
//
//	//calculate for throats
//	for (int i = 0; i < PNM.network.throats.size(); i++)
//	{
//		C1.clear(); C2.clear(); EP1.clear(), EP2.clear();
//		temp.clear();
//		p1 = PNM.network.throats[i][0];
//		p2 = PNM.network.throats[i][1];
//
//		C1 = PNM.network.pores[p1];
//		C2 = PNM.network.pores[p2];
//		EP1 = PNM.network.endpoints_throat[i][0];
//		EP2 = PNM.network.endpoints_throat[i][1];
//		Lt = PNM.network.length_throat[i];
//
//		L1 = sqrt(pow(C1[0] - EP1[0], 2) + pow(C1[1] - EP1[1], 2) + pow(C1[2] - EP1[2], 2));
//		L2 = sqrt(pow(C2[0] - EP2[0], 2) + pow(C2[1] - EP2[1], 2) + pow(C2[2] - EP2[2], 2));
//
//		Lt = max(Lt, 1e-15);
//
//		temp = { L1, Lt, L2 };
//		PNM.network.conduit_length_throat.push_back(temp);
//	}
//	return PNM;
//}
//
//classPNM classPNM::classNetwork::shape_factor_cone_stick(classPNM PNM, classUtils utils)
//{
//	PNM.network.shape_factor_throat.clear();
//	double A1, At, A2, d1, dt, d2, L1, Lt, L2; //pore and throat area, diameter and conduit length
//	int p1, p2; //pore indexes for each throat
//	double SFt, SF1, SF2; //shape folw factors for pores and throat
//	vector<double> temp; //temp storage for shape factors {pore1, hroat, pore2}
//
//	double F1, Ft, F2; //volume
//	const double pi = 3.1415926;
//
//	//calculate for throats
//	for (int i = 0; i < PNM.network.throats.size(); i++)
//	{
//		temp.clear();
//		p1 = PNM.network.throats[i][0];
//		p2 = PNM.network.throats[i][1];
//
//		//Getting equivalent areas
//		A1 = PNM.network.area_pore[p1];
//		At = PNM.network.area_troat[i];
//		A2 = PNM.network.area_pore[p2];
//		//getting diameter
//		d1 = PNM.network.diameter_pore[p1];
//		dt = PNM.network.diameter_throat[i];
//		d2 = PNM.network.diameter_pore[p2];
//		//Getting conduit lengths
//		L1 = PNM.network.conduit_length_throat[i][0];
//		Lt = PNM.network.conduit_length_throat[i][1];
//		L2 = PNM.network.conduit_length_throat[i][2];
//
//		if (L1 <= 0) { SF1 = 1; }
//		else
//		{
//			F1 = 16.0 / 3.0 * (L1 * (pow(d1, 2) + d1 * dt + pow(dt, 2)) / (pow(d1, 3)*pow(dt, 3)*pow(pi, 2)));
//			SF1 = L1 / (pow(A1, 2)*F1);
//		}
//		if (Lt <= 0) { SFt = 1; }
//		else
//		{
//			Ft = Lt / pow(At, 2);
//			SFt = Lt / (pow(At, 2)*Ft);
//		}
//		if (L2 <= 0) { SF2 = 1; }
//		else
//		{
//			F2 = 16.0 / 3.0 * (L2 * (pow(d2, 2) + d2 * dt + pow(dt, 2)) / (pow(d2, 3)*pow(dt, 3)*pow(pi, 2)));
//			SF2 = L2 / (pow(A2, 2)*F2);
//		}
//
//		temp = { SF1, SFt, SF2 };
//		PNM.network.shape_factor_throat.push_back(temp);
//	}
//
//
//	return PNM;
//}
//
//classPNM classPNM::classNetwork::generate_network_properties(classPNM PNM, classUtils utils)
//{
//	PNM = PNM.network.cross_sectional_area_simple(PNM, utils);
//	PNM = PNM.network.endpoints_spherical(PNM, utils);
//	PNM = PNM.network.length_throat_piecewise(PNM, utils);
//	PNM = PNM.network.conduit_length_throat_def(PNM, utils);
//	PNM = PNM.network.shape_factor_cone_stick(PNM, utils);
//	utils.writeLine("Network Properties calculated");
//	return PNM;
//}

//

void classPNM::classMethods::classMethodsNetwork::cross_sectional_area_simple(classPNM::classNetwork& network)
{
	network.area_pore.clear();
	network.area_troat.clear();
	//calculate for pores
	for (int i = 0; i < network.pores.size(); i++)
	{
		network.area_pore.push_back(3.1415926 / 4 * pow(network.diameter_pore[i], 2));
	}
	//calculate for throats
	for (int i = 0; i < network.throats.size(); i++)
	{
		network.area_troat.push_back(3.1415926 / 4 * pow(network.diameter_throat[i], 2));
	}

}

void classPNM::classMethods::classMethodsNetwork::endpoints_spherical(classPNM::classNetwork& network)
{
	//Calculate the coordinates of throat endpoints, assuming spherical pores.
		//This model accounts for the overlapping lens between pores and throats.

	network.endpoints_throat.clear();

	double L, L1, L2, D1, D2, Dt;
	int p1, p2; //pore indexes for each throat
	vector<double> temp_coord, unit_vec_P1T, unit_vec_P2T;
	vector<vector<double>> temp_endpoints;

	//calculate for throats
	for (int i = 0; i < network.throats.size(); i++)
	{
		temp_coord.clear();
		temp_endpoints.clear();
		unit_vec_P1T.clear();
		unit_vec_P2T.clear();

		p1 = network.throats[i][0];
		p2 = network.throats[i][1];
		Dt = network.diameter_throat[i];
		D1 = network.diameter_pore[p1];
		D2 = network.diameter_pore[p2];

		L = sqrt(pow(network.pores[p1][0] - network.pores[p2][0], 2) + pow(network.pores[p1][1] - network.pores[p2][1], 2) + pow(network.pores[p1][2] - network.pores[p2][2], 2)) + 1e-15;

		//L1
		if (Dt > D1) {
			L1 = 0.5*D1;
		}
		else {
			L1 = sqrt(pow(D1, 2) - pow(Dt, 2)) / 2;
		}
		//L2
		if (Dt > D2) {
			L2 = 0.5*D2;
		}
		else {
			L2 = sqrt(pow(D2, 2) - pow(Dt, 2)) / 2;
		}
		//unit vectors for each pore
		unit_vec_P1T.push_back((network.pores[p2][0] - network.pores[p1][0]) / L);
		unit_vec_P1T.push_back((network.pores[p2][1] - network.pores[p1][1]) / L);
		unit_vec_P1T.push_back((network.pores[p2][2] - network.pores[p1][2]) / L);

		unit_vec_P2T.push_back(-(network.pores[p2][0] - network.pores[p1][0]) / L);
		unit_vec_P2T.push_back(-(network.pores[p2][1] - network.pores[p1][1]) / L);
		unit_vec_P2T.push_back(-(network.pores[p2][2] - network.pores[p1][2]) / L);

		//create each endpoint and push into temp vector
		temp_coord.push_back(network.pores[p1][0] + L1 * unit_vec_P1T[0]);
		temp_coord.push_back(network.pores[p1][1] + L1 * unit_vec_P1T[1]);
		temp_coord.push_back(network.pores[p1][2] + L1 * unit_vec_P1T[2]);
		temp_endpoints.push_back(temp_coord);
		temp_coord.clear();
		temp_coord.push_back(network.pores[p2][0] + L2 * unit_vec_P2T[0]);
		temp_coord.push_back(network.pores[p2][1] + L2 * unit_vec_P2T[1]);
		temp_coord.push_back(network.pores[p2][2] + L2 * unit_vec_P2T[2]);
		temp_endpoints.push_back(temp_coord);

		network.endpoints_throat.push_back(temp_endpoints);

	}
}

void classPNM::classMethods::classMethodsNetwork::length_throat_piecewise(classPNM::classNetwork& network)
{
	network.length_throat.clear();

	vector<double> EP1, EP2; //temp pore cooredinates and throat endpoints
	double Lt;

	//calculate for throats
	for (int i = 0; i < network.throats.size(); i++)
	{
		EP1.clear(), EP2.clear();
		EP1 = network.endpoints_throat[i][0];
		EP2 = network.endpoints_throat[i][1];

		Lt = sqrt(pow(EP1[0] - EP2[0], 2) + pow(EP1[1] - EP2[1], 2) + pow(EP1[2] - EP2[2], 2));
		network.length_throat.push_back(Lt);
	}
}

void classPNM::classMethods::classMethodsNetwork::conduit_length_throat_def(classPNM::classNetwork& network)
{
	network.conduit_length_throat.clear();
	vector<double> C1, C2, EP1, EP2; //temp pore cooredinates and throat endpoints
	int p1, p2; //pore indexes for each throat
	double Lt, L1, L2;
	vector<double> temp; //temp storage for conduit lengths {pore1, hroat, pore2}

	//calculate for throats
	for (int i = 0; i < network.throats.size(); i++)
	{
		C1.clear(); C2.clear(); EP1.clear(), EP2.clear();
		temp.clear();
		p1 = network.throats[i][0];
		p2 = network.throats[i][1];

		C1 = network.pores[p1];
		C2 = network.pores[p2];
		EP1 = network.endpoints_throat[i][0];
		EP2 = network.endpoints_throat[i][1];
		Lt = network.length_throat[i];

		L1 = sqrt(pow(C1[0] - EP1[0], 2) + pow(C1[1] - EP1[1], 2) + pow(C1[2] - EP1[2], 2));
		L2 = sqrt(pow(C2[0] - EP2[0], 2) + pow(C2[1] - EP2[1], 2) + pow(C2[2] - EP2[2], 2));

		Lt = max(Lt, 1e-15);

		temp = { L1, Lt, L2 };
		network.conduit_length_throat.push_back(temp);
	}
}

void classPNM::classMethods::classMethodsNetwork::shape_factor_cone_stick(classPNM::classNetwork& network)
{
	network.shape_factor_throat.clear();
	double A1, At, A2, d1, dt, d2, L1, Lt, L2; //pore and throat area, diameter and conduit length
	int p1, p2; //pore indexes for each throat
	double SFt, SF1, SF2; //shape folw factors for pores and throat
	vector<double> temp; //temp storage for shape factors {pore1, hroat, pore2}

	double F1, Ft, F2; //volume
	const double pi = 3.1415926;

	//calculate for throats
	for (int i = 0; i < network.throats.size(); i++)
	{
		temp.clear();
		p1 = network.throats[i][0];
		p2 = network.throats[i][1];

		//Getting equivalent areas
		A1 = network.area_pore[p1];
		At = network.area_troat[i];
		A2 = network.area_pore[p2];
		//getting diameter
		d1 = network.diameter_pore[p1];
		dt = network.diameter_throat[i];
		d2 = network.diameter_pore[p2];
		//Getting conduit lengths
		L1 = network.conduit_length_throat[i][0];
		Lt = network.conduit_length_throat[i][1];
		L2 = network.conduit_length_throat[i][2];

		if (L1 <= 0) { SF1 = 1; }
		else
		{
			F1 = 16.0 / 3.0 * (L1 * (pow(d1, 2) + d1 * dt + pow(dt, 2)) / (pow(d1, 3)*pow(dt, 3)*pow(pi, 2)));
			SF1 = L1 / (pow(A1, 2)*F1);
		}
		if (Lt <= 0) { SFt = 1; }
		else
		{
			Ft = Lt / pow(At, 2);
			SFt = Lt / (pow(At, 2)*Ft);
		}
		if (L2 <= 0) { SF2 = 1; }
		else
		{
			F2 = 16.0 / 3.0 * (L2 * (pow(d2, 2) + d2 * dt + pow(dt, 2)) / (pow(d2, 3)*pow(dt, 3)*pow(pi, 2)));
			SF2 = L2 / (pow(A2, 2)*F2);
		}

		temp = { SF1, SFt, SF2 };
		network.shape_factor_throat.push_back(temp);
	}
}

void classPNM::classMethods::classMethodsNetwork::generate_network_properties(classPNM::classNetwork& network, classPNM::classMethods::classMethodsNetwork& methodsNetwork)
{
	methodsNetwork.cross_sectional_area_simple(network);
	methodsNetwork.endpoints_spherical(network);
	methodsNetwork.length_throat_piecewise(network);
	methodsNetwork.conduit_length_throat_def(network);
	methodsNetwork.shape_factor_cone_stick(network);
	//utils.writeLine("Network Properties calculated");
	cout << "Network Properties calculated" << endl;
}