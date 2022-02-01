#include "classPNM.h"

//
//PHYSICS
//

//classPNM classPNM::classPhysics::hydraulic_conductance_Ma_2014(classPNM PNM, classUtils utils)
//{
//	PNM.physics.hydraulic_conductance_throat.clear();
//
//	int p1, p2; //pore indices 
//	double A1, At, A2, d1, dt, d2, L1, Lt, L2, T1, Tt, T2, P1, Pt, P2, rho1, rhot, rho2, D1, Dt, D2, k1, kt, k2, SF1, SFt, SF2; //input variables
//	double g1, gt, g2; //individual conductances
//
//	//internal consts
//	double alpha = 1; //TMAC (0;1]
//	double pi = 3.1415926;
//
//	//calculate for throats
//	for (int i = 0; i < PNM.network.throats.size(); i++)
//	{
//
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
//		//Getting throat shape factors
//		SF1 = PNM.network.shape_factor_throat[i][0];
//		SFt = PNM.network.shape_factor_throat[i][1];
//		SF2 = PNM.network.shape_factor_throat[i][2];
//		//getting temperature
//		T1 = PNM.phase.temperature_pore[p1];
//		Tt = PNM.phase.temperature_throat[i];
//		T2 = PNM.phase.temperature_pore[p2];
//		//getting pressure
//		P1 = PNM.phase.pressure_pore[p1];
//		Pt = PNM.phase.pressure_throat[i];
//		P2 = PNM.phase.pressure_pore[p2];
//		//getting density
//		rho1 = PNM.phase.density_pore[p1];
//		rhot = PNM.phase.density_throat[i];
//		rho2 = PNM.phase.density_pore[p2];
//		//getting viscosity 
//		D1 = PNM.phase.viscosity_pore[p1];
//		Dt = PNM.phase.viscosity_throat[i];
//		D2 = PNM.phase.viscosity_pore[p2];
//		//getting knudsen number
//		k1 = PNM.phase.knudsen_pore[p1];
//		kt = PNM.phase.knudsen_throat[i];
//		k2 = PNM.phase.knudsen_pore[p2];
//
//		//if Li = 0, gi is infinity
//		if (L1 <= 0) { g1 = numeric_limits<double>::infinity(); }
//		else { g1 = (A1 * pow((d1 / 2), 2) / (8 * D1 * L1)) * ((1 + 1) / 2 + 16 * k1 * ((2 - alpha) / alpha) / 3 + (16 * 16 / (3 * 3 * pi)) * k1); }
//		if (Lt <= 0) { gt = numeric_limits<double>::infinity(); }
//		else { gt = (At * pow((dt / 2), 2) / (8 * Dt * Lt)) * ((1 + 1) / 2 + 16 * kt * ((2 - alpha) / alpha) / 3 + (16 * 16 / (3 * 3 * pi)) * kt); }
//		if (L2 <= 0) { g2 = numeric_limits<double>::infinity(); }
//		else { g2 = (A2 * pow((d2 / 2), 2) / (8 * D2 * L2)) * ((1 + 1) / 2 + 16 * k2 * ((2 - alpha) / alpha) / 3 + (16 * 16 / (3 * 3 * pi)) * k2); }
//
//		PNM.physics.hydraulic_conductance_throat.push_back((L1 + Lt + L2) / (L1 / g1 / SF1 + Lt / gt / SFt + L2 / g2 / SF2));
//
//	}
//
//	return PNM;
//}
//
//classPNM classPNM::classPhysics::hydraulic_conductance_Song_2018(classPNM PNM, classUtils utils)
//{
//	PNM.physics.hydraulic_conductance_throat.clear();
//
//	int p1, p2; //pore indices 
//	double A1, At, A2, d1, dt, d2, L1, Lt, L2, T1, Tt, T2, P1, Pt, P2, rho1, rhot, rho2, D1, Dt, D2, k1, kt, k2, SF1, SFt, SF2; //input variables
//	double g1, gt, g2; //individual conductances
//
//	//internal consts
//	double alpha = 1; //TMAC (0;1]
//	double pi = 3.1415926;
//
//	//calculate for throats
//	for (int i = 0; i < PNM.network.throats.size(); i++)
//	{
//
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
//		//Getting throat shape factors
//		SF1 = PNM.network.shape_factor_throat[i][0];
//		SFt = PNM.network.shape_factor_throat[i][1];
//		SF2 = PNM.network.shape_factor_throat[i][2];
//		//getting temperature
//		T1 = PNM.phase.temperature_pore[p1];
//		Tt = PNM.phase.temperature_throat[i];
//		T2 = PNM.phase.temperature_pore[p2];
//		//getting pressure
//		P1 = PNM.phase.pressure_pore[p1];
//		Pt = PNM.phase.pressure_throat[i];
//		P2 = PNM.phase.pressure_pore[p2];
//		//getting density
//		rho1 = PNM.phase.density_pore[p1];
//		rhot = PNM.phase.density_throat[i];
//		rho2 = PNM.phase.density_pore[p2];
//		//getting viscosity 
//		D1 = PNM.phase.viscosity_pore[p1];
//		Dt = PNM.phase.viscosity_throat[i];
//		D2 = PNM.phase.viscosity_pore[p2];
//		//getting knudsen number
//		k1 = PNM.phase.knudsen_pore[p1];
//		kt = PNM.phase.knudsen_throat[i];
//		k2 = PNM.phase.knudsen_pore[p2];
//
//		//if Li = 0, gi is infinity
//		/*if (L1 <= 0) { g1 = numeric_limits<double>::infinity(); }
//		else { g1 = (A1 * pow((d1 / 2), 2) / (8 * D1 * L1)) * ((1 + 1) / 2 + 16 * k1 * ((2 - alpha) / alpha) / 3 + (16 * 16 / (3 * 3 * pi)) * k1); }
//		if (Lt <= 0) { gt = numeric_limits<double>::infinity(); }
//		else { gt = (At * pow((dt / 2), 2) / (8 * Dt * Lt)) * ((1 + 1) / 2 + 16 * kt * ((2 - alpha) / alpha) / 3 + (16 * 16 / (3 * 3 * pi)) * kt); }
//		if (L2 <= 0) { g2 = numeric_limits<double>::infinity(); }
//		else { g2 = (A2 * pow((d2 / 2), 2) / (8 * D2 * L2)) * ((1 + 1) / 2 + 16 * k2 * ((2 - alpha) / alpha) / 3 + (16 * 16 / (3 * 3 * pi)) * k2); }*/
//		
//		//Song et al 2018
//		if (L1 <= 0) { g1 = numeric_limits<double>::infinity(); }
//		else 
//		{ 
//			g1 = (pi * pow(d1 / 2, 4) / (8 * D1 * L1)) * (1 + 128 * atan(4.0 * pow(k1, 0.4)) * k1 / (15 * pow(pi, 2))) * (1 + 4 * k1 / (1 + k1));
//		}
//		if (Lt <= 0) { gt = numeric_limits<double>::infinity(); }
//		else 
//		{ 
//			gt = (pi * pow(dt / 2, 4) / (8 * Dt * Lt)) * (1 + 128 * atan(4.0 * pow(kt, 0.4)) * kt / (15 * pow(pi, 2))) * (1 + 4 * kt / (1 + kt));
//		}
//		if (L2 <= 0) { g2 = numeric_limits<double>::infinity(); }
//		else 
//		{ 
//			g2 = (pi * pow(d2 / 2, 4) / (8 * D2 * L2)) * (1 + 128 * atan(4.0 * pow(k2, 0.4)) * k2 / (15 * pow(pi, 2))) * (1 + 4 * k2 / (1 + k2));
//		}
//
//		PNM.physics.hydraulic_conductance_throat.push_back((L1 + Lt + L2) / (L1 / g1 / SF1 + Lt / gt / SFt + L2 / g2 / SF2));
//
//	}
//
//	return PNM;
//}
//
//classPNM classPNM::classPhysics::generate_physics_properties(classPNM PNM, classUtils utils)
//{
//	//PNM = PNM.physics.hydraulic_conductance_Ma_2014(PNM, utils);
//	PNM = PNM.physics.hydraulic_conductance_Song_2018(PNM, utils);
//	utils.writeLine("Physics Properties calculated");
//	//utils.emptyLine();
//	return PNM;
//}

//

void classPNM::classMethods::classMethodsPhysics::hydraulic_conductance_Ma_2014(classPNM::classNetwork& network, classPNM::classMethane& phase, classPNM::classPhysics& physics)
{
	physics.hydraulic_conductance_throat.clear();

	int p1, p2; //pore indices 
	double A1, At, A2, d1, dt, d2, L1, Lt, L2, T1, Tt, T2, P1, Pt, P2, rho1, rhot, rho2, D1, Dt, D2, k1, kt, k2, SF1, SFt, SF2; //input variables
	double g1, gt, g2; //individual conductances

	//internal consts
	double alpha = 1; //TMAC (0;1]
	double pi = 3.1415926;

	//calculate for throats
	for (int i = 0; i < network.throats.size(); i++)
	{

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
		//Getting throat shape factors
		SF1 = network.shape_factor_throat[i][0];
		SFt = network.shape_factor_throat[i][1];
		SF2 = network.shape_factor_throat[i][2];
		//getting temperature
		T1 = phase.temperature_pore[p1];
		Tt = phase.temperature_throat[i];
		T2 = phase.temperature_pore[p2];
		//getting pressure
		P1 = phase.pressure_pore[p1];
		Pt = phase.pressure_throat[i];
		P2 = phase.pressure_pore[p2];
		//getting density
		rho1 = phase.density_pore[p1];
		rhot = phase.density_throat[i];
		rho2 = phase.density_pore[p2];
		//getting viscosity 
		D1 = phase.viscosity_pore[p1];
		Dt = phase.viscosity_throat[i];
		D2 = phase.viscosity_pore[p2];
		//getting knudsen number
		k1 = phase.knudsen_pore[p1];
		kt = phase.knudsen_throat[i];
		k2 = phase.knudsen_pore[p2];

		//if Li = 0, gi is infinity
		if (L1 <= 0) { g1 = numeric_limits<double>::infinity(); }
		else { g1 = (A1 * pow((d1 / 2), 2) / (8 * D1 * L1)) * ((1 + 1) / 2 + 16 * k1 * ((2 - alpha) / alpha) / 3 + (16 * 16 / (3 * 3 * pi)) * k1); }
		if (Lt <= 0) { gt = numeric_limits<double>::infinity(); }
		else { gt = (At * pow((dt / 2), 2) / (8 * Dt * Lt)) * ((1 + 1) / 2 + 16 * kt * ((2 - alpha) / alpha) / 3 + (16 * 16 / (3 * 3 * pi)) * kt); }
		if (L2 <= 0) { g2 = numeric_limits<double>::infinity(); }
		else { g2 = (A2 * pow((d2 / 2), 2) / (8 * D2 * L2)) * ((1 + 1) / 2 + 16 * k2 * ((2 - alpha) / alpha) / 3 + (16 * 16 / (3 * 3 * pi)) * k2); }

		physics.hydraulic_conductance_throat.push_back((L1 + Lt + L2) / (L1 / g1 / SF1 + Lt / gt / SFt + L2 / g2 / SF2));

	}
}

void classPNM::classMethods::classMethodsPhysics::hydraulic_conductance_Song_2018(classPNM::classNetwork& network, classPNM::classMethane& phase, classPNM::classPhysics& physics)
{
	physics.hydraulic_conductance_throat.clear();

	int p1, p2; //pore indices 
	double A1, At, A2, d1, dt, d2, L1, Lt, L2, T1, Tt, T2, P1, Pt, P2, rho1, rhot, rho2, D1, Dt, D2, k1, kt, k2, SF1, SFt, SF2; //input variables
	double g1, gt, g2; //individual conductances

	//internal consts
	double alpha = 1; //TMAC (0;1]
	double pi = 3.1415926;

	//calculate for throats
	for (int i = 0; i < network.throats.size(); i++)
	{

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
		//Getting throat shape factors
		SF1 = network.shape_factor_throat[i][0];
		SFt = network.shape_factor_throat[i][1];
		SF2 = network.shape_factor_throat[i][2];
		//getting temperature
		T1 = phase.temperature_pore[p1];
		Tt = phase.temperature_throat[i];
		T2 = phase.temperature_pore[p2];
		//getting pressure
		P1 = phase.pressure_pore[p1];
		Pt = phase.pressure_throat[i];
		P2 = phase.pressure_pore[p2];
		//getting density
		rho1 = phase.density_pore[p1];
		rhot = phase.density_throat[i];
		rho2 = phase.density_pore[p2];
		//getting viscosity 
		D1 = phase.viscosity_pore[p1];
		Dt = phase.viscosity_throat[i];
		D2 = phase.viscosity_pore[p2];
		//getting knudsen number
		k1 = phase.knudsen_pore[p1];
		kt = phase.knudsen_throat[i];
		k2 = phase.knudsen_pore[p2];

		//if Li = 0, gi is infinity
		/*if (L1 <= 0) { g1 = numeric_limits<double>::infinity(); }
		else { g1 = (A1 * pow((d1 / 2), 2) / (8 * D1 * L1)) * ((1 + 1) / 2 + 16 * k1 * ((2 - alpha) / alpha) / 3 + (16 * 16 / (3 * 3 * pi)) * k1); }
		if (Lt <= 0) { gt = numeric_limits<double>::infinity(); }
		else { gt = (At * pow((dt / 2), 2) / (8 * Dt * Lt)) * ((1 + 1) / 2 + 16 * kt * ((2 - alpha) / alpha) / 3 + (16 * 16 / (3 * 3 * pi)) * kt); }
		if (L2 <= 0) { g2 = numeric_limits<double>::infinity(); }
		else { g2 = (A2 * pow((d2 / 2), 2) / (8 * D2 * L2)) * ((1 + 1) / 2 + 16 * k2 * ((2 - alpha) / alpha) / 3 + (16 * 16 / (3 * 3 * pi)) * k2); }*/

		//Song et al 2018
		if (L1 <= 0) { g1 = numeric_limits<double>::infinity(); }
		else
		{
			g1 = (pi * pow(d1 / 2, 4) / (8 * D1 * L1)) * (1 + 128 * atan(4.0 * pow(k1, 0.4)) * k1 / (15 * pow(pi, 2))) * (1 + 4 * k1 / (1 + k1));
		}
		if (Lt <= 0) { gt = numeric_limits<double>::infinity(); }
		else
		{
			gt = (pi * pow(dt / 2, 4) / (8 * Dt * Lt)) * (1 + 128 * atan(4.0 * pow(kt, 0.4)) * kt / (15 * pow(pi, 2))) * (1 + 4 * kt / (1 + kt));
		}
		if (L2 <= 0) { g2 = numeric_limits<double>::infinity(); }
		else
		{
			g2 = (pi * pow(d2 / 2, 4) / (8 * D2 * L2)) * (1 + 128 * atan(4.0 * pow(k2, 0.4)) * k2 / (15 * pow(pi, 2))) * (1 + 4 * k2 / (1 + k2));
		}

		physics.hydraulic_conductance_throat.push_back((L1 + Lt + L2) / (L1 / g1 / SF1 + Lt / gt / SFt + L2 / g2 / SF2));

	}
}

void classPNM::classMethods::classMethodsPhysics::generate_physics_properties(classPNM::classNetwork& network, classPNM::classMethane& phase, classPNM::classPhysics& physics, classPNM::classMethods::classMethodsPhysics& methodsPhysics)
{
	//methodsPhysics.hydraulic_conductance_Ma_2014(network, phase, physics)
	methodsPhysics.hydraulic_conductance_Song_2018(network, phase, physics);
	//utils.writeLine("Physics Properties calculated");
	cout << "Physics Properties calculated" << endl;
	//utils.emptyLine();
}