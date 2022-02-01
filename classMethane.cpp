#include "classPNM.h"

double Cardano_min(double a, double b, double c, double d)
{
	double Pi = 3.1415926;
	double p = (3 * a * c - pow(b, 2)) / 3 / pow(a, 2);
	double q = 2 * pow(b, 3) / 27 / pow(a, 3) - b * c / 3 / pow(a, 2) + d / a;
	double discrim = pow((q / 2), 2) + pow((p / 3), 3);
	double cardano;

	if (discrim < 0)
	{
		double phi = acos(-q / 2 / sqrt(pow((-p / 3), 3)));
		double y1 = 2 * sqrt(-p / 3) * cos(phi / 3);
		double y2 = 2 * sqrt(-p / 3) * cos(phi / 3 + 2 * Pi / 3);
		double y3 = 2 * sqrt(-p / 3) * cos(phi / 3 + 4 * Pi / 3);

		double x1 = y1 - b / 3 / a;
		cardano = x1;
		double x2 = y2 - b / 3 / a;
		if (x2 < x1)
		{
			cardano = x2;
		}
		double x3 = y3 - b / 3 / a;
		if (x3 < x2)
		{
			cardano = x3;
		}
	}
	else if (discrim > 0)
	{
		double u = pow((-q / 2 + sqrt(discrim)), (1.0 / 3));
		double v = pow((-q / 2 - sqrt(discrim)), (1.0 / 3));

		double y1 = u + v;
		double x1 = y1 - b / 3 / a;
		cardano = x1;
	}
	else
	{
		double y1 = pow(-4 * q, 1.0 / 3);
		double y2 = pow(q / 2, 1.0 / 3);

		double x1 = y1 - b / 3 / a;
		cardano = x1;
		double x2 = y2 - b / 3 / a;
		if (x2 < x1)
		{
			cardano = x2;
		}
	}

	return cardano;
}

//
//METHANE
//

//classPNM classPNM::classMethane::density_Peng_Robinson_Methane(classPNM PNM, classUtils utils)
//{
//
//
//	PNM.phase.density_pore.clear();
//	PNM.phase.density_throat.clear();
//
//	double R = 8.3145;              //R = gas constant[J / mol / K]
//	double omega = 0.01142;       //omega value of CH4
//
//	double f_omega = 0.37464 + 1.54226 * omega - 0.26992 * (pow(omega, 2));
//
//	//calculate for pores
//	for (int i = 0; i < PNM.network.pores.size(); i++)
//	{
//		double Tr = PNM.phase.temperature_pore[i] / PNM.phase.Tc;
//		double a = 0.45724 * pow(R, 2) * pow(PNM.phase.Tc, 2) / PNM.phase.Pc * pow((1 + f_omega * (1 - pow(Tr, 0.5))), 2);
//		double b = 0.0778 * R * PNM.phase.Tc / PNM.phase.Pc;
//
//		double Af = a * PNM.phase.pressure_pore[i] / pow(R, 2) / pow(PNM.phase.temperature_pore[i], 2);
//		double Bf = b * PNM.phase.pressure_pore[i] / R / PNM.phase.temperature_pore[i];
//
//		//factors for virial equation of Z : v1 * Z ^ 3 + v2 * Z ^ 2 + v3 * Z + v4
//		double v1 = 1;
//		double v2 = -(1 - Bf);
//		double v3 = (Af - 3 * pow(Bf, 2) - 2 * Bf);
//		double v4 = -(Af * Bf - pow(Bf, 2) - pow(Bf, 3));
//
//		PNM.phase.density_pore.push_back(PNM.phase.pressure_pore[i] * PNM.phase.MW / (utils.Cardano_min(v1, v2, v3, v4) * R * PNM.phase.temperature_pore[i]));
//	}
//	//calculate for throats
//	for (int i = 0; i < PNM.network.throats.size(); i++)
//	{
//		double Tr = PNM.phase.temperature_throat[i] / PNM.phase.Tc;
//		double a = 0.45724 * pow(R, 2) * pow(PNM.phase.Tc, 2) / PNM.phase.Pc * pow((1 + f_omega * (1 - pow(Tr, 0.5))), 2);
//		double b = 0.0778 * R * PNM.phase.Tc / PNM.phase.Pc;
//
//		double Af = a * PNM.phase.pressure_throat[i] / pow(R, 2) / pow(PNM.phase.temperature_throat[i], 2);
//		double Bf = b * PNM.phase.pressure_throat[i] / R / PNM.phase.temperature_throat[i];
//
//		//factors for virial equation of Z : v1 * Z ^ 3 + v2 * Z ^ 2 + v3 * Z + v4
//		double v1 = 1;
//		double v2 = -(1 - Bf);
//		double v3 = (Af - 3 * pow(Bf, 2) - 2 * Bf);
//		double v4 = -(Af * Bf - pow(Bf, 2) - pow(Bf, 3));
//
//		PNM.phase.density_throat.push_back(PNM.phase.pressure_throat[i] * PNM.phase.MW / (utils.Cardano_min(v1, v2, v3, v4) * R * PNM.phase.temperature_throat[i]));
//	}
//
//	return PNM;
//}
//
//classPNM classPNM::classMethane::viscosity_natural_gas_Lee(classPNM PNM, classUtils utils)
//{
//	PNM.phase.viscosity_pore.clear();
//	PNM.phase.viscosity_throat.clear();
//
//	//calculate for pores
//	for (int i = 0; i < PNM.network.pores.size(); i++)
//	{
//		double K = (9.379 + 0.01607*PNM.phase.MW * 1000)*pow(PNM.phase.temperature_pore[i] * 1.8, 1.5) / (209.2 + 19.26*PNM.phase.MW * 1000 + PNM.phase.temperature_pore[i] * 1.8);
//		double X = 3.448 + 986.4 / (PNM.phase.temperature_pore[i] * 1.8) + 0.01009*PNM.phase.MW * 1000;
//		double Y = 2.447 - 0.2224*X;
//
//		PNM.phase.viscosity_pore.push_back(1e-7 * K * exp(X*pow(PNM.phase.density_pore[i] / 1000, Y)));
//	}
//	//calculate for throats
//	for (int i = 0; i < PNM.network.throats.size(); i++)
//	{
//		double K = (9.379 + 0.01607*PNM.phase.MW * 1000)*pow(PNM.phase.temperature_throat[i] * 1.8, 1.5) / (209.2 + 19.26*PNM.phase.MW * 1000 + PNM.phase.temperature_throat[i] * 1.8);
//		double X = 3.448 + 986.4 / (PNM.phase.temperature_throat[i] * 1.8) + 0.01009*PNM.phase.MW * 1000;
//		double Y = 2.447 - 0.2224*X;
//
//		PNM.phase.viscosity_throat.push_back(1e-7 * K * exp(X*pow(PNM.phase.density_throat[i] / 1000, Y)));
//	}
//	return PNM;
//}
//
//classPNM classPNM::classMethane::knudsen_number(classPNM PNM, classUtils utils)
//{
//	PNM.phase.knudsen_pore.clear();
//	PNM.phase.knudsen_throat.clear();
//
//	double kb = 1.380649e-23;  //J / K
//
//	double pi = 3.1415926;
//
//	//calculate for pores
//	for (int i = 0; i < PNM.network.pores.size(); i++)
//	{
//		PNM.phase.knudsen_pore.push_back((kb * PNM.phase.temperature_pore[i] / (1.414214 * pi * pow(PNM.phase.d_mol, 2) * PNM.phase.pressure_pore[i])) / PNM.network.diameter_pore[i]);
//	}
//	//calculate for throats
//	for (int i = 0; i < PNM.network.throats.size(); i++)
//	{
//		PNM.phase.knudsen_throat.push_back((kb * PNM.phase.temperature_throat[i] / (1.414214 * pi * pow(PNM.phase.d_mol, 2) * PNM.phase.pressure_throat[i])) / PNM.network.diameter_throat[i]);
//	}
//	return PNM;
//}
//
//classPNM classPNM::classMethane::generate_phase_properties(classPNM PNM, classUtils utils)
//{
//	PNM = PNM.phase.density_Peng_Robinson_Methane(PNM, utils);
//	PNM = PNM.phase.viscosity_natural_gas_Lee(PNM, utils);
//	PNM = PNM.phase.knudsen_number(PNM, utils);
//	utils.writeLine("Phase Properties calculated");
//	//utils.emptyLine();
//
//	return PNM;
//}

//

void classPNM::classMethods::classMethodsMethane::density_Peng_Robinson_Methane(classPNM::classMethane& phase, classPNM::classNetwork& network)
{
	phase.density_pore.clear();
	phase.density_throat.clear();

	double R = 8.3145;              //R = gas constant[J / mol / K]
	double omega = 0.01142;       //omega value of CH4

	double f_omega = 0.37464 + 1.54226 * omega - 0.26992 * (pow(omega, 2));

	//calculate for pores
	for (int i = 0; i < network.pores.size(); i++)
	{
		double Tr = phase.temperature_pore[i] / phase.Tc;
		double a = 0.45724 * pow(R, 2) * pow(phase.Tc, 2) / phase.Pc * pow((1 + f_omega * (1 - pow(Tr, 0.5))), 2);
		double b = 0.0778 * R * phase.Tc / phase.Pc;

		double Af = a * phase.pressure_pore[i] / pow(R, 2) / pow(phase.temperature_pore[i], 2);
		double Bf = b * phase.pressure_pore[i] / R / phase.temperature_pore[i];

		//factors for virial equation of Z : v1 * Z ^ 3 + v2 * Z ^ 2 + v3 * Z + v4
		double v1 = 1;
		double v2 = -(1 - Bf);
		double v3 = (Af - 3 * pow(Bf, 2) - 2 * Bf);
		double v4 = -(Af * Bf - pow(Bf, 2) - pow(Bf, 3));

		phase.density_pore.push_back(phase.pressure_pore[i] * phase.MW / (Cardano_min(v1, v2, v3, v4) * R * phase.temperature_pore[i]));
	}
	//calculate for throats
	for (int i = 0; i < network.throats.size(); i++)
	{
		double Tr = phase.temperature_throat[i] / phase.Tc;
		double a = 0.45724 * pow(R, 2) * pow(phase.Tc, 2) / phase.Pc * pow((1 + f_omega * (1 - pow(Tr, 0.5))), 2);
		double b = 0.0778 * R * phase.Tc / phase.Pc;

		double Af = a * phase.pressure_throat[i] / pow(R, 2) / pow(phase.temperature_throat[i], 2);
		double Bf = b * phase.pressure_throat[i] / R / phase.temperature_throat[i];

		//factors for virial equation of Z : v1 * Z ^ 3 + v2 * Z ^ 2 + v3 * Z + v4
		double v1 = 1;
		double v2 = -(1 - Bf);
		double v3 = (Af - 3 * pow(Bf, 2) - 2 * Bf);
		double v4 = -(Af * Bf - pow(Bf, 2) - pow(Bf, 3));

		phase.density_throat.push_back(phase.pressure_throat[i] * phase.MW / (Cardano_min(v1, v2, v3, v4) * R * phase.temperature_throat[i]));
	}
}

void classPNM::classMethods::classMethodsMethane::viscosity_natural_gas_Lee(classPNM::classMethane& phase, classPNM::classNetwork& network)
{
	phase.viscosity_pore.clear();
	phase.viscosity_throat.clear();

	//calculate for pores
	for (int i = 0; i < network.pores.size(); i++)
	{
		double K = (9.379 + 0.01607*phase.MW * 1000)*pow(phase.temperature_pore[i] * 1.8, 1.5) / (209.2 + 19.26*phase.MW * 1000 + phase.temperature_pore[i] * 1.8);
		double X = 3.448 + 986.4 / (phase.temperature_pore[i] * 1.8) + 0.01009*phase.MW * 1000;
		double Y = 2.447 - 0.2224*X;

		phase.viscosity_pore.push_back(1e-7 * K * exp(X*pow(phase.density_pore[i] / 1000, Y)));
	}
	//calculate for throats
	for (int i = 0; i < network.throats.size(); i++)
	{
		double K = (9.379 + 0.01607*phase.MW * 1000)*pow(phase.temperature_throat[i] * 1.8, 1.5) / (209.2 + 19.26*phase.MW * 1000 + phase.temperature_throat[i] * 1.8);
		double X = 3.448 + 986.4 / (phase.temperature_throat[i] * 1.8) + 0.01009*phase.MW * 1000;
		double Y = 2.447 - 0.2224*X;

		phase.viscosity_throat.push_back(1e-7 * K * exp(X*pow(phase.density_throat[i] / 1000, Y)));
	}
}

void classPNM::classMethods::classMethodsMethane::knudsen_number(classPNM::classMethane& phase, classPNM::classNetwork& network)
{
	phase.knudsen_pore.clear();
	phase.knudsen_throat.clear();

	double kb = 1.380649e-23;  //J / K

	double pi = 3.1415926;

	//calculate for pores
	for (int i = 0; i < network.pores.size(); i++)
	{
		phase.knudsen_pore.push_back((kb * phase.temperature_pore[i] / (1.414214 * pi * pow(phase.d_mol, 2) * phase.pressure_pore[i])) / network.diameter_pore[i]);
	}
	//calculate for throats
	for (int i = 0; i < network.throats.size(); i++)
	{
		phase.knudsen_throat.push_back((kb * phase.temperature_throat[i] / (1.414214 * pi * pow(phase.d_mol, 2) * phase.pressure_throat[i])) / network.diameter_throat[i]);
	}
}

void classPNM::classMethods::classMethodsMethane::generate_phase_properties(classPNM::classNetwork& network,classPNM::classMethane& phase ,classPNM::classMethods::classMethodsMethane& methodsPhase)
{
	methodsPhase.density_Peng_Robinson_Methane(phase, network);
	methodsPhase.viscosity_natural_gas_Lee(phase, network);
	methodsPhase.knudsen_number(phase, network);
	cout << "Phase Properties calculated" << endl;

}