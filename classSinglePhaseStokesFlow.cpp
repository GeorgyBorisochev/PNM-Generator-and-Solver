#include "classPNM.h"

//
//ALGORYTHM
//

template <typename output>
void printVector(vector<output> vector)
{
	for (int i = 0; i < vector.size(); i++)
	{
		cout << vector[i] << " ";
	}
	cout << endl;
}

//classPNM classPNM::classSinglePhaseStokesFlow::create_A_b(classPNM PNM, classUtils utils)
//{
//	//pass necessary properties to algorythm
//	PNM.algorithm.conductance = PNM.physics.hydraulic_conductance_throat;
//	//define temporary A and b
//	int Np = PNM.network.pores.size();
//	compressed2D<double> A(Np, Np);
//	dense_vector<double> b(Np);
//	b = 0.0;
//
//	//global P_BC vector
//	dense_vector<double> pressure_BC(Np);
//	for (int i = 0; i < Np; i++)
//	{
//		if (PNM.algorithm.inlets[i])
//		{
//			pressure_BC[i] = PNM.algorithm.Pin;
//		}
//		else if (PNM.algorithm.outlets[i])
//		{
//			pressure_BC[i] = (PNM.algorithm.Pin - PNM.algorithm.dP);
//		}
//		else
//		{
//			pressure_BC[i] = 0.0;
//		}
//	}
//	PNM.algorithm.pressure_BC = pressure_BC;
//	//cout << "pressure_BC = \n" << pressure_BC << "\n\n";
//
//	//fill A and b in a different fuction - inserter has to be destroyed
//	b = PNM.algorithm.populate_A_b(A, b, PNM);
//	//cout << "A = \n" << A << "\n\n";
//	//cout << "b = \n" << b << "\n\n";
//
//	//pass A and b to class
//	PNM.algorithm.A = A;
//	PNM.algorithm.b = b;
//
//	return PNM;
//}
//
//dense_vector<double> classPNM::classSinglePhaseStokesFlow::populate_A_b(compressed2D<double>& A, dense_vector<double> b, classPNM PNM)
//{
//	A = 0.0;
//	int Np = PNM.network.pores.size();
//	int Nt = PNM.network.throats.size();
//	int head, tail; //pore indices for each throat
//	double averageDiagonal = 0;
//	double minA, maxA, minb, maxb;
//	minA = maxA = minb = maxb = 0.0;
//	// Create incremental plus inserter for matrix A
//	//mat::inserter<compressed2D<double>, update_plus<double>> ins(A);
//	mat::inserter<compressed2D<double>> ins(A);
//	//laplacian matirx to be created first! 
//	compressed2D<double> laplacian_A(Np, Np);
//
//	//1. create boundary condition pressure vector
//	dense_vector<double> pressure_BC(Np);
//	for (int i = 0; i < Np; i++)
//	{
//		if (PNM.algorithm.inlets[i])
//		{
//			pressure_BC[i] = PNM.algorithm.Pin;
//		}
//		else if (PNM.algorithm.outlets[i])
//		{
//			pressure_BC[i] = (PNM.algorithm.Pin - PNM.algorithm.dP);
//		}
//		else
//		{
//			pressure_BC[i] = 0.0;
//		}
//	}
//	//cout << "pressure_BC = \n" << pressure_BC << "\n\n";
//
//
//	//2. complete laplacian ajaicency marix  
//	averageDiagonal = PNM.algorithm.laplacian_A(laplacian_A, PNM);
//	//cout << "average Diagonal = " << averageDiagonal << "\n";
//	//cout << "laplacian_A = \n" << laplacian_A << "\n\n";
//
//	//3. multiply laplacian by P_BC for non-BC entries of B
//	dense_vector<double> A_BC(Np);
//	A_BC = laplacian_A * pressure_BC;
//	//cout << "A_BC = \n" << A_BC << "\n\n";
//
//	//4. create A: laplacian for non-BC entries, zero for BC entries
//	//diagonal filled in a separate pass
//	for (int i = 0; i < Nt; i++)
//	{
//		head = PNM.network.throats[i][0];
//		tail = PNM.network.throats[i][1];
//		//check if non-BC for head and tail pores: 
//		if (!PNM.algorithm.inlets[head] && !PNM.algorithm.outlets[head])
//		{
//			ins[head][head] << laplacian_A[head][head];
//			//ins[head][tail] << laplacian_A[head][tail];
//			//if (laplacian_A[head][head] < minA) { minA = laplacian_A[head][head]; }
//			//if (laplacian_A[head][head] > maxA) { maxA = laplacian_A[head][head]; }
//			//if (laplacian_A[head][head] > 1.0) { void __debugbreak(); }
//		}
//		if (!PNM.algorithm.inlets[tail] && !PNM.algorithm.outlets[tail])
//		{
//			ins[tail][tail] << laplacian_A[tail][tail];
//			//ins[tail][head] << laplacian_A[tail][head];
//			//if (laplacian_A[tail][tail] < minA) { minA = laplacian_A[tail][tail]; }
//			//if (laplacian_A[tail][tail] > maxA) { maxA = laplacian_A[tail][tail]; }
//			//if (laplacian_A[tail][tail] > 1.0) { void __debugbreak(); }
//		}
//		if ((!PNM.algorithm.inlets[head] && !PNM.algorithm.outlets[head]) && (!PNM.algorithm.inlets[tail] && !PNM.algorithm.outlets[tail]))
//		{
//			ins[head][tail] << laplacian_A[head][tail];
//			ins[tail][head] << laplacian_A[tail][head];
//
//			//if (laplacian_A[head][tail] < minA) { minA = laplacian_A[head][tail]; }
//			//if (laplacian_A[head][tail] > maxA) { maxA = laplacian_A[head][tail]; }
//
//			//if (laplacian_A[tail][head] < minA) { minA = laplacian_A[tail][head]; }
//			//if (laplacian_A[tail][head] > maxA) { maxA = laplacian_A[tail][head]; }
//
//			//if (laplacian_A[head][tail]  > 1.0) { void __debugbreak(); }
//			//if (laplacian_A[tail][head]  > 1.0) { void __debugbreak(); }
//
//		}
//		//if (maxA > 1) { cout << i << "\n"; }
//	}
//
//	//4. cont-d: fill BC diagonal values
//	//5. create b: average_diagonal*P_BC for BC entries, -laplacian*P_BC for non-BC entries
//	for (int i = 0; i < Np; i++)
//	{
//		if (PNM.algorithm.inlets[i])
//		{
//			ins[i][i] << averageDiagonal;
//			b[i] = averageDiagonal * PNM.algorithm.Pin;
//
//			//if (b[i] < minb) { minb = b[i]; }
//			//if (b[i] > maxb) { maxb = b[i]; }
//		}
//		else if (PNM.algorithm.outlets[i])
//		{
//			ins[i][i] << averageDiagonal;
//			b[i] = averageDiagonal * (PNM.algorithm.Pin - PNM.algorithm.dP);
//
//			//if (b[i] < minb) { minb = b[i]; }
//			//if (b[i] > maxb) { maxb = b[i]; }
//		}
//		else
//		{
//			b[i] = -A_BC[i];
//
//			//if (b[i] < minb) { minb = b[i]; }
//			//if (b[i] > maxb) { maxb = b[i]; }
//		}
//	}
//	//debug
//	/*cout << "minA = " << minA << "\n";
//	cout << "maxA = " << maxA << "\n";
//	cout << "minb = " << minb << "\n";
//	cout << "maxb = " << maxb << "\n";*/
//
//	return b;
//}
//
//double classPNM::classSinglePhaseStokesFlow::laplacian_A(compressed2D<double>& laplacian_A, classPNM PNM)
//{
//	laplacian_A = 0.0;
//	int Np = PNM.network.pores.size();
//	int Nt = PNM.network.throats.size();
//	int head, tail; //pore indices for each throat
//	double averageDiagonal = 0;
//	// Create inserter for matrix A
//	mat::inserter<compressed2D<double>, update_plus<double>> insert(laplacian_A);
//
//	//complete A in laplacian form without BC 
//	for (int i = 0; i < Nt; i++)
//	{
//		head = PNM.network.throats[i][0];
//		tail = PNM.network.throats[i][1];
//
//		insert[head][head] << PNM.algorithm.conductance[i];
//		insert[head][tail] << -PNM.algorithm.conductance[i];
//
//		insert[tail][tail] << PNM.algorithm.conductance[i];
//		insert[tail][head] << -PNM.algorithm.conductance[i];
//
//		averageDiagonal += 2 * PNM.algorithm.conductance[i];
//	}
//	averageDiagonal /= Np;
//	return averageDiagonal;
//}
//
//classPNM classPNM::classSinglePhaseStokesFlow::solve_A_b(classPNM PNM, classUtils utils)
//{
//	//function for solving matirx equation once
//	//local variables
//	int Np = PNM.network.pores.size();
//	compressed2D<double> A = PNM.algorithm.A;
//	dense_vector<double> b = PNM.algorithm.b;
//	dense_vector<double> X(Np, 0.0);
//	//cout << "b = \n" << b << "\n\n";
//
//	// Create an ILU(0) preconditioner
//	pc::ilu_0<compressed2D<double>> P(A);
//
//	// Termination criterion: r < 1e-6 * b or N iterations
//	cyclic_iteration<double>       iter(b, 5000, 1.0e-8, 0, 100);
//
//	// Solve Ax == b with left preconditioner P
//	bicgstab(A, X, b, P, iter);
//	//cout << "X = \n" << X << "\n\n";
//	PNM.algorithm.X = X;
//
//	return PNM;
//}
//
//classPNM classPNM::classSinglePhaseStokesFlow::pass_result_to_phase(classPNM PNM, classUtils utils)
//{
//	int Np = PNM.network.pores.size();
//	int Nt = PNM.network.throats.size();
//	int head, tail;
//
//	//pass solution vector X to phase pressure
//	//check for pressure values out of bounds - due to unconnected regions
//	for (int i = 0; i < Np; i++)
//	{
//		if (PNM.algorithm.X[i] > PNM.algorithm.Pin) { PNM.phase.pressure_pore[i] = PNM.algorithm.Pin; }
//		else if (PNM.algorithm.X[i] < (PNM.algorithm.Pin - PNM.algorithm.dP)) { PNM.phase.pressure_pore[i] = (PNM.algorithm.Pin - PNM.algorithm.dP); }
//		else
//		{
//			PNM.phase.pressure_pore[i] = PNM.algorithm.X[i];
//		}
//	}
//
//	//interpolate pore pressure for throats
//	for (int i = 0; i < Nt; i++)
//	{
//		head = PNM.network.throats[i][0];
//		tail = PNM.network.throats[i][1];
//		PNM.phase.pressure_throat[i] = 0.5*(PNM.phase.pressure_pore[head] + PNM.phase.pressure_pore[tail]);
//	}
//	return PNM;
//}
//
//classPNM classPNM::classSinglePhaseStokesFlow::calculate_flow_rate(classPNM PNM, classUtils utils)
//{
//	int Np = PNM.network.pores.size();
//	int Nt = PNM.network.throats.size();
//	int head, tail;
//	double Qt;
//
//	PNM.algorithm.flow_rate_pore.clear();
//	PNM.algorithm.flow_rate_throat.clear();
//	vector<double> Qp(Np, 0);
//
//	//calculate flow rate for every pore
//	for (int i = 0; i < Nt; i++)
//	{
//		head = PNM.network.throats[i][0];
//		tail = PNM.network.throats[i][1];
//		Qt = PNM.algorithm.conductance[i] * (PNM.phase.pressure_pore[tail] - PNM.phase.pressure_pore[head]);
//
//		Qp[head] -= Qt;
//		Qp[tail] += Qt;
//		PNM.algorithm.flow_rate_throat.push_back(Qt);
//	}
//	PNM.algorithm.flow_rate_pore = Qp;
//	return PNM;
//}
//
//classPNM classPNM::classSinglePhaseStokesFlow::calculate_effective_permeability(classPNM PNM, classUtils utils)
//{
//	double iterator = 0;
//	double domain_area = pow(PNM.settings.domain, 2);
//
//	int Np = PNM.network.pores.size();
//	//calculate Q*viscosity for each inlet pore
//	for (int i = 0; i < Np; i++)
//	{
//		if (PNM.algorithm.inlets[i])
//		{
//			iterator += PNM.algorithm.flow_rate_pore[i] * PNM.phase.viscosity_pore[i];
//		}
//	}
//
//	//calculate effective permeability in nano-Darcy
//	double temp = iterator * PNM.settings.domain / (domain_area * PNM.algorithm.dP) / 9.86e-12 * 1e9;
//	cout << "Effective Permeability: " << temp << " [nD]";
//	PNM.algorithm.effective_permeability.push_back(temp);
//	utils.emptyLine();
//
//	return PNM;
//}
//
//classPNM classPNM::classSinglePhaseStokesFlow::run_once(classPNM PNM, classUtils utils)
//{
//	PNM = PNM.algorithm.create_A_b(PNM, utils);
//	utils.writeLine("Algorythm Parameters Initialised");
//
//	PNM = PNM.algorithm.solve_A_b(PNM, utils);
//	PNM = PNM.algorithm.pass_result_to_phase(PNM, utils);
//	/*utils.writeLine("Linear equation solved");
//	utils.emptyLine();*/
//	PNM = PNM.algorithm.calculate_flow_rate(PNM, utils);
//	PNM = PNM.algorithm.calculate_effective_permeability(PNM, utils);
//	return PNM;
//}
//
//classPNM classPNM::classSinglePhaseStokesFlow::run_linear_iterative(classPNM PNM, classUtils utils)
//{
//	PNM.algorithm.effective_permeability.clear();
//	int n = 5; //iterations
//	for (int i = 0; i < n; i++)
//	{
//		cout << "\n"<<"Solver Iteration " << i << ":\n";
//		PNM.methods.methodsPhase.generate_phase_properties(PNM.network, PNM.phase, PNM.methods.methodsPhase);
//		PNM.methods.methodsPhysics.generate_physics_properties(PNM.network, PNM.phase, PNM.physics, PNM.methods.methodsPhysics);
//		PNM = PNM.algorithm.run_once(PNM, utils);
//		
//	}
//	return PNM;
//}

//

void classPNM::classMethods::classMethodsSinglePhaseStokesFlow::
			create_A_b(classPNM::classNetwork& network,
						classPNM::classPhysics& physics,
						classPNM::classSinglePhaseStokesFlow& algorithm,
						classPNM::classMethods::classMethodsSinglePhaseStokesFlow& methodsAlgorithm)
{
	//pass necessary properties to algorythm
	algorithm.conductance = physics.hydraulic_conductance_throat;
	//define temporary A and b
	int Np = network.pores.size();
	compressed2D<double> A(Np, Np);
	dense_vector<double> b(Np);
	b = 0.0;

	//global P_BC vector
	dense_vector<double> pressure_BC(Np);
	for (int i = 0; i < Np; i++)
	{
		if (algorithm.inlets[i])
		{
			pressure_BC[i] = algorithm.Pin;
		}
		else if (algorithm.outlets[i])
		{
			pressure_BC[i] = (algorithm.Pin - algorithm.dP);
		}
		else
		{
			pressure_BC[i] = 0.0;
		}
	}
	algorithm.pressure_BC = pressure_BC;
	//cout << "pressure_BC = \n" << pressure_BC << "\n\n";

	//fill A and b in a different fuction - inserter has to be destroyed
	b = methodsAlgorithm.populate_A_b(A, b, network, algorithm, methodsAlgorithm);
	//cout << "A = \n" << A << "\n\n";
	//cout << "b = \n" << b << "\n\n";

	//pass A and b to class
	algorithm.A = A;
	algorithm.b = b;
}

dense_vector<double> classPNM::classMethods::classMethodsSinglePhaseStokesFlow::
						populate_A_b(compressed2D<double>& A, 
										dense_vector<double>& b, 
										classPNM::classNetwork& network,
										classPNM::classSinglePhaseStokesFlow& algorithm,
										classPNM::classMethods::classMethodsSinglePhaseStokesFlow& methodsAlgorithm)
{
	A = 0.0;
	int Np = network.pores.size();
	int Nt = network.throats.size();
	int head, tail; //pore indices for each throat
	double averageDiagonal = 0;
	double minA, maxA, minb, maxb;
	minA = maxA = minb = maxb = 0.0;
	// Create incremental plus inserter for matrix A
	//mat::inserter<compressed2D<double>, update_plus<double>> ins(A);
	mat::inserter<compressed2D<double>> ins(A);
	//laplacian matirx to be created first! 
	compressed2D<double> laplacian_A(Np, Np);

	//1. create boundary condition pressure vector
	dense_vector<double> pressure_BC(Np);
	for (int i = 0; i < Np; i++)
	{
		if (algorithm.inlets[i])
		{
			pressure_BC[i] = algorithm.Pin;
		}
		else if (algorithm.outlets[i])
		{
			pressure_BC[i] = (algorithm.Pin - algorithm.dP);
		}
		else
		{
			pressure_BC[i] = 0.0;
		}
	}
	//cout << "pressure_BC = \n" << pressure_BC << "\n\n";


	//2. complete laplacian ajaicency marix  
	averageDiagonal = methodsAlgorithm.laplacian_A(laplacian_A, network, algorithm);
	//cout << "average Diagonal = " << averageDiagonal << "\n";
	//cout << "laplacian_A = \n" << laplacian_A << "\n\n";

	//3. multiply laplacian by P_BC for non-BC entries of B
	dense_vector<double> A_BC(Np);
	A_BC = laplacian_A * pressure_BC;
	//cout << "A_BC = \n" << A_BC << "\n\n";

	//4. create A: laplacian for non-BC entries, zero for BC entries
	//diagonal filled in a separate pass
	for (int i = 0; i < Nt; i++)
	{
		head = network.throats[i][0];
		tail = network.throats[i][1];
		//check if non-BC for head and tail pores: 
		if (!algorithm.inlets[head] && !algorithm.outlets[head])
		{
			ins[head][head] << laplacian_A[head][head];
			//ins[head][tail] << laplacian_A[head][tail];
			//if (laplacian_A[head][head] < minA) { minA = laplacian_A[head][head]; }
			//if (laplacian_A[head][head] > maxA) { maxA = laplacian_A[head][head]; }
			//if (laplacian_A[head][head] > 1.0) { void __debugbreak(); }
		}
		if (!algorithm.inlets[tail] && !algorithm.outlets[tail])
		{
			ins[tail][tail] << laplacian_A[tail][tail];
			//ins[tail][head] << laplacian_A[tail][head];
			//if (laplacian_A[tail][tail] < minA) { minA = laplacian_A[tail][tail]; }
			//if (laplacian_A[tail][tail] > maxA) { maxA = laplacian_A[tail][tail]; }
			//if (laplacian_A[tail][tail] > 1.0) { void __debugbreak(); }
		}
		if ((!algorithm.inlets[head] && !algorithm.outlets[head]) && (!algorithm.inlets[tail] && !algorithm.outlets[tail]))
		{
			ins[head][tail] << laplacian_A[head][tail];
			ins[tail][head] << laplacian_A[tail][head];

			//if (laplacian_A[head][tail] < minA) { minA = laplacian_A[head][tail]; }
			//if (laplacian_A[head][tail] > maxA) { maxA = laplacian_A[head][tail]; }

			//if (laplacian_A[tail][head] < minA) { minA = laplacian_A[tail][head]; }
			//if (laplacian_A[tail][head] > maxA) { maxA = laplacian_A[tail][head]; }

			//if (laplacian_A[head][tail]  > 1.0) { void __debugbreak(); }
			//if (laplacian_A[tail][head]  > 1.0) { void __debugbreak(); }

		}
		//if (maxA > 1) { cout << i << "\n"; }
	}

	//4. cont-d: fill BC diagonal values
	//5. create b: average_diagonal*P_BC for BC entries, -laplacian*P_BC for non-BC entries
	for (int i = 0; i < Np; i++)
	{
		if (algorithm.inlets[i])
		{
			ins[i][i] << averageDiagonal;
			b[i] = averageDiagonal * algorithm.Pin;

			//if (b[i] < minb) { minb = b[i]; }
			//if (b[i] > maxb) { maxb = b[i]; }
		}
		else if (algorithm.outlets[i])
		{
			ins[i][i] << averageDiagonal;
			b[i] = averageDiagonal * (algorithm.Pin - algorithm.dP);

			//if (b[i] < minb) { minb = b[i]; }
			//if (b[i] > maxb) { maxb = b[i]; }
		}
		else
		{
			b[i] = -A_BC[i];

			//if (b[i] < minb) { minb = b[i]; }
			//if (b[i] > maxb) { maxb = b[i]; }
		}
	}
	//debug
	/*cout << "minA = " << minA << "\n";
	cout << "maxA = " << maxA << "\n";
	cout << "minb = " << minb << "\n";
	cout << "maxb = " << maxb << "\n";*/

	return b;
}

double classPNM::classMethods::classMethodsSinglePhaseStokesFlow::
		laplacian_A(compressed2D<double>& laplacian_A, 
					classPNM::classNetwork& network, 
					classPNM::classSinglePhaseStokesFlow& algorithm)
{
	laplacian_A = 0.0;
	int Np = network.pores.size();
	int Nt = network.throats.size();
	int head, tail; //pore indices for each throat
	double averageDiagonal = 0;
	// Create inserter for matrix A
	mat::inserter<compressed2D<double>, update_plus<double>> insert(laplacian_A);

	//complete A in laplacian form without BC 
	for (int i = 0; i < Nt; i++)
	{
		head = network.throats[i][0];
		tail = network.throats[i][1];

		insert[head][head] << algorithm.conductance[i];
		insert[head][tail] << -algorithm.conductance[i];

		insert[tail][tail] << algorithm.conductance[i];
		insert[tail][head] << -algorithm.conductance[i];

		averageDiagonal += 2 * algorithm.conductance[i];
	}
	averageDiagonal /= Np;
	return averageDiagonal;
}

void classPNM::classMethods::classMethodsSinglePhaseStokesFlow::
			solve_A_b(classPNM::classNetwork& network,
						classPNM::classSinglePhaseStokesFlow& algorithm)
{
	cout << "MTL BICGSTAB" << endl;
	//function for solving matirx equation once
	//local variables
	int Np = network.pores.size();
	compressed2D<double> * A = &algorithm.A; 
	dense_vector<double> * b = &algorithm.b; 
	dense_vector<double> X(Np, 0.0);
	//cout << "b = \n" << b << "\n\n";

	// Create an ILU(0) preconditioner
	pc::ilu_0<compressed2D<double>> P(*A);

	// Termination criterion: r < 1e-6 * b or N iterations
	cyclic_iteration<double>       iter(*b, 5000, 1.0e-8, 0, 100);

	// Solve Ax == b with left preconditioner P
	bicgstab(*A, X, *b, P, iter);

	
	//cout << "X = \n" << X << "\n\n";
	algorithm.X = X;
}

void classPNM::classMethods::classMethodsSinglePhaseStokesFlow::
		solve_A_b_pypardiso(classPNM::classNetwork& network,
							classPNM::classSinglePhaseStokesFlow& algorithm)
{
	cout << "PARDISO" << endl;
	//function for solving matirx equation once
	//local variables
	int Np = network.pores.size();
	compressed2D<double> * A = &algorithm.A;
	dense_vector<double> * bb = &algorithm.b;
	dense_vector<double> X(Np, 0.0);

	//
	//pypradiso example reworked
	//

	#if !defined(MKL_ILP64)
	#define IFORMAT "%i"
	#else
	#define IFORMAT "%lli"
	#endif

	vector<int> ia_vec, ja_vec;
	vector<double> b_vec;
	vector<double> x_vec(Np);
	vector<double> bs_vec(Np);

	//convert starts and indices to int vectors
	for (int i = 0; i < A->starts.size(); i++)
	{
		ia_vec.push_back((int)A->starts[i] + 1);
	}
	for (int i = 0; i < A->indices.size(); i++)
	{
		ja_vec.push_back((int)A->indices[i] + 1);
	}
	for (int i = 0; i < Np; i++)
	{
		b_vec.push_back(algorithm.b[i]);
	}

	MKL_INT n = Np;
	MKL_INT* ia = &ia_vec[0];
	MKL_INT* ja = &ja_vec[0];
	double* a = &A->data[0];

	/* RHS and solution vectors. */
	double* b = &b_vec[0];
	double* x = &x_vec[0];
	double* bs = &bs_vec[0];
	double res, res0;

	MKL_INT mtype = 11;       /* Real unsymmetric matrix */
  // Descriptor of main sparse matrix properties
	struct matrix_descr descrA;
	// Structure with sparse matrix stored in CSR format
	sparse_matrix_t       csrA;
	sparse_operation_t    transA;

	MKL_INT nrhs = 1;     /* Number of right hand sides. */
	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void *pt[64];
	/* Pardiso control parameters. */
	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error, msglvl;
	/* Auxiliary variables. */
	MKL_INT i, j;
	double ddum;          /* Double dummy */
	MKL_INT idum;         /* Integer dummy. */
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)
	{
		iparm[i] = 0;
	}
	iparm[0] = 1;         /* No solver default */
	iparm[1] = 2;         /* Fill-in reordering from METIS */
	iparm[3] = 0;         /* No iterative-direct algorithm */
	iparm[4] = 0;         /* No user fill-in reducing permutation */
	iparm[5] = 0;         /* Write solution into x */
	iparm[6] = 0;         /* Not in use */
	iparm[7] = 2;         /* Max numbers of iterative refinement steps */
	iparm[8] = 0;         /* Not in use */
	iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0;        /* Conjugate transposed/transpose solve */
	iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
	iparm[13] = 0;        /* Output: Number of perturbed pivots */
	iparm[14] = 0;        /* Not in use */
	iparm[15] = 0;        /* Not in use */
	iparm[16] = 0;        /* Not in use */
	iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;       /* Output: Mflops for LU factorization */
	iparm[19] = 0;        /* Output: Numbers of CG Iterations */
	maxfct = 1;           /* Maximum number of numerical factorizations. */
	mnum = 1;         /* Which factorization to use. */
	msglvl = 0;           /* Print statistical information  */
	error = 0;            /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)
	{
		pt[i] = 0;
	}
	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("ERROR during symbolic factorization: " IFORMAT, error);
		exit(1);
	}
	//printf("Reordering completed ... ");
	//printf("\nNumber of nonzeros in factors = " IFORMAT, iparm[17]);
	//printf("\nNumber of factorization MFLOPS = " IFORMAT, iparm[18]);
	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during numerical factorization: " IFORMAT, error);
		exit(2);
	}
	//printf("\nFactorization completed ... ");
	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;

	descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
	descrA.mode = SPARSE_FILL_MODE_UPPER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ONE, n, n, ia, ia + 1, ja, a);

	//solving  Ax=b
	
	iparm[11] = 0;        /* Conjugate transposed/transpose solve */
	transA = SPARSE_OPERATION_NON_TRANSPOSE;
		
	printf("\nSolving system  ...\n", iparm[11]);
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
	if (error != 0)
	{
		printf("\nERROR during solution: " IFORMAT, error);
		exit(3);
	}

	//put solution into MTL vector X
	for (j = 0; j < n; j++)
	{
		X[j] = x[j];
	}
	//cout << endl << X << endl;
	algorithm.X = X;
	
	// Compute residual
	mkl_sparse_d_mv(transA, 1.0, csrA, descrA, x, 0.0, bs);
	res = 0.0;
	res0 = 0.0;
	for (j = 1; j <= n; j++)
	{
		res += (bs[j - 1] - b[j - 1]) * (bs[j - 1] - b[j - 1]);
		res0 += b[j - 1] * b[j - 1];
	}
	res = sqrt(res) / sqrt(res0);
	//printf("\nRelative residual = %e", res);
	//cout << "\n\n";
	// Check residual
	if (res > 1e-10)
	{
		printf("Error: residual is too high!\n\n");
		exit(10 + i);
	}
	
	mkl_sparse_destroy(csrA);

	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1;           /* Release internal memory. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, &ddum, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
}

void classPNM::classMethods::classMethodsSinglePhaseStokesFlow::
			pass_result_to_phase(classPNM::classNetwork& network,
									classPNM::classMethane& phase,
									classPNM::classSinglePhaseStokesFlow& algorithm)
{
	int Np = network.pores.size();
	int Nt = network.throats.size();
	int head, tail;

	//pass solution vector X to phase pressure
	//check for pressure values out of bounds - due to unconnected regions
	for (int i = 0; i < Np; i++)
	{
		if (algorithm.X[i] > algorithm.Pin) { phase.pressure_pore[i] = algorithm.Pin; }
		else if (algorithm.X[i] < (algorithm.Pin - algorithm.dP)) { phase.pressure_pore[i] = (algorithm.Pin - algorithm.dP); }
		else
		{
			phase.pressure_pore[i] = algorithm.X[i];
		}
	}

	//interpolate pore pressure for throats
	for (int i = 0; i < Nt; i++)
	{
		head = network.throats[i][0];
		tail = network.throats[i][1];
		phase.pressure_throat[i] = 0.5*(phase.pressure_pore[head] + phase.pressure_pore[tail]);
	}
}

void classPNM::classMethods::classMethodsSinglePhaseStokesFlow::
			calculate_flow_rate(classPNM::classNetwork& network,
								classPNM::classMethane& phase,
								classPNM::classSinglePhaseStokesFlow& algorithm)
{
	int Np = network.pores.size();
	int Nt = network.throats.size();
	int head, tail;
	double Qt;

	algorithm.flow_rate_pore.clear();
	algorithm.flow_rate_throat.clear();
	vector<double> Qp(Np, 0);

	//calculate flow rate for every pore
	for (int i = 0; i < Nt; i++)
	{
		head = network.throats[i][0];
		tail = network.throats[i][1];
		Qt = algorithm.conductance[i] * (phase.pressure_pore[tail] - phase.pressure_pore[head]);

		Qp[head] -= Qt;
		Qp[tail] += Qt;
		algorithm.flow_rate_throat.push_back(Qt);
	}
	algorithm.flow_rate_pore = Qp;
}

void classPNM::classMethods::classMethodsSinglePhaseStokesFlow::
		calculate_effective_permeability(classPNM::classNetwork& network,
											classPNM::classSettings& settings,
											classPNM::classMethane& phase,
											classPNM::classSinglePhaseStokesFlow& algorithm)
{
	double iterator = 0;
	double domain_area = pow(settings.domain, 2);

	int Np = network.pores.size();
	//calculate Q*viscosity for each inlet pore
	for (int i = 0; i < Np; i++)
	{
		if (algorithm.inlets[i])
		{
			iterator += algorithm.flow_rate_pore[i] * phase.viscosity_pore[i];
		}
	}

	//calculate effective permeability in nano-Darcy
	double temp = iterator * settings.domain / (domain_area * algorithm.dP) / 9.86e-12 * 1e9;
	cout << "\nEffective Permeability: " << temp << " [nD]" << endl;
	algorithm.effective_permeability.push_back(temp);
}

void classPNM::classMethods::classMethodsSinglePhaseStokesFlow::
			run_once(classPNM::classNetwork& network,
						classPNM::classSettings& settings,
						classPNM::classMethane& phase,
						classPNM::classPhysics& physics,
						classPNM::classSinglePhaseStokesFlow& algorithm,
						classPNM::classMethods::classMethodsSinglePhaseStokesFlow& methodsAlgorithm)
{
	methodsAlgorithm.create_A_b(network, physics, algorithm, methodsAlgorithm);
	cout << "Algorythm Parameters Initialised" << endl;

	//methodsAlgorithm.solve_A_b(network, algorithm);
	methodsAlgorithm.solve_A_b_pypardiso(network, algorithm);
	methodsAlgorithm.pass_result_to_phase(network, phase, algorithm);

	methodsAlgorithm.calculate_flow_rate(network, phase, algorithm);
	methodsAlgorithm.calculate_effective_permeability(network, settings, phase, algorithm);
}

void classPNM::classMethods::classMethodsSinglePhaseStokesFlow::
		run_linear_iterative(classPNM::classNetwork& network,
								classPNM::classSettings& settings,
								classPNM::classMethane& phase,
								classPNM::classPhysics& physics,
								classPNM::classSinglePhaseStokesFlow& algorithm,
								classPNM::classMethods& methods)
{
	algorithm.effective_permeability.clear();
	int n = settings.n_linear_iterations; //iterations
	for (int i = 0; i < n; i++)
	{
		cout << "\n" << "Solver Iteration " << i << ":\n";
		methods.methodsPhase.generate_phase_properties(network, phase, methods.methodsPhase);
		methods.methodsPhysics.generate_physics_properties(network, phase, physics, methods.methodsPhysics);
		methods.methodsAlgorithm.run_once(network, settings, phase, physics, algorithm, methods.methodsAlgorithm);

	}
}