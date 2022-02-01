#pragma once
#define BOOST_TYPEOF_EMULATION
#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
// Windows Header Files
//#include <windows.h>
// C RunTime Header Files
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <tchar.h>

//#include <shobjidl.h> 
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <utility>
#include <stdexcept> 
#include <sstream> 
#include <iterator>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include "chrono"
#include <math.h>
#include <csignal>
//#include <Python.h>
#include <stdio.h>
#include <conio.h>
#include <iomanip>
#include <list>
#include<map>
#include <random>


using namespace std;
using namespace mtl;
using namespace itl;