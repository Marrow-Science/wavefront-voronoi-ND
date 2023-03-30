// MPI imports
#include <mpi.h>
#include <mpi-ext.h>
#include "mapreduce.h"
#include "keyvalue.h"
using namespace MAPREDUCE_NS;

// C++ imports
#include <cmath>
#include <iostream>


// C imports
#include <stdio.h>

// Voronoi imports
#include "itpsolv.h"
using namespace itp;


// TODO: remove this code
double linfunctest(double x) { return 2*x - 12; }
double quadfunctest(double x) { return x*x - 6*x + 10; }
double cubfunctest(double x) { return x*x*x + 2*x*x; }

void printmap(int itask, KeyValue* kv, void* ptr)
{
	printf("Hello world: %i\n",itask);
}

void testITPsolv(int itask, KeyValue* kv, void* ptr)
{
	double res = 10.0;
	double eps = 0.00000001;
	switch(itask)
	{
		case 0:
		res = ITPsolvSafe(linfunctest,{false,false,-20.0,40.0},eps);
		break;
		case 1:
		res = ITPsolvSafe(quadfunctest,{false,false,-20.0,20.0},eps);
		break;
		case 2:
		res = ITPsolvSafe(cubfunctest,{false,false,-20.0,30.0},eps);
		break;
		default: break;
	}
	printf("Solution: %f\n",res);
}

// Each wave stores all data needed to calculate wave intersections.
//class wave
//{


//}

// Initialize all the MPI datatypes, for use in messages
void bindMPIDatatype(MPI_Datatype** handle)
{
	// The solver bound datatype
	MPI_Datatype basetype[1];
	MPI_Aint offset[1];
	int count[1];
	basetype[0] = MPI_DOUBLE;
	count[0] = 2;
	offset[0] = 0;
	MPI_Type_create_struct(1,count,offset,basetype,handle[0]);
	MPI_Type_commit(handle[0]);
	// The waveform datatype
	// TODO!
}

int main(int argc, char **argv)
{
	try {
	// Initialize the MPI runtime
	MPI_Init(&argc,&argv);
	int initialized;
	MPI_Initialized(&initialized);
	if(!initialized) { return 1; }
	// Bind all the datatypes we need to the appropriate handlers
	MPI_Datatype solvBound;
	MPI_Datatype* handle[1] = {&solvBound};
	bindMPIDatatype(handle);
	// Execute our ITP code
	MapReduce *mr = new MapReduce(MPI_COMM_WORLD);
		uint64_t kvupdate = mr->map(4,&testITPsolv,NULL);
	delete mr;
	MPI_Finalize();
	} catch(itp::BoundError e)
	{
		printf("Caught error %s\n",e.what());
	}
	return 0;
}
