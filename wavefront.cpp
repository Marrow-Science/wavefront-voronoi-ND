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
double quadfunctest(double x) { return x*x - 8; }
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
		res = ITPsolvSafe(quadfunctest,{false,false,-0.0,20.0},eps);
		break;
		case 2:
		res = ITPsolvSafe(cubfunctest,{false,false,-20.0,30.0},eps);
		break;
		default: break;
	}
	printf("Solution: %f\n",res);
}

void testDataSending(int rank, MPI_Datatype** handle)
{
	printf("Rank is %d\n",rank);
	SOLVEBOUND msg = {false, false, -10.0*(rank+2), 20.0*(rank+2)};
	MPI_Status status;
	if(rank == 0)
	{
		MPI_Bcast(&msg, 1, *handle[0], 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		printf("%d: Sent message %d %d, %f %f\n",rank,msg.verified, msg.safe,msg.lower,msg.upper);
	}
	if(rank != 0)
	{
		printf("%d: Before msg is %d %d %f %f\n",rank,msg.verified,msg.safe,msg.lower,msg.upper);

		MPI_Bcast(&msg, 1, *handle[0], 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		printf("%d: After msg is %d %d %f %f\n",rank,msg.verified,msg.safe,msg.lower,msg.upper);
	}
}

// Initialize all the MPI datatypes, for use in messages
void bindMPIDatatype(MPI_Datatype** handle)
{
	// The solver bound datatype
	MPI_Datatype basetype[4] = {
		MPI_CXX_BOOL,
		MPI_CXX_BOOL,
		MPI_DOUBLE,
		MPI_DOUBLE
	};
	MPI_Aint offset[4] = {
		offsetof(SOLVEBOUND,verified),
		offsetof(SOLVEBOUND,safe),
		offsetof(SOLVEBOUND,lower),
		offsetof(SOLVEBOUND,upper)
	};
	int count[4] = {1,1,1,1};
	MPI_Type_create_struct(4,count,offset,basetype,handle[0]);
	MPI_Type_commit(handle[0]);
	// Datatype for a geometric form
	MPI_Datatype formvec;
	
	// The waveform datatype
	// TODO!
}

int main(int argc, char **argv)
{
	// Initialize the MPI runtime
	MPI_Init(&argc,&argv);
	int initialized;
	MPI_Initialized(&initialized);
	if(!initialized) { return 1; }
	// Bind all the datatypes we need to the appropriate handlers
	int rank,size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Datatype solvBound;
	MPI_Datatype* handle[1] = {&solvBound};
	bindMPIDatatype(handle);
	// Execute our ITP code
	MapReduce *mr = new MapReduce(MPI_COMM_WORLD);
	uint64_t kvupdate = mr->map(4,&testITPsolv,NULL);
	delete mr;
	// Check sending all our structures
	testDataSending(rank, handle);
	// Finalize everything
	MPI_Finalize();
	return 0;
}
