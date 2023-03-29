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

// Use the Interpolate, Truncate, Project algorithm to solve the given function
// on the given interval. Is not safe, will not do checks

typedef struct {
	double lower, upper, eps;
} SOLVEBOUND;

const double k1 = 0.1;
const double k2 = 2.56567330897; // (0.98*(phi+1)), adjusted as per the paper

double ITPsolv(double(*function)(double),SOLVEBOUND inbound, double eps)
{
	double window = inbound.upper - inbound.lower;
	double nbisection = std::log2(window/(2.0*eps));
	double nmax = nbisection + 1.0; // n0 is one, based on prev work
	int iteration = 0;
	SOLVEBOUND bound = inbound;
	while(bound.upper - bound.lower > 2*eps)
	{
		// Interpolation: bisection and regula falsi points
		double fa = function(bound.lower);
		double fb = function(bound.upper);
		double bisection = (bound.lower + bound.upper) / 2;
		double regulafalsi = (bound.upper*fa - bound.lower*fb)/(fa-fb);
		// Truncation, perturb towards center
		bool alph = std::signbit(bisection - regulafalsi);
		double d1 = k1*std::pow((bound.upper - bound.lower),k2);
		double d2 = std::abs(bisection - regulafalsi);
		double del = d1 < d2 ? d1 : d2;
		double xt = alph ? regulafalsi + del : regulafalsi - del;
		// Projection, project the estimator to new minmax interval
		double win = (bound.upper - bound.lower);
		double pr = eps*std::pow(2.0,nmax - iteration)-(win/2.0);
		double pf = std::abs(xt - bisection);
		double proj = pr < pf ? pr : pf;
		// Find the new bound!
		double xITP = alph ? bisection - proj : bisection + proj;
		double yITP = function(xITP);
		if(yITP > 0) { bound.upper = xITP; }
		else if(yITP < 0) { bound.lower = xITP; }
		else { bound.lower = xITP; bound.upper = xITP; }
		iteration++;
	}
	return (bound.upper + bound.lower) / 2.0;
}

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
		res = ITPsolv(linfunctest,{-20.0,40.0},eps);
		break;
		case 1:
		res = ITPsolv(quadfunctest,{-20.0,20.0},eps);
		break;
		case 2:
		res = ITPsolv(cubfunctest,{-20.0,30.0},eps);
		break;
		default: break;
	}
	printf("Solution: %f\n",res);
}

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
	return 0;
}
