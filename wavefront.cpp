// MPI imports
#include <mpi.h>
#include <mpi-ext.h>
#include "mapreduce.h"
#include "keyvalue.h"
using namespace MAPREDUCE_NS;

// C++ imports
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>


// C imports
#include <stdio.h>
#include <unistd.h>

// Voronoi imports
#include "itpsolv.h"
using namespace itp;

// Namespaced "globals" for optimization reasons
namespace wevoND
{

class Form
{
	private:
	std::vector<double> dat;
	// Const data for compilation efficiency
	const uint32_t dimension;
	const uint32_t manifold;
	public:
	Form() : dimension{0}, manifold{0}
	{
		dat.resize(dimension*manifold,0.0);
	}
	Form(uint32_t dim) : dimension{dim}, manifold{dim}
	{
		dat.resize(dimension*manifold,0.0);
	}
	Form(uint32_t dim, uint32_t man) : dimension{dim}, manifold{man}
	{
		dat.resize(dimension*manifold,0.0);
	}
	double& at(uint32_t i,uint32_t j) { return dat[dimension*i+j]; }
	const uint32_t space() { return dimension; }
	const uint32_t dim() { return manifold; }
};

class Vector
{
	private:
	std::vector<double> dat = std::vector<double>(dimension,0.0);
	const uint32_t dimension;
	public:
	Vector() : dimension{0} {}
	Vector(std::initializer_list<double> list) : dimension{uint32_t(list.size())}
	{
		size_t i = 0;
		for(auto&& x : list) { dat[i] = x; ++i; }
	}
	double& at(uint32_t i) { return dat[i]; }
	const uint32_t space() { return dimension; }
};

struct wave {
	// Main data, need these for a unique wave
	Vector center[3];
	Form form;
	double span[3];
	// Derivative, can be recalculated but are cached for speed
	bool isflat;
	double radius;
	// Const data for compile time easy access.
	// Yes these somewhat bloat MPI bandwidth but the ease of use is #worth
	const uint32_t dim = form.dim();
	const uint32_t space = form.space();
};

// Just a way to store the main runtime data, all the base waves their ids,
// and such
class wavefront
{
	const uint32_t dimension;
	const uint32_t manifoldDimension;
	public:
	wavefront(uint32_t dim, uint32_t man) :
		dimension{dim}, manifoldDimension{man}
	{
		
	}
	uint32_t dim() { return dimension; }
	void addBaseWave(Vector v)
	{
	}
};

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

struct wave inline initBaseWave(Vector center)
{
	wave ret = {{center,{},{}},Form{},{0.0,-1.0,-1.0}};
	ret.isflat = true;
	return ret;
}

void inline printWave(struct wave wave)
{
	Form testform{2};
	printf("Testform %d\n",testform.dim());
	printf("Wave print (%d,%d)...\n",wave.space,wave.dim);
	printf("%s\n", wave.isflat ? "infinite" : "ends");
	printf("Center:\n\t");
	for(int i = 0; i < 3; i++)
	{
		printf("(");
		for(int s = 0; s < wave.center[i].space(); s++)
		{
			printf("%f",wave.center[i].at(s));
		}
		printf(") ");
	}
	printf("\n");
	printf("Span:\n\t");
	for(int i = 0; i < 3; i++)
	{
		printf("(%f) ",wave.span[i]);
	}
	printf("\n");
	printf("FORM:\n");
	for(int d = 0; d < wave.dim; d++)
	{
		printf("\t(");
		for(uint32_t s = 0; s < wave.space; s++)
		{
			printf("%f ",wave.form.at(d,s));
		}
		printf(")\n");
	}
	printf("\n");
}

void testDataSending(int rank, MPI_Datatype** handle)
{

	// TODO: actually test the send
	wave basewave = initBaseWave(Vector{0.0,1.0,2.0,3.0});

	printWave(basewave);

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
void bindMPIDatatype(MPI_Datatype** handle, uint32_t space, uint32_t manifold)
{
	// The solver bound datatype
	MPI_Datatype solvbase[4] = {
		MPI_CXX_BOOL,
		MPI_CXX_BOOL,
		MPI_DOUBLE,
		MPI_DOUBLE
	};
	MPI_Aint solvoffset[4] = {
		offsetof(SOLVEBOUND,verified),
		offsetof(SOLVEBOUND,safe),
		offsetof(SOLVEBOUND,lower),
		offsetof(SOLVEBOUND,upper)
	};
	int solvcount[4] = {1,1,1,1};
	MPI_Type_create_struct(4,solvcount,solvoffset,solvbase,handle[0]);
	MPI_Type_commit(handle[0]);
	// Datatype for a vector, single block of space size
	MPI_Type_vector(1,space,0,MPI_FLOAT,handle[1]);
	MPI_Type_commit(handle[1]);
	// Datatype for a geometric form, square block
	MPI_Type_vector(manifold,space,0,MPI_FLOAT,handle[2]);
	MPI_Type_commit(handle[2]);
	// The waveform datatype
	MPI_Datatype wavebase[3] = {
		*handle[1],
		*handle[2],
		MPI_DOUBLE
	};
	MPI_Aint waveoffset[3] = {
		offsetof(wave,center),
		offsetof(wave,form),
		offsetof(wave,span)
	};
	int wavecount[3] = {3,1,3};
	MPI_Type_create_struct(3,wavecount,waveoffset,wavebase,handle[3]);
	MPI_Type_commit(handle[3]);
}

}

wevoND::Vector readInputLine(std::string instring, uint32_t dim)
{
	std::istringstream linestream(instring);
	return wevoND::Vector{0.0,0.0,0.0};
}

void readInputFile(wevoND::wavefront mesh, std::ifstream infile)
{
	std::string line;
	while(std::getline(infile,line))
	{
		mesh.addBaseWave(readInputLine(line,mesh.dim()));
	}
}

int main(int argc, char **argv)
{
	// Initialize the MPI runtime
	MPI_Init(&argc,&argv);
	int initialized;
	MPI_Initialized(&initialized);
	if(!initialized) { return 1; }
	int rank,size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// Find the data file input
	uint32_t dimension = -1;
	uint32_t manifoldDimension = -1;
	std::string infile = "lattice.dat";
	std::string outfile = "partition.dat";
	int arg = 0;
	// optarg is a global defined in unistd.h, so use requires annotation
	while((arg = getopt(argc, argv, "iodm:")) != -1)
	{
		switch(arg)
		{
			case 'i':
				infile = optarg;
			break;
			case 'o':
				outfile = optarg;
			break;
			case 'd':
				dimension = atoi(optarg);
			break;
			case 'm':
				manifoldDimension = atoi(optarg);
			break;
			default:
			continue;
		}
	}
	// Find and broadcast the dimensionality(ies) we are working in
	if(rank == 0)
	{
		std::ifstream indat;
		indat.open(infile);
		std::string instr;
		std::getline(indat, instr);
		std::istringstream linestream(instr);
		indat.close();
		double var;
		uint32_t newdim = 0;
		while(linestream >> var) { newdim++; }
		if(dimension == -1) { dimension = newdim; }
		if(manifoldDimension == -1) {manifoldDimension = dimension; }
	}
	
	MPI_Bcast(&dimension,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&manifoldDimension,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Dimensions are: %d,%d\n",dimension,manifoldDimension);
	MPI_Finalize();
	return 0;

	MPI_Datatype solvBoundMsg;
	MPI_Datatype vecMsg;
	MPI_Datatype formMsg;
	MPI_Datatype waveMsg;
	MPI_Datatype* handle[4] = {
		&solvBoundMsg,
		&vecMsg,
		&formMsg,
		&waveMsg
	};
	wevoND::bindMPIDatatype(handle, dimension, manifoldDimension);
	
	MPI_Finalize();
	return 0;
	// Create and broadcast the wavefront
	// TODO: multicore file read (needs binary input)
	//auto wavefront = new wevoND::wavefront(dimension,manifoldDimension);
	// Read and construct all the base waves from the input file
	// Execute our ITP code
	MapReduce *mr = new MapReduce(MPI_COMM_WORLD);
	uint64_t kvupdate = mr->map(4,&wevoND::testITPsolv,NULL);
	delete mr;
	// Check sending all our structures
	wevoND::testDataSending(rank, handle);
	// Finalize everything
	MPI_Finalize();
	return 0;
}
