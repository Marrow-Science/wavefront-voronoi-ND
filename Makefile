INCLUDE=$(realpath lib/MR-MPI/)
SOURCE=$(wildcard $(realpath obj/Obj_linux/)/*.o)

all:
	mpicxx.openmpi wavefront.cpp $(SOURCE) -o voronoi -I $(INCLUDE)

run:
	srun ./voronoi
