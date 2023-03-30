INCLUDE=$(realpath lib/MR-MPI/)
SOURCE=$(wildcard $(realpath obj/Obj_linux/)/*.o)
OBJ=itpsolv.o
CC=mpicxx.openmpi 


itpsolv:
	$(CC) -c itpsolv.cpp

all: itpsolv
	$(CC) wavefront.cpp $(SOURCE) $(OBJ) -o voronoi -I $(INCLUDE)

run:
	srun ./voronoi
