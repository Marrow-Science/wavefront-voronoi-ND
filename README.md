# wavefront-voronoi-ND

This is an extension of work on wavefront voronoi in 2D space as exemplified in
wevo (https://github.com/cgalab/wevo.git).
By abstracting wave intersections,
the wavefront algorithm can be extended to ND space for use in high-dimensional
applications such as materials science and machine learning.

# Dependencies

The basic algorithm runs on ```scipy```.
The visualization code ```vis.py``` depends on the ```Panda3d``` library
for 3d rendering (https://www.panda3d.org).
Install these dependencies using ```pip```:

```pip install panda3d```
```pip install scipy```

# Running

```wevoNd.py``` contains a simple command-line runtime for the ```voronoi.py```
library functions.
Execute by using python3.

```python3 wevoNd.py```

# Command line arguments

-o, --output	<filename> designate an output file

-i, --input	<filename> designate an input file

-b, --binary	use binary input/output encodings

-v, --vis	run a panda3d visualization of the algorithm

-d, --debug	print debug output

-h, --help	display a help message

Without designated input and output files, inputs are read from stdin and
outputs are printed to stdout.

# Input file format

Human-readable input consists of lines.
Each line contains a space-seperated list of floats designating a point
location, with an optional weight designated by a "w=".
Points without a weight are assumed to have a weight of 1.0.
Comments after a \# on a line are ignored.
IDs may be manually set using an exclamation point (!).
For example:

```
1.0 2.0 3.0 w=4.0
3.0 2.0 1.0 w=2.0 #3D point at 3.0 2.0 1.0 with weight 2.0
!Foo 2.0 3.0 1.0 #3D point at 2.0 3.0 1.0 with weight 1.0 and the ID "Foo"
```

# Output file format

Human-readable output consists of a list of parent points with IDs
followed by a list of verticies.
Each parent consists of an ID followed by a point and then a weight.
For example:

```
!1 1.0 2.0 3.0 w=4.0
!2 2.0 3.0 4.0 w=2.0
```

Each vertex contains an ID, a location in space, a list of parent nodes,
and an arc center plus form subspace for computing full vertex subsurfaces.
IDs start with an exclamation point (!).
Locations are space-separated.
Parents are space-separated within brackets ([]).
Arc center follows a "c=" and is contained in parentheses (()).
Forms follow a "f=" and are a space-separated 2-d lists of vectors enclosed
in brackets ([[][]]).
For example:

```
!12 2.0 1.0 3.0 [1 2] c=(2.0 3.0 1.0) f=[[1.0 2.0 3.0] [3.0 2.0 1.0]]
```

# Binary file format NOT YET IMPLEMENTED

The library functions "readbinI" and "readbinO" reads binary input and output
files into internal ```wavefront``` and ```partition``` objects.

Binary formats are not yet implemented.

# MPI Code NOT YET IMPLEMENTED

The MPI Code is currently being developed and tested.
It depends on the OpenMPI c++ compiler
and the Sandia National Laboratories Map Reduce MPI ```MR-MPI``` code.
There is a provided make script for checking and installing these dependencies.

```make configure```

Compile the library object files using make.
By default they will be placed in the build/ directory

```make```

This package provides a few examples for using the mpi code.
They depend on the wevoNd library files. Compile them using Make.

```make examples```

It is recommened to run the included mpi examples using ```slurm```

```sudo apt-get install slurmctld```
```sudo apt-get install slurmd```

Slurm needs detailed setup, see https://slurm.schedmd.com/quickstart_admin.html

After the detailed setup, run the examples using slurm.

```srun voronoi-test```

Or use the provided MPI runtime directly.
The MPI runtime is identical to the python runtime with a few restrictions.
It must read from files, not stdin, and must write to files, not stdout.

```srun wevoNd -i <input_file> -o <output_file>```

# WARNING

The algorithm implementation is not yet thouroughly tested, although it does produce correct output
in many cases.

# Future work

Stay tuned for an MPI implementation, more testing, and further algorithm optimizations.

# Love you!
I hope you like it!
