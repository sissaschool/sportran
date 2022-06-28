# Example inputs for LAMMPS simulation of glass silica

In the folder [examples/data/Silica/lammps](https://github.com/sissaschool/sportran/tree/develop/examples/data/Silica/lammps) you can find an example input file to generate some data that you can analyze with sportran.
To run it you need a working LAMMPS installation and to run the file with, for example

	mpirun -np 12 lmp_mpi -in silica.in -log silica.out

After this it is possible to extract a simple table file by running the script `convert_lammps_log.sh`.
This input generated the data used in the silica examples.
