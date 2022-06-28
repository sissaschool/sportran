# Example inputs files for LAMMPS simulation of silica glass

In the folder [examples/data/Silica/lammps](https://github.com/sissaschool/sportran/tree/develop/examples/data/Silica/lammps) you can find an example LAMMPS input file, to generate heat-flux data that can be analyzed with SporTran.
To run it you need a working LAMMPS installation and to launch the simulation with, for example:

	mpirun -np 12 lmp_mpi -in silica.in -log silica.out

After this it is possible to extract the data and export it into a simple *table file* and a *numpy binary data file* by running the script `convert_lammps_log.sh`.
This input generated the data used in the silica examples.
