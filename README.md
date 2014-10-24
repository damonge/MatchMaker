# MatchMaker

MatchMaker is a simple MPI-parallel FoF halo finder.


## Algorithm

The algorithm is as follows:

1. Read particles and separate amongst nodes according
to their position along the x-axis in equi-spaced
intervals.
2. Communicate particles in a small buffer at the
beginning of each node's interval to the node on
the left.
3. Each node runs the FoF starting only with particles
within its own interval, but allowing them to be
friends with particles from the right node's buffer.
4. Resolve redundant halos. A halo found by node A with
particles in the buffer of its right node (A+1) has
precedence over the one found by the A+1
(because it also contains particles from the node's
interval). For this algorithm to succeed the size
of the buffer has to be larger than the size of the
largest halo in the simulation. For normal runs
a buffer of 5-10 Mpc/h should therefore be safe.
5. Compute halo properties: mass, center of mass,
size, velocity, rms velocity.


## Compiling

MatchMaker doesn't depend on any external libraries
(besides MPI), just type "make". Look at the Makefile in
order to tweak compilation flags.


## Parameter file

The behavior of MatchMaker is controlled by a param file
that is fed to the code as a single command-line argument.
A sample param file is provided (param_sample.ini). The
comments within it should be sufficient to understand the
role of each parameter.


## Running

Once compiled just type:

> mpirun -np number-of-processors ./MatchMaker param-file


## Input and output format

In its current version MatchMaker only reads snapshot files
in GADGET-1 format, and assumes periodic boundary conditions.
MatchMaker reads off a lot of information from the snapshot's
header, so make sure the contents of this header are correct!
Furthermore, the code currently assumes that all positions
are given in units of Mpc/h, velocities in km/s and masses in
10^10 M_sun/h. The results are also provided in these units.

The only supported output format is ASCII files with a number
of columns for different halo properties. The first line in
the output file explains the contents of each column.


## License

MatchMaker is distributed under the GNU Public License v3
(see COPYING for details).


## Author

For questions or queries e-mail the author David Alonso:

   david dot alonso at astro dot ox dot ac dot uk
