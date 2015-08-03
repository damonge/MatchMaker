# MatchMaker - Friends-of-Friends halo finder.

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

Three output formats are provided, given by the parameter
"output_format" in the parameter file. These are :
 - ASCII : normal text file (ASCII).
 - FITS : FITS binary table.
 - BINARY : MatchMaker's own binary format.

ASCII and FITS will contain a 29 columns, corresponding to the
quantities:

 - 1  "ID"     - Unique Halo ID
 - 2  "NP"     - Number of particles
 - 3  "MASS"   - Mass (in units of M_sun/h)
 - 4  "PX_CM"  - Center of mass position (in units of Mpc/h)
 - 5  "PY_CM"
 - 6  "PZ_CM"
 - 7  "PX_RMS" - Standard deviation of the particle positions
 - 8  "PY_RMS"   (in units of Mpc/h)
 - 9  "PZ_RMS"
 - 10 "VX_CM"  - Center-of-mass velocity (in km/s)
 - 11 "VY_CM"
 - 12 "VZ_CM"
 - 13 "VX_RMS" - Velocity dispersion (in km/s)
 - 14 "VY_RMS"
 - 15 "VZ_RMS"
 - 16 "LX"     - Angular momentum (in units of Mpc/h * km/s)
 - 17 "LY"
 - 18 "LZ"
 - 19 "B"      - Ratio between the second and first (largest)
                  eigenvalues of the inertia tensor
 - 20 "C"      - Ratio between the third (smallest) and first
                  (largest) eigenvalues of the inertia tensor
 - 21 "EAX"    - Eigenvector for the first (largest)
 - 22 "EAY"      eigenvalue of the inertia tensor
 - 23 "EAZ"
 - 24 "EBX"    - Eigenvector for the second eigenvalue of the
 - 25 "EBY"      inertia tensor
 - 26 "EBZ"
 - 27 "ECX"    - Eigenvector for the third (smallest) eigenvalue
 - 28 "ECY"      of the inertia tensor
 - 29 "ECZ"

The halos are ordered by mass (largest to smallest). An example
showing how to interpret the MatchMaker binary format is given
in sample/test_read_binary.c. The directory "sample" also
contains a sample MatchMaker output file for each of the three
output formats (.txt -> ASCII, .dat -> BINARY, .fits -> FITS).


## License

MatchMaker is distributed under the GNU Public License v3
(see COPYING for details).


## Author

For questions or queries e-mail the author David Alonso:

   david dot alonso at astro dot ox dot ac dot uk
