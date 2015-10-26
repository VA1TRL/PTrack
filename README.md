# Off-Line Lagrangian Particle Tracking

## About

The particle tracking model this program implements (PTrack) was originally developed as a module embedded in the Unstructured Grid Finite Volume Community Ocean Model (FVCOM), first developed by Dr. Changsheng Chen et al. at the University of Massachusetts - Dartmouth. The current iteration of the particle tracking code found here is an almost complete rewrite of this module, designed to run independently of FVCOM. Presently FVCOM is being developed and maintained by the UMASSD Marine Ecosystem Dynamics Modeling Laboratory (MEDML), and more information may be found on their website: http://fvcom.smast.umassd.edu/

## Running the model

PTrack requires the creation of three different configuration and data files.

* Model configuration file (text)
* Particle seed file (text)
* Tidal forcing and model domain file (NetCDF)

The tidal forcing file contains the flow-field data to drives the movement of the particles being simulated, as well as the data needed to reconstruct the physical domain (mesh, bathymetry, etc.), and is produced as the primary output of FVCOM. The structure and content of the text files will be elaborated on in the following sections. After compiling the input, the model may be run from the command-line using the following syntax:

    $> .\PTrack run.dat

Where *$>* represents the command prompt and *run.dat* is the model configuration file.


### Configuration file

PARAMETER | DESCRIPTION
--------- | -----------------------------------------------
DTI       | Internal simulation time step (seconds)
NDRFT     | Number of particles in tracking simulation
INSTP     | NetCDF input file time step (seconds)
DTOUT     | Output interval (seconds)
DAYST     | Delay before particle tracking begins (relative to beginning of NetCDF file)
P_SIGMA   | Run vertical simulation over sigma levels, instead of meters
F_DEPTH   | Run simulation holding particle depth constant
P_REL_B   | Particle positions relative to the bottom (instead of surface)
OUT_SIGMA | Output particle z position as sigma depth instead of meters
GRIDFN    | NetCDF input file name, containing grid and flow field data
OUTFN     | Output file name
STARTSEED | Particle location input filename


### Particle seed file

(By column, each row defines a different particle)

 COLUMN         | TYPE | DESCRIPTION
--------------- | ---- | ----------------------------------
 ID             | INT  | Arbitrary identifier
 X              | REAL | Domain co-ordinates (meters)
 Y              | REAL | Domain co-ordinates (meters)
 Z              | REAL | Particle depth (meters)
 Release Time   | REAL | Time to release particle into the simulation (mjd)
 Track End Time | REAL | Time to remove particle from the simulation (mjd)


### Output

(By column, each row defines a different particle)

 PARAMETER      | TYPE | DESCRIPTION
--------------- | ---- | ----------------------------------
 ID             | INT  | Arbitrary identifier
 X              | REAL | Domain co-ordinates (meters)
 Y              | REAL | Domain co-ordinates (meters)
 Z              | REAL | Particle depth (meters)
 TIME           | REAL | Time of particle record (mjd)

