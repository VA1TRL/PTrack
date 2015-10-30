# PTrack
Off-Line Lagrangian Particle Tracking

# About
The particle tracking model this program implements (PTrack) was originally developed as a module embedded in the Unstructured Grid Finite Volume Community Ocean Model (FVCOM), first developed by Dr. Changsheng Chen et al. at the University of Massachusetts - Dartmouth. The current iteration of the particle tracking code found here is an almost complete rewrite of this module, designed to run independently of FVCOM. Presently FVCOM is being developed and maintained by the UMASSD Marine Ecosystem Dynamics Modelling Laboratory (MEDML), and more information may be found on their website: http://fvcom.smast.umassd.edu/

# Running the model
PTrack requires the creation of three different configuration and data files.

* Model configuration file (text)
* Particle seed file (text)
* Tidal forcing and model domain file (NetCDF)

The tidal forcing file contains the flow-field data to drives the movement of the particles being simulated, as well as the data needed to reconstruct the physical domain (mesh, bathymetry, etc.), and is produced as the primary output of FVCOM. The structure and content of the text files will be elaborated on in the following sections. After compiling the input, the model may be run from the command-line using the following syntax:
```sh
.\PTrack run.dat
```
Where `run.dat` is the model configuration file.

## Configuration file
 Parameter    | Units     | Description
 ------------ | --------- | -----------
 `DTI`        | seconds   | Simulation time step
 `DTOUT`      | seconds   | Output interval
 `F_DEPTH`    | boolean   | Holding particle depth constant
 `P_REL_B`    | boolean   | Particle positions are relative to the bottom (instead of surface)
 `OUT_SIGMA`  | boolean   | Output particle z position as sigma depth instead of meters
 `GRIDFN`     | file path | NetCDF grid and flow-field file
 `OUTFN`      | file path | Output file
 `STARTSEED`  | file path | Particle initial position (seed) file
 `P_RND_WALK` | boolean   | Apply random walk behaviour to the active particles
 `K_XY`       | m^2/s     | Horizontal particle diffusivity (if `P_RND_WALK` is **T**)
 `K_Z`        | m^2/s     | Vertical particle diffusivity (if `P_RND_WALK` is **T**)

## Particle seed file
(By column, each row defines a different particle)

 Column | Name           | Type | Units  | Description
 ------ | -------------- | ---- | ------ | -----------
 1      | ID             | INT  | *NA*   | Arbitrary identifier
 2      | X              | REAL | meters | Domain co-ordinates
 3      | Y              | REAL | meters | Domain co-ordinates
 4      | Z              | REAL | meters | Particle depth
 5      | Release Time   | REAL | MJD    | Time to release particle into the simulation
 6      | Track End Time | REAL | MJD    | Time to remove particle from the simulation

## Output
(By column, each row defines a different particle)

 Column | Name | Type | Units  | Description
 ------ | ---- | ---- | ------ | -----------
 1      | ID   | INT  | *NA*   | Arbitrary identifier
 2      | X    | REAL | meters | Domain co-ordinates
 3      | Y    | REAL | meters | Domain co-ordinates
 4      | Z    | REAL | meters | Particle depth
 5      | Time | REAL | MJD    | Time of particle record

