# Stochastic parameterization in CROCO

This directory contains modules to introduce stochastic parameterizations in CROCO.
This includes the modules to generate stochatic processes with the required statistical properties,
the interface with CROCO, and the stochastic parameterizations themselves.

### Stochastic modules

These modules are meant to generate maps of stochastic processes
with the requested statistical properties
(time and space correlation, marginal distribution,...).
They are fully independent of CROCO and can follow updates
of the stogen package without interference with anything else.

- **_stoarray_** :
    This is the data module, where all stochastic fields are defined and stored:
    - it receives the requests from the users
      (with routines `sto_array_request_size` and `sto_array_request_new`),
    - it allocates the required arrays according to requests,
      and check the consistency of the requested features
      (with routine `sto_array_init`).

- **_stopar_** :
    This is the time evolution module, where all stochastic fields are updated at each time step:
    - initialization phase (routine `sto_par_init`):
      - seed random number generator (according to subdomain and member indices),,
      - initialize methods to generate spatially correlated random fields
        (calls to `sto_diff_init`, `sto_kernel_init`,...),
      - initialize transformations to requested marginal distributions
        (call to `sto_marginal_init`),
      - initialize parameters of autoregressive processes,
      - initialize random fields (from restart or from the requested method
        to generate spatially correlated random fields: `sto_diff`, `sto_kernel`,...);
    - time update (routine `sto_par`):
      - forward the autoregressive process in time  (or
        interpolate between a past and future state of the autoregressive process),
      - perform the transformation to the requested marginal distribution.

- **_stowhite_** :
    Generate a map of Gaussian white noise, with zero mean and unit standard deviation.

- **_stodiff_** :
    Generate a map of spatiallye correlated noise, with zero mean and unit standard deviation,
    using the recursive application of a diffusion operator on a white noise.

- **_stokernel_** :
    Generate a map of spatiallye correlated noise, with zero mean and unit standard deviation,
    using the convolution of a white noise with a filtering kernel.

    The convolution integral is computed using a Quasi Monte Carlo approximation,
    by summing over a limited number of kernel locations.

    The Quasi Monte Carlo sequence of kernel locations
    is obtained from a 2D random Sobol sequence (with module `stosobolseq`).

    Options for the filtering kernel include: Gaussian kernel, Laplacian kernel,
    box kernel, triangle kernel, Mexican hat wavelet (Ricker wavelet),
    Morlet wavelet a (with specific choice of frequency, adjust if needed).

    Options for computing distances include: grid coordinates, Cartesian coordinates,
    spherical coordinates (more expensive).

- **_stosobolseq_** :
    Module to generate mutlidimensional Sobol sequences
    (obtained from https://github.com/DaanVanVugt).

- **_stomarginal_** :
    Transform the Gaussian process to the requested marginal distribution.

- **_storng_kiss_** :
    Random number generator. This includes the kiss32 and kiss64 random number generators
    and the polar method to transform the integer sequence into Gaussian numbers.

- **_storng_ziggurat_** :
    Random number generator. This includes the shr3 random number generator
    and the ziggurat method to transform the integer sequence into Gaussian numbers.

### Interface with CROCO

- **_stomod_** :
    Main stochastic module embedding all stochastic parameterizations.
    This module includes the two main routines that are called from within CROCO:
    - initialization phase (routine `sto_mod_init`):
      - initialization of every dynamical stochastic parameterizations (e.g. `sto_bulk_init`),
      - initialization of the structure of the stochastic arrays (call to `sto_array_init`),
      - initialization of the time iteration of the stochastic arrays (call to `sto_par_init`);
    - time update (routine `sto_mod`):
      - update stochastic fields (call to `sto_par`),
      - apply dynamical stochastic parameterization (e.g., call to `sto_bulk`).

    The call of the routines may need to be organized differently depending on
    where the stochastic parameterization code must be used CROCO.

- **_stoexternal_** :
    This module is used by the stochastic code to collect all information
    required from CROCO: type of variables, description of the model grid,
    ensemble parameters, lateral boundary conditions (or connection between subdomains).
    This is the only place where model data go to the stochastic modules,
    so that this can be easily identified and possibly upgraded. This is model dependent.

### Stochastic parameterizations in CROCO

There must an additional module for each type of parameterization introduced in CROCO.
At this stage, we have introduced only one:

- **_stobulk_** :
    Stochastic parameterization of the bulk formulation for the air-sea fluxes.
    At this stage, we have only introduced simple perturbations of the drag coefficient.

### Modification in the CROCO code (in the OCEAN directory)

These modifications are all included with the **CPP key STOGEN**.
Deactivating this key should remove all stochastic features.
They include:

- The definition of the CPP key STOGEN in `cppdefs.h`  

- The compilation STOGEN scripts in `jobcomp`

- The initialization of the stochastic modules in `main.F`.

- The update of the stochastic processes in `step.F`.

- The use of the stochastic arrays if and where needed.

This last kind of modification should be avoided as much as possible.
The effect of the stochastic arrays in CROCO should be embedded
in each specific stochastic parameteriation module (e.g. stobulk) whenever possible.
In any case, each stochastic parameterization should never represents
more than a few lines of code (e.g. a call to a routine) in the main CROCO code.

### What remains to be done

- There can be mistakes in the size of the arrays,
  and the exchange of information between subdomains?
  This really needs to be checked.
  
- Grid an mask are made to point to CROCO arrays in memory (does it work?).
  At this stage, the global grid arrays are left unallocated in stoexternal.
  This means that all options that require them
  (i.e. stokernel in world coordinates) will just fail.

- Read parameters in stobulk from namelist or parameter file.

- Introduce restart files for the stochastic arrays.

- Do we need to deal with the OpenMP parallelization of CROCO 
  or are we happy with MPI paralelization only ?
  If we need OpenMP, this would require modifying the range of the variations
  of the i and j indices in the whole stochastic code (as it is done for PISCES).
  This would also require a testcase with OpenMP.
