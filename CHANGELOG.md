# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [Unreleased] - xxxx-xx-xx
### Added

- ABL1D : 
  - TODO add description
  - Add KILPATRICK test-case

- SUBMASSBALANCE : Add submassbalance feature as in MARS, see issue 
  [#65](https://gitlab.inria.fr/croco-ocean/croco/-/issues/65) ; 
  check merge request and documentation if interested

- PISCES :
  - TODO add description
  - Update PISCES (XIOS on-the-fly, ...)

- SED_DENS : Add suspended sediment contribution to density whith cpp key 
  #SED_DENS. #SED_DENS applies to #SEDIMENT and #MUSTANG cases. For MUSTANG 
  case, variable using #key_sand2D are exclude of the contribution. See 
  issue [#124](https://gitlab.inria.fr/croco-ocean/croco/-/issues/124)

- CPL : management of WCHPC (fix issue #149)

- OBSTRUCTION : Add a process-based model for 3-dimensional simulation of 
  flow in presence of various obstructions, activated with cpp key #OBSTRUCTION. 
  Obstructions can be rigid or flexible,  and of 3 types : 
  - upward (like seagrass), 
  - downward (like mussel long-line),
  - 3D (like oyster tables)
  This module requires the keys #SOLVE3D, #GLS_MIXING and #GLS_KEPSILON. See issue
  [#123](https://gitlab.inria.fr/croco-ocean/croco/-/issues/123); 
  check merge request and documentation if interested.

- GPU : [#142]
  - **MAJOR :   Python is now required to compile croco. Even if no GPU.**
    - A second stage of pre-processing is added before compilation
    - Dependency :  os, re, sys, argparse, tempfile, shutil
  - OPENACC : Active GPU computation by adding OPENACC cpp keys in cppdefs.h
  - M3FAST_HIS, MP_M3FAST_SEDLAYERS : Keys for acoustic sediment layers
  - OPENACCNUMERIC : help on numeric validation cpu vs gpu 
    with RVTK_DEBUG (work in progress)
  - Limitations :
    - Test with Nvidia gpu and nvfortan (pgi also ok)
    - No GPU with : AGRIF MUSTANG OBSTRUCTION   PISCES modules
    - BAND_DEBUG : Check Mpi redundancy, don't work with OpenACC directives


### Fixed
- Correction of the vertical transformation function (in the NEW_S_COORD case) 
  to better handle negative values of h when using wet/dry. 
  The correction guarantees monotonicity of vertical coordinates 
  for all stretching parameters.

- Correct unitialized variables and divisions by zero. 
  [#158](https://gitlab.inria.fr/croco-ocean/croco/-/issues/158)

- DIAGNOSTICS_VRT : misspelling in key name
  [#159](https://gitlab.inria.fr/croco-ocean/croco/-/issues/159)

- Averaged files not written correctly if ntsavg=0 (time-steps missing)
  [#157](https://gitlab.inria.fr/croco-ocean/croco/-/issues/157)

- Avoid unnecessary compiling of partit.F when PARALLEL_FILES is not defined 
  [#125](https://gitlab.inria.fr/croco-ocean/croco/-/issues/125)

- Fix river in run_croco_inter.bash
  [#145](https://gitlab.inria.fr/croco-ocean/croco/-/issues/145)

- TKE : Changing the coefficient 'invG' used in the tridiagonal solving of 
  the tke and gls variables in gls_mixing.F. Using the Patankar trick in the 
  tridiagonal solving consists of moving the Bprod term from the forcing 
  vector S to the diagonal D when the total forcing Sprod + Bprod is 
  negative. Doing this change needs to add to the term Bprod a 'minus' sign 
  and a division by the variable (tke for the tke loop and psi for the gls 
  loop). This was done for the gls loop (invG = 1/psi in this case) but not 
  for the tke loop (ingG = 1). Hence, invG is set to 1/tke in the tke loop.

- EDDY_DIAGS : fixed see issues: 
  [#107](https://gitlab.inria.fr/croco-ocean/croco/-/issues/107), 
  [#117](https://gitlab.inria.fr/croco-ocean/croco/-/issues/117)

- KPP when ANA_DIURNAL_SW is activated : see issue 
  [#115](https://gitlab.inria.fr/croco-ocean/croco/-/issues/115)

- MRL_WCI : correct references to Qiao et al. non-breaking wave diffusion
  (thanks to ChangShui Xia )

- CPL : AGRIF+CPL fixed / compliant with all OASIS versions (fix issue #155),
  corrections for create oasis from pre-existing conditions, 
  correction for LEFTRARU machine (fix issue #149, #135)

- Run issue of RIP and SANDBAR test cases : see issue
  [#130](https://gitlab.inria.fr/croco-ocean/croco/-/issues/130)

- Compilation issue with ifort and openmp : see issue
  [#140](https://gitlab.inria.fr/croco-ocean/croco/-/issues/140)

- Clean files header (fix issue #165)

- Paths in VILAINE test case (fix issue #179)

### Changed

- CI runs through docker container from croco-ocean/croco_docker project
  with both gfortran and intel compilers
  
- OUTPUT : rename Cs_r > Cs_rho, delete redundancy with sc_r and sc_w , see issue 
  [#126](https://gitlab.inria.fr/croco-ocean/croco/-/issues/126) 

- Use of separate vname(s) vector to avoid risk of overlapping 
  [#133](https://gitlab.inria.fr/croco-ocean/croco/-/issues/133)

- CPL : scripts have been modified for easier use, and improved logs 
  and error tracking (fix issue #153) 
  + add management of WW3 extra outputs (spec) 
  + add management of OW_COUPLING_FULL and WAVE_SMFLUX croco cppkeys 
  (solve issues #150 and #168)

### Deprecated

### Removed

### Other
