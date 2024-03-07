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

- OBSTRUCTION : Add a process-based model for 3-dimensional simulation of 
  flow in presence of various obstructions, activated with cpp key #OBSTRUCTION. Obstructions can be rigid or flexible,  and of 3 types : 
  - upward (like seagrass), 
  - downward (like mussel long-line),
  - 3D (like oyster tables)
  This module requires the keys #SOLVE3D, #GLS_MIXING and #GLS_KEPSILON. See issue
  [#123](https://gitlab.inria.fr/croco-ocean/croco/-/issues/123); 
  check merge request and documentation if interested.

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

### Changed

- OUTPUT : rename Cs_r > Cs_rho, delete redundancy with sc_r and sc_w , see issue 
  [#126](https://gitlab.inria.fr/croco-ocean/croco/-/issues/126) 

- Use of separate vname(s) vector to avoid risk of overlapping 
  [#133](https://gitlab.inria.fr/croco-ocean/croco/-/issues/133)

### Deprecated

### Removed

### Other
