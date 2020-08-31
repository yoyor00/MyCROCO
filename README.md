# Online Analysis Module
## Spectral and Wavelet analysis 
### Objective : implement an efficient interface between Croco and the Online Analysis Module.
**Authors** F. Auclair and others.

- Preliminary version of the interface between Croco and the Online Analysis code. 
-- The interface is first based on calls to the subroutine croco_oa within the main.F and step.F Croco subroutines (calls are performed at the initialization step - 3Dfast now - and at subsequent 3D-fast time steps). The croco_oa.F90 routine includes several Croco header files (grid.h, ocean2d/3d.h,...) so that the domain and vertical grid variables and parameters (lon-lat geographical or curvilinear coordinates, ocean thickness, sub-domain indices,...) can be passed as arguments to the Stand-alone Online Analysis subroutines init_oa and main_oa.
-- Within the Stand Alone module, the calculation of the correlation product between the psi_oa function (Fourier/Wavelet) and the requested ocean field value (cf namelist_oa and requested variable-configuration-analysis) is performed thanks to the var_oa external real function which returns (possibly at each time step) the value of the Croco ocean field value at a given domain location. Indeed, the var_oa external real function also includes several Croco header files (grid.h, ocean2d/3d.h,...) since most of the Croco ocean fields are not passed as arguments to the Stand-alone OA code.
-- Finally, the Online Analysis data is output by the mean of the XIOS2 facility thanks to few changes implmented in the send send_xios_diags.F routine (eg, use of the specific Online Analysis modules to see the OA data).

- The Test case with an ad hoc test function, e.g., utest = 5. * cos(omega_S2 * t) + 3. cos( omega_S4 * t ) is now working (see bug corrected in croco_oa and comments related to the stand-alone Online Analysis code where rhp_t argument is only used for tracking isopycnal movement)

- The IGW test case using pseudo S2 tidal forcings (12h) has been launched over 16 days performing Morlet S2 and S4 Online Analysis every 1h (correlation product sampling every time step over a convolution window which is 4 times the Wavelet period). 

###TODO
- It is necessary to improve the current Croco - Online Analysis Interface :
-- the objective is to replace croco_oa by a "loop on tile" construct specific to Croco with dedicated modules for the Croco to OA interface and OA to Croco- XIOS2 inteface. This will solve the duplication of the Croco arrays in memory with the current interface (Croco arrays are currently passed as arguments to the Online Analysis routines init_oa and main_oa with reduced dimensions).
- A check regarding the treatment of the precision of complex numbers (intel/gnu compilation flag -r8, croco mpc.F preproc applied to F90 sources, F90 fixed format sources or not)
- Tests must be conducted with GNU compilers (at this stage only only Intel compilers were tested). 
