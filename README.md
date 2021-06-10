# Online Analysis Module
## Spectral and Wavelet analysis 
### Objective : implement an efficient interface between Croco and the Online Analysis Module.
**Authors** B. Lemieux-Dudon, F. Auclair and C. Nguyen

The Online Analysis module enables to perform online Wavelet and spectral analysis targetting specific simulated ocean fields.
A preliminary version of this interface was set up in Spring/Summer 2020 for the Croco ocean model.

This preliminary version of the Online Analysis - Croco interface provides appropriate results when applied to the two Test Cases :
1. A uniform field varyrying as a linear combination of cosine functions, e.g.,
F(i,j,k,t) = 5 cos( W<sub>S2</sub> t ) + 3 cos( W<sub>S4</sub> t )

[Validation slides](TEST_CASES/IGW_OA/MorletWaveletAnalysis/OnlineAnalysis_Croco_Interface_validation.pdf) ( right-click, 'Save link as...' on the name )

2. The IGW academic test case with a continental slope at its Eastern boundary and with a pseudo-S2 (12h) tidal forcing at its Western boundary.

![IGW Internal Tides - Online Analysis of rho with Morlet S2+S4](TEST_CASES/IGW_OA/MorletWaveletAnalysis/IGW_tidal_forcing_S2/IGW_S2_rho_oa_MorletS2S4_200-240.gif)

<!-- [![IGW Internal Tides](https://img.youtube.com/vi/AB84oYuNDKA/0.jpg)](https://www.youtube.com/watch?v=AB84oYuNDKA) -->
<!-- [![IGW Internal Tides](https://img.youtube.com/vi/tsKXy_j93yA/0.jpg)](https://www.youtube.com/watch?v=tsKXy_j93yA) -->
<!--
[![IGW Internal Tides - Full Online Analysis of rho with Morlet S2+S4 16 days (youtube video) ]()](https://www.youtube.com/watch?v=tsKXy_j93yA)
-->
- IGW Internal Tides - Full Online Analysis of rho with Morlet S2+S4 16 days (youtube video) : https://www.youtube.com/watch?v=tsKXy_j93yA
- Plots were done installing and modifying the pyroms tools for Croco : https://github.com/ESMG/pyroms.

**Current status**

1. Part of the Online Analysis are now operational within the Croco model but some options still need to be tested.
So far :
- only the online Wavelet transform with the Morlet atom was fully tested. 
- The Croco code and Online Analysis module were compiled with the intel compilers (18.0.2) only.
- The IGW test case ran over 16 days performing Morlet S2 and S4 Online Analysis every 1h or every 600s using 8 MPI subdomains.
- The simulation fields and the Online Analysed fields have been output using XIOS2 with a modified version of send_xios_diags.F routine and properly parametrized xml files. 

2. The test case is now working (see bug corrected in croco_oa and comments related to the stand-alone Online Analysis code where rhp_t argument is only used for tracking isopycnal movement).


**The preliminary version of the Online Analysis - Croco interface must be improved to solve several issues.**

The interface is currently based :
1. on calls to the subroutine croco_oa within the main.F and step.F Croco subroutines (calls are performed at the initialization step - 3D now - and at subsequent 3D time steps). The croco_oa.F90 routine includes several Croco header files (grid.h, ocean2d/3d.h,...) so that the horizontal domain and vertical grid variables and parameters (lon-lat geographical or curvilinear coordinates, ocean thickness, sub-domain indices,...) can be passed as arguments to the Stand-alone Online Analysis subroutines init_oa and main_oa.
2. Within the Stand Alone module, the calculation of the correlation product between the atom function (e.g., Fourier, Wavelet,...) and the requested ocean field (see the user namelist and requested variable-configuration-analysis) is performed thanks to the var_oa external real function which returns the value of the Croco ocean field at a given location (possibly at each time step). Indeed, the var_oa external real function also includes several Croco header files (grid.h, ocean2d/3d.h,...) since most of the Croco ocean fields are not passed as arguments to the Stand-alone OA code.
3. Finally, the Online Analysis data is output using the XIOS2 facility (few changes were implemented in the send send_xios_diags.F routine with the use of specific Online Analysis modules).

The current interface involve the duplication of the Croco arrays in memory (see the Stand-alone algo. with Croco arrays passed as arguments to init_oa and main_oa with reduced dimensions).

**TODO**
- Improve the current Croco - Online Analysis Interface to solve the memory issue and to conform with some Croco pre-requisite (double OpenMPI - OpenMP paralelization, shared header files, the way to handle precision with source pre-processing ) :
  -- the objective is to replace croco_oa by a "loop on tile" construct specific to Croco with dedicated modules for the Croco to OA interface and OA to Croco-XIOS2 interface.
- A check regarding the treatment of the precision of complex numbers (intel/gnu compilation flag -r8, croco mpc.F preproc applied to F90 sources, F90 fixed format sources or not).
- Tests must be conducted with GNU compilers (at this stage only only Intel compilers were tested).
- Tests must be conducted with mask options,...etc.
