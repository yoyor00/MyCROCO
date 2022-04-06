Functions description
=====================

.. _functions_description:

Standard functions
------------------

In coupleur_BIOLink, you will find the following subroutines.


   * BIOLink_initialization 

   * BIOLink_init             
                                                                        
   * BIOLink_alloc            
   
   * BIOLink_update
               
   * BIOLink_convarray

   * BIOLink2hydro

   * BIOLink_updateconc_BIO
   
They are used to initialize the arrays or call functions from the biogeochemical models. They also are used to transfer informations from the hydro model and the biogeochemical model.

.. f:subroutine:: BIOLink_initialization(icall)
  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_initialization  ***
  !&E
  !&E ** Purpose : intialization of the biological/tracer model
  !&E
  !&E ** Description : Works in two calls. The first call, with the 
  !&E                  icall = 1 reads the parameter files of the different
  !&E                  models. The second calls use internal functions of 
  !&E                  those models to initialize their variables
  !&E
  !&E ** Called by : main.F90
  !&E
  !&E ** External calls : peptic_param,peptic_alloc_var 
  !&E                     from peptic_initdefine, 
  !&E                     bloom_param,bloom_init_iv 
  !&E                     from peptic_initdefine or 
  !&E                     meteor_param from meteor_initdefine
  !&E
  !&E ** Reference : I do not know
  !&E
  !&E ** History :
  !&E       !  2019-08  (B.Thouvenin) :  creation 
  !&E       !  2022-02  (G. Koenig)   : Commenting 
  !&E---------------------------------------------------------------------
     


.. f:subroutine :: BIOLink_init(ifirst,ilast,jfirst,jlast)
  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_init  ***
  !&E
  !&E ** Purpose : Initialization of tables of variables in biological 
  !&E              models and of table for the helping functions 
  !&E
  !&E ** Description : First it sets the time at the initial time of 
  !&E                  the hydro model + one slow time step, then it
  !&E                  initializes the settling velocities and the
  !&E                  conservation log file. Then it uses biological/
  !&E                  tracer routines to initialize the table of those
  !&E                  models and finally it initializes the tables of 
  !&E                  the helping functions.
  !&E
  !&E ** Called by : main
  !&E
  !&E ** External calls : bloom_userinit of bloom_initdefine,
  !&E                     bloom_wavefile_MANGAbio from bloom,
  !&E                     meteor_read_react from meteor_initdefine
  !&E
  !&E ** Reference : I do not know
  !&E
  !&E ** History :
  !&E       !  2015-07  : Creation by B. Thouvenin I guess 
  !&E       !  2022-02  : Commenting (G. Koenig)
  !&E
  !&E---------------------------------------------------------------------
  
.. f:subroutine :: BIOLink_alloc()
  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_alloc ***
  !&E
  !&E ** Purpose : allocation of common variables tables. Those variables 
  !&E              were declared earlier but their dimensions depends on 
  !&E              the reading of some parametrization files that are read  
  !&E              later in the CROCO model. Once those files are read, we 
  !&E              can allocate the tables their final shapes.
  !&E                            
  !&E
  !&E ** Description : The entire function consist of a succession of 
  !&E                  memory allocation followed by value initialization 
  !&E                  of the tables to dummy values.
  !&E
  !&E ** Called by : coupleur_BIOLink in the second call of the subroutine 
  !&E                   "BIOLink_initialization"
  !&E
  !&E ** External calls :  All the tables are declared as public and 
  !&E                      therefore are accessible even if they are not
  !&E                      given as argument of the subroutine. 
  !&E
  !&E ** Reference : None
  !&E
  !&E ** History : Not sure. Probably created by B. Thouvenin.
  !&E
  !&E--------------------------------------------------------------------- 
  
  
.. f:subroutine :: BIOLink_update(ifirst,ilast,jfirst,jlast, CELL_SURF)
  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_update  ***
  !&E
  !&E ** Purpose : Update of the BIOLink concentration, sources and sinks
  !&E              and helping variables (settling velocities, number of 
  !&E              oyster, etc...) at each time step.
  !&E
  !&E ** Description : 
  !&E
  !&E ** Called by : init
  !&E
  !&E ** External calls : peptic_sksc_wat, peptic_SPMtot_Chla from peptic
  !&E                     bloom_sksc_wat, bloom_eval_diag2d, 
  !&E                     bloom_SPMtot_Chla and bloom_extinction_avg
  !&E                     from bloom 
  !&E                     bloom_wavefile_MANGAbio from bloom
  !&E                     COUPLEUR_PHY_BIO from ECO3M
  !&E                     meteor_sksc_wat, meteor_reac_equi from METeOR
  !&E                     meteor_sksc_wat, meteor_reac_equi from comvars2d
  !&E
  !&E ** Reference : There is first an evaluation to know if the time step 
  !&E                the hydrodynamic model has evolved enough. If yes,
  !&E                there is an update of the biological time, concentration
  !&E                , height of the water column and helping variables. 
  !&E                And finally sometimes we reinitialize the oysters.
  !&E
  !&E ** History :
  !&E       !  2019-08  (B.Thouvenin)  creation 
  !&E       !  2022-02  (G. Koenig) Commenting
  !&E---------------------------------------------------------------------
  
.. f:subroutine :: BIOLink_convarray(ifirst,ilast,jfirst,jlast,WAT_SETTL)
  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_convarray ***
  !&E
  !&E ** Purpose : conversion of 3D or 4D array from hydro model to BIOLink
  !&E
  !&E ** Description : Everything is said in the purpose
  !&E
  !&E ** Called by : BIOLink_update
  !&E
  !&E ** External calls : None
  !&E
  !&E ** Reference :
  !&E
  !&E ** History : Creation by B. Thouvenin ( date unknown)
  !&E              Commenting and editing by G. Koenig (February 2022)
  !&E
  !&E---------------------------------------------------------------------
  
.. f:subroutine:: BIOLink2hydro(ifirst,ilast,jfirst,jlast,WAT_SETTL)     
  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink2hydro ***
  !&E
  !&E ** Purpose : Conversion of setling velocities array from BIOLink 
  !&E              to hydro model 
  !&E
  !&E ** Description : The array in BIOLink are first indexed with the depth,
  !&E                  while the ones of hydro models are first indexed by
  !&E                  horizontal positions. Here we convert the array of 
  !&E                  BIOLink to the hydro model format
  !&E
  !&E ** Called by :  BIOLink_update
  !&E
  !&E ** External calls : None
  !&E
  !&E ** Reference :
  !&E
  !&E ** History : ! Created by B. Thouvenin 
  !&E              ! Commented by G. Koenig ( february 2022)
  !&E
  !&E---------------------------------------------------------------------
 
 
.. f:subroutine:: BIOLink_updateconc_BIO(ifirst,ilast,jfirst,jlast)     
  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_updateconc_BIO ***
  !&E
  !&E ** Purpose : Updating tracers concentration
  !&E
  !&E ** Description : First there is a loop to determine if the conce
  !&E                  ntrations are going to become negative. If is the 
  !&E                  case we limit the sources and sink terms. And then
  !&E                  we compute the advective flux of tracers in the cell
  !&E                  with the hydro time steps.
  !&E
  !&E ** Called by : BIOLink_update
  !&E
  !&E ** External calls :
  !&E
  !&E ** Reference :
  !&E
  !&E ** History : ! Created by B. Thouvenin
  !&E              ! Commented by G. Koenig ( february 2022 )
  !&E
  !&E---------------------------------------------------------------------
  
  
Physics functions
-----------------

In coupleur_BIOLink_physics, you will find the following subroutines.


   * BIOLink_water_column
   
It is used to determine the height of each cell in the water column.

.. f:subroutine:: BIOLink_water_column(ifirst,ilast,jfirst,jlast)

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_water_column ***
  !&E
  !&E ** Purpose : Computation of total water thickness and thickness of 
  !&E              each vertical layer
  !&E
  !&E ** Description : We read the depth and water elevation from the 
  !&E                  hydro model and use it to compute the total water
  !&E                  thickness. Then we use the difference of height
  !&E                  of each cell to compute the thickness of each cell
  !&E
  !&E ** Called by : BIOLink_update
  !&E
  !&E ** External calls : limitation of the subgrid given by MPI ifirst,
  !&E                     ilast, jfirst and jlast
  !&E
  !&E ** Reference :
  !&E
  !&E ** History : Creation by B. Thouvenin ( date unknown)
  !&E              Commenting by G. Koenig ( february 2022)
  !&E
  !&E---------------------------------------------------------------------
  
Helping functions
-----------------

In coupleur_BIOLink_helping, you will find the following subroutines.

   * BIOLink_read_vardiag

   * BIOLink_sinking_rate

   * BIOLink_eval_PAR

   * BIOLink_substance

   * BIOLink_SPMsat_file
   
They are used for auxiliary functions for each submodel.

.. f:subroutine:: BIOLink_read_vardiag
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE BIOLink_read_vardiag  ***
   !&E
   !&E ** Purpose : Reading of the file describing the diagnostic variables
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : BIOLink_update
   !&E
   !&E ** External calls : bloom_init_id,bloom_create_var_diagtracer 
   !&E                     from bloom_initdefine
   !&E                     nb_var_tracerN,nb_var_tracerP from parameters
   !&E                    
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2019-08 (B. Thouvenin) issued from sub_read_vardiag 
   !&E       !  2022-02 (G. Koenig) commenting 
   !&E
   !&E---------------------------------------------------------------------
   
.. f:subroutine:: BIOLink_sinking_rate(ifirst,ilast,jfirst,jlast)

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE BIOLink_sinking_rate  ***
   !&E
   !&E ** Biologic dynamics:  - limiting sinking rate for each variables
   !&E
   !&E ** Purpose :  Updating sinking rates
   !&E              
   !&E ** Description :  It first computes the settling velocity from 
   !&E                   the provided files. If not, it computes it
   !&E                   from the ratio of the cell thickness and BIOLink 
   !&E                   time step. In the case we use PEPTIC, it computes
   !&E                   it from PEPTIC 
   !&E
   !&E ** Called by : BIOLink_update
   !&E
   !&E ** External calls : None
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2019-08 (B. Thouvenin) issued from peptic_dynzwat - verti_quota 
   !&E                  (A. Menesguen, M. Sourrisseau)
   !&E       !  2022-02 (G. Koenig) commenting and editing
   !&E
   !&E---------------------------------------------------------------------
   
.. f:subroutine:: BIOLink_substance(icall)
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE coupl_BIOLink_dim  ***
   !&E
   !&E ** Purpose : Remplacement of substance module for the settlement
   !&E              and interaction with the sediment
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : BIOLink_update
   !&E
   !&E ** External calls : None
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2016-11  (B. Thouvenin) 
   !&E       !  2022-02  (G. Koenig)
   !&E--------------------------------------------------------------------------
   
.. f:subroutine:: BIOLink_eval_PAR(ifirst,ilast,jfirst,jlast,cdate)

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_eval_PAR ***
  !&E
  !&E ** Purpose : Evaluation of solar radiation extinction, attenuation and PAR 
  !&E
  !&E ** Description : 
  !&E
  !&E ** Called by : BIOLink_update
  !&E
  !&E ** External calls : None
  !&E
  !&E ** Reference :
  !&E
  !&E ** History : B. Thouvenin ( date unknown), creation from verti_quota
  !&E              G. Koenig ( February 2022), commenting
  !&E
  !&E---------------------------------------------------------------------
  
.. f:subroutine:: BIOLink_SPMsat_file(ifirst,ilast,jfirst,jlast,icall,forcSPM)
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE BIOLink_SPMsat_file  ***
   !&E
   !&E ** Purpose : special reading SPM satellite file for key_messat
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : BIOLink_update
   !&E
   !&E ** External calls : 
   !&E
   !&E ** Reference :
   !&E
   !&E ** History : !Creation by B. Thouvenin 
   !&E              ! Commenting by G. Koenig ( February 2022)
   !&E
   !&E--------------------------------------------------------------------------
