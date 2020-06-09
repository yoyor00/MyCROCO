:=============================================================================
:         ENREGISTREMENT DES PARAMETRES ET DES SUBSTANCES SOUS FORME DE NAMELIST
:         INPUT PARAMETERS AND SUBSTANCES LIST NAMELIST
:
:      with  CPP key : SUBSTANCE (without CPP key MUSTANG)
:=============================================================================
:
:-----------------------------------------------------------------------------
:
: nmlnbvar : nombres de substances de chaque type a definir (autre que T et S)
:          : number of each type of substance to be defined (other than T & S)
:    other variables can be created in addition with keys : contaminant, peptic, P_tracer, N_tracer
:-----------------------------------------------------------------------------
:    nv_dis : number of dissolved susbtances
:    nv_ncp : number of  Non Constitutive Particulate subtances
:    nv_fix : number of  fixed susbtances (not advected)
:    nv_bent : number of benthic susbtances
:   

 &nmlnbvar
    nv_dis=2
    nv_ncp=2
    nv_fix=1
    nv_bent=0 /

:-----------------------------------------------------------------------------
: Namelist pour chaque type de substance  :
: Namelist for each type of substance 
:-----------------------------------------------------------------------------

:-----------------------------------------------------------------------------
:
: nmlpartnc : NON constitutive substances with NOT defined key_sed_MUSTANG
:-----------------------------------------------------------------------------
:    name_var_n() : name of variable
:    long_name_var_n() : long name of variable
:    standard_name_var_n() : standard name of variable
:    unit_var_n() : unit of concentration of variable
:    flx_atm_n() : uniform atmospherical deposition (unit/m2/s) 
:    cv_rain_n() : concentration in rainwater (unit/m3 of water)
:    cini_wat_n() : initial concentration in water column (unit/m3)
:    cobc_wat_n() : boundaries uniform and constant concentration  (unit/m3)
:    cini_air_n() : initial concentration in air
:    l_out_subs_n() : saving in output file if TRUE
:    init_cv_name_n() : name of substance read from initial condition file
:    obc_cv_name_n() : name of substance read from obc file
:    ws_free_min_n() : minimum setling velocity (m/s)
:    ws_free_max_n() : maximum setling velocity (m/s)
:

&nmlpartnc
    name_var_n(1)='MES1'
       long_name_var_n(1)='MES1'
       standard_name_var_n(1)='MES1'
       unit_var_n(1)='g.l-1'
       flx_atm_n(1)=0.0
       cv_rain_n(1)=0.0
       cini_wat_n(1)=1.
       cobc_wat_n(1)=1.
       cini_air_n(1)=0.0
       l_out_subs_n(1)=.TRUE.
       init_cv_name_n(1)=''
       obc_cv_name_n(1)=''
       ws_free_min_n(1)=0.0004
       ws_free_max_n(1)=0.001 
          
    name_var_n(2)='MES2'
       long_name_var_n(2)='MES2'
       standard_name_var_n(2)='MES2'
       unit_var_n(2)='micromole.l-1'
       flx_atm_n(2)=0.0
       cv_rain_n(2)=0.0
       cini_wat_n(2)=10.
       cobc_wat_n(2)=1.
       cini_air_n(2)=0.0
       l_out_subs_n(2)=.TRUE.
       init_cv_name_n(2)=''
       obc_cv_name_n(2)=''
       ws_free_min_n(2)=1.e-7
       ws_free_max_n(2)=1.e-7 /


 : nmlvardiss : DISSOLVED SUBSTANCES
:    name_var_n() : name of variable
:    long_name_var_n() : long name of variable
:    standard_name_var_n() : standard name of variable
:    unit_var_n() : unit of concentration of variable
:    flx_atm_n() : uniform atmospherical deposition (unit/m2/s) 
:    cv_rain_n() : concentration in rainwater (unit/m3 of water)
:    cini_wat_n() : initial concentration in water column (unit/m3)
:    cobc_wat_n() : boundaries uniform and constant concentration  (unit/m3)
:    cini_air_n() : initial concentration in air
:    l_out_subs_n() : saving in output file if TRUE
:    init_cv_name_n() : name of substance read from initial condition file
:    obc_cv_name_n() : name of substance read from obc file
:-----------------------------------------------------------------------------

&nmlvardiss
    name_var_n(1)='ammonium'
       long_name_var_n(1)='ammonium'
       standard_name_var_n(1)='mole_concentration_of_ammonium_in_sea_water'
       unit_var_n(1)='micromole.l-1'
       flx_atm_n(1)=0.0
       cv_rain_n(1)=0.0
       cini_wat_n(1)=100
       cobc_wat_n(1)=100.
       cini_air_n(1)=0.0
       l_out_subs_n(1)=.TRUE.
       init_cv_name_n(1)='ammonium'
       obc_cv_name_n(1)='ammonium'
       
    name_var_n(2)='nitrate'
       long_name_var_n(2)='nitrate'
       standard_name_var_n(2)='mole_concentration_of_nitrate_in_sea_water'
       unit_var_n(2)='micromole.l-1'
       flx_atm_n(2)=0.0
       cv_rain_n(2)=0.0
       cini_wat_n(2)=1.0
       cobc_wat_n(1)=100.
       cini_air_n(2)=0.0
       l_out_subs_n(2)=.FALSE.
       init_cv_name_n(2)='nitrate'
       obc_cv_name_n(2)='nitrate'  /

 : nmlvarfix : FIXED SUBSTANCES (not advected)
:-----------------------------------------------------------------------------
:     

 &nmlvarfix
    name_var_fix(1)='cumulative_nanoflagellate_carbon_production'
       long_name_var_fix(1)='cumulative_nanoflagellate_carbon_production'
       standard_name_var_fix(1)='cumulative_nanoflagellate_production_expressed_as_carbon_in_sea_water'
       unit_var_fix(1)='g.m-3_from_01jan'
       cini_wat_fix(1)=0.0
       l_out_subs_fix(1)=.TRUE.
       init_cv_name_fix(1)='none'  /


  : nmlbent : used only if defined key_benthic 
:-----------------------------------------------------------------------------
: & nmlvarbent


