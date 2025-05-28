 MODULE comionc4

   !!======================================================================
   !!                   ***  MODULE comionc4  ***
   !!
   !! ** Description : global declarations for I/O netcdf
   !!
   !! ** History :
   !!       !  2002-01  (J.-F. Leroux)  Original code
   !!       !  2011-01  (R. Ramel from Alyotech) update to fortran90/95
   !!
   !!======================================================================

   IMPLICIT NONE
      
   REAL,PARAMETER :: ionc_valmanque=1.e+37
   INTEGER,PARAMETER ::  ionc_longtabfich=10000
   INTEGER,PARAMETER ::  ionc_longnomfich=200
   CHARACTER(LEN=*),PARAMETER :: ionc_originet='01-JAN-1900 00:00:00'
   CHARACTER(LEN=*),PARAMETER :: ionc_longnamet='time in seconds (UT)'
   CHARACTER(LEN=*),PARAMETER :: ionc_unitst='seconds since 1900-01-01T00:00:00Z'
   CHARACTER(LEN=*),PARAMETER :: ionc_pnomlat  ='nj'
   CHARACTER(LEN=*),PARAMETER :: ionc_pnomlon  ='ni'
   CHARACTER(LEN=*),PARAMETER :: ionc_pnomlat_u='nj_u'
   CHARACTER(LEN=*),PARAMETER :: ionc_pnomlon_u='ni_u'
   CHARACTER(LEN=*),PARAMETER :: ionc_pnomlat_v='nj_v'
   CHARACTER(LEN=*),PARAMETER :: ionc_pnomlon_v='ni_v'
   CHARACTER(LEN=*),PARAMETER :: ionc_pnomlat_f='nj_f'
   CHARACTER(LEN=*),PARAMETER :: ionc_pnomlon_f='ni_f'
   CHARACTER(LEN=*),PARAMETER :: ionc_pnomtime='time'
   CHARACTER(LEN=*),PARAMETER :: ionc_pnomz='level'
   CHARACTER(LEN=*),PARAMETER :: ionc_pnomz_w='level_w'
   CHARACTER(LEN=*),PARAMETER :: ionc_pnomz_s='level_s'
   CHARACTER(LEN=*),PARAMETER :: ionc_pnomtraj='traj'

   INTEGER :: ionc_limi,ionc_lima,ionc_pasi,ionc_global_imin,ionc_global_imax
   INTEGER :: ionc_ljmi,ionc_ljma,ionc_pasj,ionc_global_jmin,ionc_global_jmax
   INTEGER :: ionc_lkmi,ionc_lkma,ionc_pask
   INTEGER :: ionc_nfich
   INTEGER :: ionc_init_ok
   CHARACTER(LEN=ionc_longnomfich),DIMENSION(ionc_longtabfich) :: ionc_nomfich
   CHARACTER(LEN=30) :: ionc_rout
   CHARACTER(LEN=30) :: ionc_nomlat
   CHARACTER(LEN=30) :: ionc_nomlon
   CHARACTER(LEN=30) :: ionc_nomtime
   CHARACTER(LEN=30) :: ionc_nomz
!   CHARACTER(LEN=30) :: ionc_nomz_w
   CHARACTER(LEN=30) :: ionc_nomtraj
   CHARACTER(LEN=30) :: ionc_nomsta
   CHARACTER(LEN=30) :: ionc_nomlattraj
   CHARACTER(LEN=30) :: ionc_nomlontraj

   CHARACTER(LEN=30),PARAMETER :: ionc_conv_data_type='OCO oriented grid'
   CHARACTER(LEN=30),PARAMETER :: ionc_conv_format_version='1.3.1'
   CHARACTER(LEN=30),PARAMETER :: ionc_conv_conventions='CF-1.6 OCO-1.3.1 COMODO-1.0'
   CHARACTER(LEN=30),PARAMETER :: ionc_conv_netcdf_version='4.1.2'
   CHARACTER(LEN=30),PARAMETER :: ionc_conv_product_version='1.0'
   CHARACTER(LEN=30),PARAMETER :: ionc_conv_software_version='MARS V9.06'
   CHARACTER(LEN=30),PARAMETER :: ionc_conv_references='http://www.previmer.org/'
   CHARACTER(LEN=30),PARAMETER :: ionc_conv_easting='longitude'
   CHARACTER(LEN=30),PARAMETER :: ionc_conv_northing='latitude'
   CHARACTER(LEN=30),PARAMETER :: ionc_conv_grid_projection='n/a'
   CHARACTER(LEN=30),PARAMETER :: ionc_conv_contact='cdoco-exploit@ifremer.fr'
   CHARACTER(LEN=30),PARAMETER :: ionc_conv_minimum_altitude='0.0 m'
   CHARACTER(LEN=30),PARAMETER :: ionc_conv_maximum_altitude='0.0 m'
   CHARACTER(LEN=30),PARAMETER :: ionc_conv_altitude_resolution='n/a'
   CHARACTER(LEN=44),PARAMETER :: ionc_conv_distribution='Data restrictions: for registered users only'
   CHARACTER(LEN=30),PARAMETER :: ionc_conv_quality_index='1'

   INTEGER,DIMENSION(ionc_longtabfich) :: ionc_idfich

   ! codes d'erreur:

   INTEGER,PARAMETER :: ionc_errcreation = -1000
   INTEGER,PARAMETER :: ionc_errfich     = -1001 
   INTEGER,PARAMETER :: ionc_err_i       = -1002
   INTEGER,PARAMETER :: ionc_err_j       = -1003
   INTEGER,PARAMETER :: ionc_err_k       = -1004
   INTEGER,PARAMETER :: ionc_err_t       = -1005

   ! number of bits of the packed
   INTEGER,        PARAMETER   :: nb_bits = 16

 END MODULE comionc4      
