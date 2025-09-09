module ionc4
#include "cppdefs.h"
   !!======================================================================
   !!                   ***  MODULE ionc4  ***
   !!
   !! ** Description : Library for I/O netcdf
   !!
   !! ** History :
   !!       !  2002-01  (J.-F. Le Roux)  Original code io_netcdf.F90
   !!       !  2011-01  (R. Ramel from Alyotech) upgrate to fortran90/95
   !!       !  2012-01  (J.-F. Le Roux) add subroutines to read traj
   !!       !  2013-07 (A. Thevenin - CERFACS) coupling using OASIS3-MCT 
   !!
   !!======================================================================

!   use module_lagrangian
   use netcdf
   use typeSizes
   use comionc4

!#include "set_global_definitions.h"

   implicit none

#ifdef key_oasis
   INTEGER, PUBLIC :: MPI_COMM_MARS
#endif


   interface ionc4_createfile
      module procedure ionc4_createfile_1d_double
      module procedure ionc4_createfile_1d_real
      module procedure ionc4_createfile_2d_double
      module procedure ionc4_createfile_2d_real
   end interface ionc4_createfile

   interface ionc4_createvar
      module procedure ionc4_createvar_int
      module procedure ionc4_createvar_real
      module procedure ionc4_createvar_double
   end interface ionc4_createvar
   
   interface ionc4_write_zxyt
      module procedure  ionc4_write_zxyt_real
      module procedure  ionc4_write_zxyt_double
   end interface ionc4_write_zxyt
   
   interface ionc4_read_zxyt
      module procedure ionc4_read_zxyt_real
      module procedure ionc4_read_zxyt_double
   end interface ionc4_read_zxyt
   
   interface ionc4_write_xyt
      module procedure  ionc4_write_xyt_int
      module procedure  ionc4_write_xyt_real
      module procedure  ionc4_write_xyt_double
   end interface ionc4_write_xyt
   
   interface ionc4_read_xyt
      module procedure ionc4_read_xyt_int
      module procedure ionc4_read_xyt_real
      module procedure ionc4_read_xyt_double
   end interface ionc4_read_xyt
   
   interface ionc4_write_zxy
      module procedure ionc4_write_zxy_real
      module procedure ionc4_write_zxy_double
   end interface ionc4_write_zxy
   
   interface ionc4_read_zxy
      module procedure ionc4_read_zxy_real
      module procedure ionc4_read_zxy_double
   end interface ionc4_read_zxy
   
   interface ionc4_write_xy
      module procedure  ionc4_write_xy_real
      module procedure  ionc4_write_xy_double
   end interface ionc4_write_xy
   
   interface ionc4_read_xy
      module procedure ionc4_read_xy_real
      module procedure ionc4_read_xy_double
   end interface ionc4_read_xy
   
   interface ionc4_read_zxt
      module procedure ionc4_read_zxt_real
   end interface ionc4_read_zxt
   
   interface ionc4_write_xt
      module procedure ionc4_write_xt_real
   end interface ionc4_write_xt
   
   interface ionc4_read_xt
      module procedure ionc4_read_xt_real
   end interface ionc4_read_xt
   
   interface ionc4_write_lon
      module procedure ionc4_write_lon_real
      module procedure ionc4_write_lon_double
   end interface ionc4_write_lon
   
   interface ionc4_read_lon
      module procedure ionc4_read_lon_real
      module procedure ionc4_read_lon_double
   end interface ionc4_read_lon
   
   interface ionc4_write_lat
      module procedure ionc4_write_lat_real
      module procedure ionc4_write_lat_double
   end interface ionc4_write_lat
   
   interface ionc4_read_lat
      module procedure ionc4_read_lat_real
      module procedure ionc4_read_lat_double
   end interface ionc4_read_lat
   
   interface ionc4_write_sig
      module procedure  ionc4_write_sig_real
      module procedure  ionc4_write_sig_double
   end interface ionc4_write_sig
   
   interface ionc4_write_sigw
      module procedure  ionc4_write_sigw_real
      module procedure  ionc4_write_sigw_double
   end interface ionc4_write_sigw
   
   interface ionc4_write_z
      module procedure  ionc4_write_z_real
   end interface ionc4_write_z
   
   interface ionc4_read_z
      module procedure  ionc4_read_z_real
   end interface ionc4_read_z
   
   interface ionc4_createfile_sta
      module procedure ionc4_createfile_sta_double
      module procedure ionc4_createfile_sta_real
   end interface ionc4_createfile_sta
   
   interface ionc4_createvar_traj
      module procedure ionc4_createvar_traj_int
      module procedure ionc4_createvar_traj_real
      module procedure ionc4_createvar_traj_double
   end interface ionc4_createvar_traj

   interface ionc4_write_traj
      module procedure  ionc4_write_traj_int
      module procedure  ionc4_write_traj_real
      module procedure  ionc4_write_traj_double
   end interface ionc4_write_traj
   
   interface ionc4_write_trajt
      module procedure  ionc4_write_trajt_int
      module procedure  ionc4_write_trajt_real
      module procedure  ionc4_write_trajt_double
   end interface ionc4_write_trajt

   interface ionc4_read_traj
      module procedure  ionc4_read_traj_int
      module procedure  ionc4_read_traj_real
      module procedure  ionc4_read_traj_double
   end interface ionc4_read_traj
   
   interface ionc4_read_trajt
      module procedure  ionc4_read_trajt_int
      module procedure  ionc4_read_trajt_real
      module procedure  ionc4_read_trajt_double
   end interface ionc4_read_trajt
   
   interface ionc4_write_t
      module procedure  ionc4_write_t_real
      module procedure  ionc4_write_t_double
   end interface ionc4_write_t
   
   interface ionc4_read_t
      module procedure  ionc4_read_t_real
      module procedure  ionc4_read_t_double
   end interface ionc4_read_t
   
   interface ionc4_write
      module procedure  ionc4_write_int
      module procedure  ionc4_write_real
      module procedure  ionc4_write_double
   end interface ionc4_write
   
   interface ionc4_read
      module procedure  ionc4_read_int
      module procedure  ionc4_read_real
      module procedure  ionc4_read_double
   end interface ionc4_read
   
   interface ionc4_vatt_fill
      module procedure ionc4_vatt_fill_real
      module procedure ionc4_vatt_fill_double
   end interface ionc4_vatt_fill
   
   interface ionc4_vatt_valid_min
      module procedure ionc4_vatt_valid_min_real
      module procedure ionc4_vatt_valid_min_double
   end interface ionc4_vatt_valid_min
   
   interface ionc4_vatt_valid_max
      module procedure ionc4_vatt_valid_max_real
      module procedure ionc4_vatt_valid_max_double
   end interface ionc4_vatt_valid_max
   
   interface ionc4_gatt
      module procedure ionc4_gatt_char
      module procedure ionc4_gatt_int
   end interface ionc4_gatt

   interface ionc4_read_subzxyt
      module procedure ionc4_read_subzxyt_real
      module procedure ionc4_read_subzxyt_double
   end interface ionc4_read_subzxyt

   interface ionc4_read_sublzxyt
      module procedure ionc4_read_sublzxyt_real
   end interface ionc4_read_sublzxyt

   interface ionc4_read_subxyt
      module procedure ionc4_read_subxyt_real
      module procedure ionc4_read_subxyt_double
   end interface ionc4_read_subxyt

   interface ionc4_read_subxy
      module procedure ionc4_read_subxy_real
      module procedure ionc4_read_subxy_double
   end interface ionc4_read_subxy
   
   interface ionc4_read_sublxyt
      module procedure ionc4_read_sublxyt_real
   end interface ionc4_read_sublxyt
   
   interface ionc4_createvar_sta
      module procedure ionc4_createvar_sta_real
      module procedure ionc4_createvar_sta_char
   end interface ionc4_createvar_sta
   
   interface ionc4_write_sta
      module procedure ionc4_write_sta_char
      module procedure ionc4_write_sta_real
      module procedure ionc4_write_sta_double
   end interface ionc4_write_sta

   interface ionc4_write_stat
      module procedure ionc4_write_stat_real
      module procedure ionc4_write_stat_double
   end interface ionc4_write_stat

   interface ionc4_write_zsta
      module procedure ionc4_write_zsta_real
      module procedure ionc4_write_zsta_double
   end interface ionc4_write_zsta

   interface ionc4_write_zstat
      module procedure ionc4_write_zstat_real
      module procedure ionc4_write_zstat_double
   end interface ionc4_write_zstat
 

   contains
   
! ***************************************************************
! * subroutine ionc4_init                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    21/02/01                                *
! * derniere modif :    21/02/01 par jfleroux                   *
! *                     30/07/13 par A. Thevenin (CERFACS)      *
! *      *
! * Role : initialisation du module ionc       *
! *      *
! ***************************************************************
#ifdef key_oasis
   subroutine ionc4_init(id_lcomm)

   INTEGER, INTENT(IN) :: id_lcomm     ! local communicator provided by OASIS
                                       ! to be used instead of MPI_COMM_WORLD
   MPI_COMM_MARS = id_lcomm 

#else
   subroutine ionc4_init()
#endif
   
   ionc_nfich = 0
   ionc_nomlat = "latitude"
   ionc_nomlon = "longitude"
   ionc_nomtime = ionc_pnomtime
   ionc_nomz = ionc_pnomz
   ionc_nomtraj = ionc_pnomtraj
   ionc_nomsta = 'sta'
   ionc_init_ok = 1
   
   return
   end subroutine ionc4_init
   
   
! ***************************************************************
! * SUBROUTINE ionc4_createfile_2d_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    20/02/01 par jfleroux                   *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)          *
!                       30/07/13 par A. Thevenin (CERFACS)      *                       
! *      *
! * Role : creer le fichier Netcdf qui contiendra les donnees      *
! *         Initialiser les donnees ne dependant pas du time :*
! *      - lon,lat*
! *      - SIG(k)*
! *         initialise les variables correspondants a la loupe*
! *         initialise les unites de temps et l'origine          *
! *        aux valeurs par defaut                               *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier : chemin complet du fichier a *
! *      creer*
! *      - dg   :               *
! *      - g0   :               *
! *      - dfi  :               *
! *      - fi0  :               *
! *      - limi : limite min de l'indice i de la loupe   *
! *      - lima : limite max de l'indice i de la loupe   *
! *      - pasi : pas de l'indice i de la loupe          *
! *      - ljmi : limite min de l'indice j de la loupe   *
! *      - ljma : limite max de l'indice j de la loupe*
! *      - pasj : pas de l'indice j de la loupe          *
! *      - lkmi : limite min de l'indice k de la loupe   *
! *      - lkma : limite max de l'indice k de la loupe   *
! *      - pask : pas de l'indice k de la loupe          *
! *      - kmax : borne max pour k                       *
! *      0 si 2D*
! *      - dim_time : nombre de pas de time que le *
! *           fichier contiendra si dim_time =0  *
! *           la dimension time est unlimited    *
! *           si dim_time = -1 pas de dim temps  *
! *      - SIG : niveaux           *
! ***************************************************************
   subroutine ionc4_createfile_2d_real(nom_fichier,lon2d,lat2d, &
                                    limi,lima,pasi,ljmi,ljma,pasj, &
                                    lkmi,lkma,pask,kmax,dim_time,  &
                                    l_out_nc4par,comm_active,zaxis_var)
#ifdef MPI   
   include 'mpif.h'
#endif
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*), intent(in) :: nom_fichier
   integer,         intent(in)  :: dim_time
   integer,         intent(in)  :: limi,lima,pasi,ljmi,ljma,pasj,lkmi,lkma,pask,kmax
   integer,intent(in),optional  :: comm_active
   real(kind=4),dimension(INT((lima-limi)/pasi)+1,INT((ljma-ljmi)/pasj)+1),intent(in) :: lon2d,lat2d
   real(kind=4),dimension(kmax),intent(in),optional :: zaxis_var
   logical,intent(in),optional  :: l_out_nc4par
   
! ******** VARIABLES DE TRAVAIL **************

   logical :: l_useless
   integer :: i,communicator
   integer :: dim_lon,dim_lat,dim_k
   
   integer :: nc_id,lat_id,lon_id,time_id,sig_id
   
! ** Declaration des IDS identifiant les dimensions des variables.
   integer :: dim_lon_id,dim_lat_id,dim_time_id,dim_k_id,dim_w_id,ni_id,nj_id, &
              ni_u_id,nj_u_id,ni_v_id,nj_v_id,ni_f_id,nj_f_id,sig_w_id, &
              dim_lon_u_id,dim_lat_u_id,dim_lon_v_id,dim_lat_v_id,dim_lon_f_id,dim_lat_f_id
   
! ** Declarations des tableaux qui contiendront les IDS des dimensions des variables
   integer,dimension(1) :: tab_dim_time,tab_dim_sigma
   integer,dimension(2) :: tab_dim_lon,tab_dim_lat
   integer  :: nc_err,ik,indice,cmode,cache_size
   real(kind=4),dimension(:),allocatable :: ini,inj
   real(kind=4),dimension(:,:),allocatable :: rlon,rlat
   real(kind=4),dimension(:),  allocatable :: rsig

! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_createfile_2d_real"
   
! *** Recuperation des dimensions
   dim_lon = INT((lima-limi)/pasi) +1
   dim_lat = INT((ljma-ljmi)/pasj) +1
   allocate(ini(dim_lon))
   allocate(inj(dim_lat))
   allocate(rlon(dim_lon,dim_lat))
   allocate(rlat(dim_lon,dim_lat))
   cache_size = 10*dim_lon*dim_lat

   if (kmax .GT. 0) then
      dim_k   = (lkma-lkmi)/pask +1
      allocate(rsig(dim_k))
      cache_size = cache_size*dim_k
   endif

   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      
      ionc_nfich = ionc_nfich + 1
      if (ionc_nfich .gt. ionc_longtabfich) then
         call ionc4_err(ionc_errcreation,ionc_rout, ' ' ,nom_fichier)
         stop
      else
         ionc_nomfich(ionc_nfich) = nom_fichier
         cmode = nf90_share
#ifdef MPI
         if (present(comm_active)) then
            communicator=comm_active
         else
#ifdef key_oasis
            communicator=MPI_COMM_MARS
#else
            communicator=MPI_COMM_WORLD
#endif
         endif
         if (present(l_out_nc4par)) then
            if (l_out_nc4par) then
               !print*,'Using parallel I/O features !!!' 
               cmode = ior(NF90_NETCDF4, NF90_CLASSIC_MODEL)
               cmode = ior(cmode, NF90_MPIIO)
               nc_err = nf90_create(nom_fichier, cmode, nc_id, cache_size=cache_size, &
                                    comm=communicator, info=MPI_INFO_NULL)
            else
               nc_err = nf90_create(nom_fichier, cmode, nc_id)
            endif
         else
            nc_err = nf90_create(nom_fichier, cmode, nc_id)
         endif
#else
         if (present(comm_active)) communicator=comm_active  ! useless for portability only
         if (present(l_out_nc4par)) l_useless = l_out_nc4par  ! useless for portability only
         nc_err = nf90_create(nom_fichier, cmode, nc_id)
#endif
         call ionc4_err(nc_err,ionc_rout, 'nf90_create',nom_fichier)
         
         ionc_idfich(ionc_nfich) = nc_id
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlon),dim_lon,dim_lon_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_pnomlon)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlat),dim_lat,dim_lat_id)
         call ionc4_err(nc_err,ionc_rout, 'nf90_def_dim',ionc_pnomlat)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlon_u),dim_lon,dim_lon_u_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_pnomlon_u)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlat_u),dim_lat,dim_lat_u_id)
         call ionc4_err(nc_err,ionc_rout, 'nf90_def_dim',ionc_pnomlat_u)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlon_v),dim_lon,dim_lon_v_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_pnomlon_v)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlat_v),dim_lat,dim_lat_v_id)
         call ionc4_err(nc_err,ionc_rout, 'nf90_def_dim',ionc_pnomlat_v)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlon_f),dim_lon,dim_lon_f_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_pnomlon_f)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlat_f),dim_lat,dim_lat_f_id)
         call ionc4_err(nc_err,ionc_rout, 'nf90_def_dim',ionc_pnomlat_f)
         
         
         if (kmax .GT. 0) then
            nc_err = nf90_def_dim(nc_id,trim(ionc_nomz),dim_k,dim_k_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_nomz)
    !        nc_err = nf90_def_dim(nc_id,trim(ionc_pnomz_w),dim_k,dim_w_id)
    !        call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_pnomz_w)
         endif
         
         if (dim_time .ne. -1) then
            nc_err = nf90_def_dim(nc_id,trim(ionc_nomtime),dim_time,dim_time_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_nomtime)
         endif
         
         
! ******** On definit les variables  **********
        
         tab_dim_lon(1) = dim_lon_id
         tab_dim_lat(1) = dim_lat_id
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlon), nf90_real, tab_dim_lon(1), ni_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlon))
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlat), nf90_real, tab_dim_lat(1), nj_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlat))
         
         tab_dim_lon(1) = dim_lon_u_id
         tab_dim_lat(1) = dim_lat_u_id
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlon_u), nf90_real, tab_dim_lon(1), ni_u_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlon_u))
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlat_u), nf90_real, tab_dim_lat(1), nj_u_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlat_u))
         
         tab_dim_lon(1) = dim_lon_v_id
         tab_dim_lat(1) = dim_lat_v_id
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlon_v), nf90_real, tab_dim_lon(1), ni_v_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlon_v))
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlat_v), nf90_real, tab_dim_lat(1), nj_v_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlat_v))
         
         tab_dim_lon(1) = dim_lon_f_id
         tab_dim_lat(1) = dim_lat_f_id
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlon_f), nf90_real, tab_dim_lon(1), ni_f_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlon_f))
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlat_f), nf90_real, tab_dim_lat(1), nj_f_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlat_f))
         
         tab_dim_lat(1) = dim_lon_id
         tab_dim_lat(2) = dim_lat_id
         nc_err = nf90_def_var(nc_id, trim(ionc_nomlat), nf90_real, tab_dim_lat, lat_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomlat))
         
         tab_dim_lon(1) = dim_lon_id
         tab_dim_lon(2) = dim_lat_id
         nc_err = nf90_def_var(nc_id, trim(ionc_nomlon), nf90_real, tab_dim_lon, lon_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomlon))
         
         if (dim_time .ne. -1) then
            tab_dim_time(1) = dim_time_id
            nc_err = nf90_def_var(nc_id, trim(ionc_nomtime), nf90_real, tab_dim_time, time_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomtime))
#ifdef MPI
            if (present(l_out_nc4par)) then
              if (l_out_nc4par) then
                nc_err = nf90_var_par_access(nc_id,time_id,nf90_collective)
                call ionc4_err(nc_err,ionc_rout,'nf90_var_par_access',trim(ionc_nomtime))
              end if
            end if
#endif
         endif
         
         if (kmax .GT. 0) then
            tab_dim_sigma(1) = dim_k_id
            nc_err = nf90_def_var(nc_id, trim(ionc_nomz), nf90_real, tab_dim_sigma, sig_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomz))
     !       tab_dim_sigma(1) = dim_w_id
     !       nc_err = nf90_def_var(nc_id, trim(ionc_pnomz_w), nf90_real, tab_dim_sigma, sig_w_id)
     !       call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomz_w))
         endif
         
! ******** On definit les attributs de la loupe  **********
         nc_err = nf90_put_att(nc_id, nf90_global, 'limi', limi)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','limi')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'lima', lima)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','lima')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'pasi', pasi)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','pasi')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'ljmi', ljmi)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','ljmi')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'ljma', ljma)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','ljma')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'pasj', pasj)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','pasi')
         
         if (kmax .GT. 0) then
            nc_err = nf90_put_att(nc_id, nf90_global, 'lkmi', lkmi)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','lkmi')
            
            nc_err = nf90_put_att(nc_id, nf90_global, 'lkma', lkma)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','lkma')
            
            nc_err = nf90_put_att(nc_id, nf90_global, 'pask', pask)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','pask')
            
         endif
         
! ******** On definit les attributs des variables  **********
         nc_err = nf90_put_att(nc_id,nj_id,'long_name','y-dimension of the grid')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:long_name')
         
         nc_err = nf90_put_att(nc_id,nj_id,'standard_name','y_grid_index')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:standard_name')
         
         nc_err = nf90_put_att(nc_id,nj_u_id,'long_name','y-dimension of the grid at u location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:long_name')
         
         nc_err = nf90_put_att(nc_id,nj_u_id,'standard_name','y_grid_index_at_u_location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:standard_name')
         
         nc_err = nf90_put_att(nc_id,nj_v_id,'long_name','y-dimension of the grid at v location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:long_name')
         
         nc_err = nf90_put_att(nc_id,nj_v_id,'standard_name','y_grid_index_at_v_location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:standard_name')
         
         nc_err = nf90_put_att(nc_id,nj_f_id,'long_name','y-dimension of the grid at f location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:long_name')
         
         nc_err = nf90_put_att(nc_id,nj_f_id,'standard_name','y_grid_index_at_f_location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:standard_name')
         
         !nc_err = nf90_put_att(nc_id, nj_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:units')

         !nc_err = nf90_put_att(nc_id, nj_u_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:units')

         !nc_err = nf90_put_att(nc_id, nj_v_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:units')
         
         !nc_err = nf90_put_att(nc_id, nj_f_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:units')
         
         nc_err = nf90_put_att(nc_id, nj_id,'axis','Y')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:axis')
         
         nc_err = nf90_put_att(nc_id, nj_u_id,'axis','Y')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:axis')
         
         nc_err = nf90_put_att(nc_id, nj_v_id,'axis','Y')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:axis')
         
         nc_err = nf90_put_att(nc_id, nj_f_id,'axis','Y')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:axis')
         
         nc_err = nf90_put_att(nc_id,ni_id,'long_name','x-dimension of the grid')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:long_name')
         
         nc_err = nf90_put_att(nc_id,ni_id,'standard_name','x_grid_index')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:standard_name')
         
         nc_err = nf90_put_att(nc_id,ni_u_id,'long_name','x-dimension of the grid at u location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:long_name')
         
         nc_err = nf90_put_att(nc_id,ni_u_id,'standard_name','x_grid_index_at_u_location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:standard_name')
         
         nc_err = nf90_put_att(nc_id,ni_v_id,'long_name','x-dimension of the grid at v location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:long_name')
         
         nc_err = nf90_put_att(nc_id,ni_v_id,'standard_name','x_grid_index_at_v_location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:standard_name')
         
         nc_err = nf90_put_att(nc_id,ni_f_id,'long_name','x-dimension of the grid at f location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:long_name')
         
         nc_err = nf90_put_att(nc_id,ni_f_id,'standard_name','x_grid_index_at_f_location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:standard_name')

         nc_err = nf90_put_att(nc_id,ni_id,'c_grid_axis_shift',0.)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:c_grid_axis_shift')

         nc_err = nf90_put_att(nc_id,ni_u_id,'c_grid_axis_shift',0.5)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI_U:c_grid_axis_shift')

         nc_err = nf90_put_att(nc_id,ni_v_id,'c_grid_axis_shift',0.)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI_V:c_grid_axis_shift')

         nc_err = nf90_put_att(nc_id,ni_f_id,'c_grid_axis_shift',0.5)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI_U:c_grid_axis_shift')

         nc_err = nf90_put_att(nc_id,nj_id,'c_grid_axis_shift',0.)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:c_grid_axis_shift')

         nc_err = nf90_put_att(nc_id,nj_u_id,'c_grid_axis_shift',0.)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ_U:c_grid_axis_shift')

         nc_err = nf90_put_att(nc_id,nj_v_id,'c_grid_axis_shift',0.5)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ_V:c_grid_axis_shift')

         nc_err = nf90_put_att(nc_id,nj_f_id,'c_grid_axis_shift',0.5)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ_V:c_grid_axis_shift')

         !nc_err = nf90_put_att(nc_id, ni_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:units')

         !nc_err = nf90_put_att(nc_id, ni_u_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:units')

         !nc_err = nf90_put_att(nc_id, ni_v_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:units')

         !nc_err = nf90_put_att(nc_id, ni_f_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:units')

         nc_err = nf90_put_att(nc_id, ni_id,'axis','X')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:axis')
         
         nc_err = nf90_put_att(nc_id, ni_u_id,'axis','X')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:axis')
         
         nc_err = nf90_put_att(nc_id, ni_v_id,'axis','X')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:axis')
         
         nc_err = nf90_put_att(nc_id, ni_f_id,'axis','X')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:axis')
         
         nc_err = nf90_put_att(nc_id,lat_id,'long_name','latitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:long_name')
         
         nc_err = nf90_put_att(nc_id,lat_id,'standard_name','latitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:standard_name')
         
         nc_err = nf90_put_att(nc_id, lat_id,'units', 'degrees_north')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:units')

         !nc_err = nf90_put_att(nc_id, lat_id,'content','yx')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:content')

         !nc_err = nf90_put_att(nc_id, lat_id,'coordinates','latitude longitude')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:coordinates')

         !nc_err = nf90_put_att(nc_id, lat_id,'associate','latitude longitude')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:associate')
         
         !nc_err = nf90_put_att(nc_id, lat_id,'axis','Y')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:axis')
         
         nc_err = nf90_put_att(nc_id,lat_id,'valid_min',REAL(-90.0,4))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',trim(ionc_nomlat))
         
         nc_err = nf90_put_att(nc_id,lat_id,'valid_max',REAL(90.0,4))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',trim(ionc_nomlat))
         
         nc_err = nf90_put_att(nc_id,lat_id,'_FillValue',REAL(999.0,4))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue','LAT:_FillValue')

         nc_err = nf90_put_att(nc_id,lon_id,'long_name','longitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:long_name')
         
         nc_err = nf90_put_att(nc_id,lon_id,'standard_name','longitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:standard_name')
         
         nc_err = nf90_put_att(nc_id, lon_id,'units', 'degrees_east')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:units')
         
         !nc_err = nf90_put_att(nc_id, lon_id,'axis','X')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:axis')
        
         !nc_err = nf90_put_att(nc_id, lon_id,'content','yx')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:content')

         !nc_err = nf90_put_att(nc_id, lon_id,'coordinates','latitude longitude')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:coordinates')

         !nc_err = nf90_put_att(nc_id, lon_id,'associate','latitude longitude')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:associate')
 
         nc_err = nf90_put_att(nc_id,lon_id,'valid_min',REAL(-180.0,4))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',trim(ionc_nomlon))
         
         nc_err = nf90_put_att(nc_id,lon_id,'valid_max',REAL(180.0,4))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',trim(ionc_nomlon))
         
         nc_err = nf90_put_att(nc_id,lon_id,'_FillValue',REAL(999.0,4))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue','LON:_FillValue')

         if (kmax .GT. 0) then
            
            !nc_err = nf90_put_att(nc_id,sig_id,'units','level')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:units')
            
            nc_err = nf90_put_att(nc_id,sig_id, 'long_name', 'sigma level')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:long_name')

            nc_err = nf90_put_att(nc_id,sig_id,'axis','Z')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')

            nc_err = nf90_put_att(nc_id,sig_id,'c_grid_axis_shift',0.)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:c_grid_axis_shift')

            !nc_err = nf90_put_att(nc_id,sig_id,'content','z')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:content')

            !nc_err = nf90_put_att(nc_id,sig_id,'coordinates','level')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:coordinates')

            !nc_err = nf90_put_att(nc_id,sig_id,'associate','level')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:associate')

       !     nc_err = nf90_put_att(nc_id,sig_w_id, 'long_name', 'sigma level at w location')
       !     call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:long_name')
            
            !nc_err = nf90_put_att(nc_id,sig_w_id,'standard_name','z_grid_index_at_w_location')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:standard_name')
            
       !     nc_err = nf90_put_att(nc_id,sig_w_id,'axis','Z')
       !     call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')
            
       !     nc_err = nf90_put_att(nc_id,sig_w_id,'c_grid_axis_shift',0.5)
       !     call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')

         endif
         
         if (dim_time .ne. -1) then
            nc_err = nf90_put_att(nc_id, time_id, 'long_name',ionc_longnamet)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:long_name')
            nc_err = nf90_put_att(nc_id, time_id, 'standard_name','time')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:standard_name')
            nc_err = nf90_put_att(nc_id, time_id, 'units',ionc_unitst)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:units')
            !nc_err = nf90_put_att(nc_id, time_id,'content','t')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:content')
            !nc_err = nf90_put_att(nc_id, time_id,'coordinates','time')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:coordinate')
            !nc_err = nf90_put_att(nc_id, time_id,'associate','time')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:associate')
            nc_err = nf90_put_att(nc_id, time_id,'axis','T')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:axis')
            nc_err = nf90_put_att(nc_id, time_id, 'time_origin',ionc_originet)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:time_origin')
         endif
         
! ******** On peut maintenant quitter le mode definition *****************
         
         nc_err = nf90_enddef(nc_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_enddef',' ')
         
! ******** Ajout des ni, nj ********************************************

         ini(:) = (/ (i,i=0,dim_lon-1) /)
         nc_err = nf90_put_var(nc_id,ni_id,ini)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlon))
         nc_err = nf90_put_var(nc_id,ni_v_id,ini)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlon_v))
         ini(:) = ini(:) + 0.5
         nc_err = nf90_put_var(nc_id,ni_u_id,ini)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlon_u))
         nc_err = nf90_put_var(nc_id,ni_f_id,ini)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlon_f))

         inj(:) = (/ (i,i=0,dim_lat-1) /)
         nc_err = nf90_put_var(nc_id,nj_id,inj)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlat))
         nc_err = nf90_put_var(nc_id,nj_u_id,inj)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlat_u))
         inj(:) = inj(:) + 0.5
         nc_err = nf90_put_var(nc_id,nj_v_id,inj)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlat_v))
         nc_err = nf90_put_var(nc_id,nj_f_id,inj)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlat_f))

! ******** Ajout des lon, lat ********************************************
         
         rlon(:,:) = lon2d(:,:)
         nc_err = nf90_put_var(nc_id,lon_id,rlon)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomlon))
         
         rlat(:,:) = lat2d(:,:)
         nc_err = nf90_put_var(nc_id,lat_id,rlat)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomlat))
         
! ******* Ajout de SIG(k)  ***************************************************
         
         if (kmax .GT. 0) then
           indice = 1
           if (present(zaxis_var)) then
             do ik = lkmi,lkma,pask
               rsig(indice) = zaxis_var(ik)
               indice = indice + 1
             enddo
           else
             do ik = lkmi,lkma,pask
               rsig(indice) = ik
               indice = indice + 1
             enddo
           endif
           nc_err = nf90_put_var(nc_id,sig_id,rsig)
           call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomz))
         endif
         
      endif
      
   else
      call ionc4_err(ionc_errcreation,ionc_rout, ' ' ,nom_fichier)
   endif
  
   deallocate(rlon,rlat)
   if(allocated(rsig)) deallocate(rsig)
 
   return
   end subroutine ionc4_createfile_2d_real

! ***************************************************************
! * SUBROUTINE ionc4_createfile_2d_double                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    20/02/01 par jfleroux                   *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)          *
!                       30/07/13 par A. Thevenin (CERFACS)      *
! *      *
! * Role : creer le fichier Netcdf qui contiendra les donnees      *
! *         Initialiser les donnees ne dependant pas du time :*
! *      - lon,lat*
! *      - SIG(k)*
! *         initialise les variables correspondants a la loupe*
! *         initialise les unites de temps et l'origine          *
! *        aux valeurs par defaut                               *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier : chemin complet du fichier a *
! *      creer*
! *      - dg   :               *
! *      - g0   :               *
! *      - dfi  :               *
! *      - fi0  :               *
! *      - limi : limite min de l'indice i de la loupe   *
! *      - lima : limite max de l'indice i de la loupe   *
! *      - pasi : pas de l'indice i de la loupe          *
! *      - ljmi : limite min de l'indice j de la loupe   *
! *      - ljma : limite max de l'indice j de la loupe*
! *      - pasj : pas de l'indice j de la loupe          *
! *      - lkmi : limite min de l'indice k de la loupe   *
! *      - lkma : limite max de l'indice k de la loupe   *
! *      - pask : pas de l'indice k de la loupe          *
! *      - kmax : borne max pour k                       *
! *      0 si 2D*
! *      - dim_time : nombre de pas de time que le *
! *           fichier contiendra si dim_time =0  *
! *           la dimension time est unlimited    *
! *           si dim_time = -1 pas de dim temps  *
! *      - SIG : niveaux*
! ***************************************************************
   subroutine ionc4_createfile_2d_double(nom_fichier,lon2d,lat2d, &
                                      limi,lima,pasi,ljmi,ljma,pasj,    &
                                      lkmi,lkma,pask,kmax,dim_time,     &
                                      l_out_nc4par,comm_active,zaxis_var)
#ifdef MPI 
   include 'mpif.h'
#endif

! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*), intent(in) :: nom_fichier
   integer,         intent(in)  :: dim_time
   integer,         intent(in)  :: limi,lima,pasi,ljmi,ljma,pasj,lkmi,lkma,pask,kmax
   integer,intent(in),optional  :: comm_active
   real(kind=8),dimension(INT((lima-limi)/pasi)+1,INT((ljma-ljmi)/pasj)+1),intent(in) :: lon2d,lat2d
   real(kind=4),dimension(kmax),intent(in),optional :: zaxis_var
   logical,intent(in),optional :: l_out_nc4par
   
! ******** VARIABLES DE TRAVAIL **************
   
   logical :: l_useless
   integer :: i,communicator
   integer :: dim_lon,dim_lat,dim_k
   
   integer :: nc_id,lat_id,lon_id,time_id,sig_id,ni_id,nj_id, &
              ni_u_id,nj_u_id,ni_v_id,nj_v_id,ni_f_id,nj_f_id,sig_w_id
   
! ** Declaration des IDS identifiant les dimensions des variables.
   integer :: dim_lon_id,dim_lat_id,dim_time_id,dim_k_id,dim_w_id, &
              dim_lon_u_id,dim_lat_u_id,dim_lon_v_id,dim_lat_v_id,dim_lon_f_id,dim_lat_f_id
   
! ** Declarations des tableaux qui contiendront les IDS des dimensions des variables
   integer,dimension(1) :: tab_dim_time,tab_dim_sigma
   integer,dimension(2) :: tab_dim_lon,tab_dim_lat
   integer  :: nc_err,ik,indice,cmode,cache_size
   real(kind=4),dimension(:),allocatable :: ini,inj
   real(kind=8),dimension(:,:),allocatable :: rlon,rlat
   real(kind=4),dimension(:),  allocatable :: rsig

! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   
   ionc_rout = "ionc4_createfile_2d_double"
   
! *** Recuperation des dimensions
   dim_lon = INT((lima-limi)/pasi) +1
   dim_lat = INT((ljma-ljmi)/pasj) +1
   allocate(ini(dim_lon))
   allocate(inj(dim_lat))
   allocate(rlon(dim_lon,dim_lat))
   allocate(rlat(dim_lon,dim_lat))
   cache_size = 10*dim_lon*dim_lat
    
   if (kmax .GT. 0) then
      dim_k   = (lkma-lkmi)/pask +1
      allocate(rsig(dim_k))
      cache_size = cache_size*dim_k
   end if
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      
      ionc_nfich = ionc_nfich + 1
      if (ionc_nfich .gt. ionc_longtabfich) then
         call ionc4_err(ionc_errcreation,ionc_rout, ' ' ,nom_fichier)
         stop
      else
         ionc_nomfich(ionc_nfich) = nom_fichier
         cmode = nf90_share
#ifdef MPI
         if (present(comm_active)) then
            communicator=comm_active
         else
#ifdef key_oasis
            communicator=MPI_COMM_MARS
#else
            communicator=MPI_COMM_WORLD
#endif
         endif
         if (present(l_out_nc4par)) then
            if (l_out_nc4par) then
               !print*,'Using parallel I/O features !!!' 
               cmode = ior(NF90_NETCDF4, NF90_CLASSIC_MODEL)
               cmode = ior(cmode, NF90_MPIIO)
               nc_err = nf90_create(nom_fichier, cmode, nc_id, cache_size=cache_size, &
                                    comm=communicator, info=MPI_INFO_NULL)
            else
               nc_err = nf90_create(nom_fichier, cmode, nc_id)
            endif
         else
            nc_err = nf90_create(nom_fichier, cmode, nc_id)
         endif
#else
         if (present(comm_active)) communicator=comm_active  ! useless for portability only
         if (present(l_out_nc4par)) l_useless = l_out_nc4par  ! useless for portability only
         nc_err = nf90_create(nom_fichier, cmode, nc_id)
#endif
         call ionc4_err(nc_err,ionc_rout, 'nf90_create',nom_fichier)
         
         ionc_idfich(ionc_nfich) = nc_id
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlon),dim_lon,dim_lon_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_pnomlon)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlat),dim_lat,dim_lat_id)
         call ionc4_err(nc_err,ionc_rout, 'nf90_def_dim',ionc_pnomlat)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlon_u),dim_lon,dim_lon_u_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_pnomlon_u)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlat_u),dim_lat,dim_lat_u_id)
         call ionc4_err(nc_err,ionc_rout, 'nf90_def_dim',ionc_pnomlat_u)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlon_v),dim_lon,dim_lon_v_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_pnomlon_v)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlat_v),dim_lat,dim_lat_v_id)
         call ionc4_err(nc_err,ionc_rout, 'nf90_def_dim',ionc_pnomlat_v)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlon_f),dim_lon,dim_lon_f_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_pnomlon_f)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_pnomlat_f),dim_lat,dim_lat_f_id)
         call ionc4_err(nc_err,ionc_rout, 'nf90_def_dim',ionc_pnomlat_f)
         
         
         if (kmax .GT. 0) then
            nc_err = nf90_def_dim(nc_id,trim(ionc_nomz),dim_k,dim_k_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_nomz)
            !nc_err = nf90_def_dim(nc_id,trim(ionc_pnomz_w),dim_k,dim_w_id)
            !call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_pnomz_w)
         end if
         
         if (dim_time .ne. -1) then
            nc_err = nf90_def_dim(nc_id,trim(ionc_nomtime),dim_time,dim_time_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_nomtime)
         endif
         
! ******** On definit les variables  **********
         
         tab_dim_lon(1) = dim_lon_id
         tab_dim_lat(1) = dim_lat_id
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlon), nf90_real, tab_dim_lon(1), ni_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlon))
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlat), nf90_real, tab_dim_lat(1), nj_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlat))
         
         tab_dim_lon(1) = dim_lon_u_id
         tab_dim_lat(1) = dim_lat_u_id
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlon_u), nf90_real, tab_dim_lon(1), ni_u_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlon_u))
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlat_u), nf90_real, tab_dim_lat(1), nj_u_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlat_u))
         
         tab_dim_lon(1) = dim_lon_v_id
         tab_dim_lat(1) = dim_lat_v_id
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlon_v), nf90_real, tab_dim_lon(1), ni_v_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlon_v))
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlat_v), nf90_real, tab_dim_lat(1), nj_v_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlat_v))
         
         tab_dim_lon(1) = dim_lon_f_id
         tab_dim_lat(1) = dim_lat_f_id
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlon_f), nf90_real, tab_dim_lon(1), ni_f_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlon_f))
         nc_err = nf90_def_var(nc_id, trim(ionc_pnomlat_f), nf90_real, tab_dim_lat(1), nj_f_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomlat_f))
         
         tab_dim_lat(1) = dim_lon_id
         tab_dim_lat(2) = dim_lat_id
         nc_err = nf90_def_var(nc_id, trim(ionc_nomlat), nf90_double, tab_dim_lat, lat_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomlat))
         
         tab_dim_lon(1) = dim_lon_id
         tab_dim_lon(2) = dim_lat_id
         nc_err = nf90_def_var(nc_id, trim(ionc_nomlon), nf90_double, tab_dim_lon, lon_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomlon))
         
         if (dim_time .ne. -1) then
            tab_dim_time(1) = dim_time_id
            nc_err = nf90_def_var(nc_id, trim(ionc_nomtime), nf90_double, tab_dim_time, time_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomtime))
#ifdef MPI
            if (present(l_out_nc4par)) then
              if (l_out_nc4par) then
                nc_err = nf90_var_par_access(nc_id,time_id,nf90_collective)
                call ionc4_err(nc_err,ionc_rout,'nf90_var_par_access',trim(ionc_nomtime))
              end if
            end if
#endif
         endif
         
         if (kmax .GT. 0) then
            tab_dim_sigma(1) = dim_k_id
            nc_err = nf90_def_var(nc_id, trim(ionc_nomz), nf90_real, tab_dim_sigma, sig_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomz))
            !tab_dim_sigma(1) = dim_w_id
            !nc_err = nf90_def_var(nc_id, trim(ionc_pnomz_w), nf90_real, tab_dim_sigma, sig_w_id)
            !call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_pnomz_w))
         endif
         
! ******** On definit les attributs de la loupe  **********
         nc_err = nf90_put_att(nc_id, nf90_global, 'limi', limi)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','limi')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'lima', lima)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','lima')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'pasi', pasi)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','pasi')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'ljmi', ljmi)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','ljmi')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'ljma', ljma)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','ljma')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'pasj', pasj)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','pasi')
         
         if (kmax .GT. 0) then
            nc_err = nf90_put_att(nc_id, nf90_global, 'lkmi', lkmi)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','lkmi')
            
            nc_err = nf90_put_att(nc_id, nf90_global, 'lkma', lkma)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','lkma')
            
            nc_err = nf90_put_att(nc_id, nf90_global, 'pask', pask)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','pask')
            
         endif
         
! ******** On definit les attributs des variables  **********
         nc_err = nf90_put_att(nc_id,nj_id,'long_name','y-dimension of the grid')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:long_name')
         
         nc_err = nf90_put_att(nc_id,nj_id,'standard_name','y_grid_index')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:standard_name')
         
         nc_err = nf90_put_att(nc_id,nj_u_id,'long_name','y-dimension of the grid at u location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:long_name')
         
         nc_err = nf90_put_att(nc_id,nj_u_id,'standard_name','y_grid_index_at_u_location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:standard_name')
         
         nc_err = nf90_put_att(nc_id,nj_v_id,'long_name','y-dimension of the grid at v location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:long_name')
         
         nc_err = nf90_put_att(nc_id,nj_v_id,'standard_name','y_grid_index_at_v_location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:standard_name')
         
         nc_err = nf90_put_att(nc_id,nj_f_id,'long_name','y-dimension of the grid at f location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:long_name')
         
         nc_err = nf90_put_att(nc_id,nj_f_id,'standard_name','y_grid_index_at_f_location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:standard_name')
         
         !nc_err = nf90_put_att(nc_id, nj_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:units')

         !nc_err = nf90_put_att(nc_id, nj_u_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:units')

         !nc_err = nf90_put_att(nc_id, nj_v_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:units')

         !nc_err = nf90_put_att(nc_id, nj_f_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:units')

         nc_err = nf90_put_att(nc_id, nj_id,'axis','Y')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:axis')
         
         nc_err = nf90_put_att(nc_id, nj_u_id,'axis','Y')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:axis')
         
         nc_err = nf90_put_att(nc_id, nj_v_id,'axis','Y')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:axis')
         
         nc_err = nf90_put_att(nc_id, nj_f_id,'axis','Y')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:axis')
         
         nc_err = nf90_put_att(nc_id,ni_id,'long_name','x-dimension of the grid')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:long_name')
         
         nc_err = nf90_put_att(nc_id,ni_id,'standard_name','x_grid_index')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:standard_name')
         
         nc_err = nf90_put_att(nc_id,ni_u_id,'long_name','x-dimension of the grid at u location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:long_name')
         
         nc_err = nf90_put_att(nc_id,ni_u_id,'standard_name','x_grid_index_at_u_location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:standard_name')
         
         nc_err = nf90_put_att(nc_id,ni_v_id,'long_name','x-dimension of the grid at v location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:long_name')
         
         nc_err = nf90_put_att(nc_id,ni_v_id,'standard_name','x_grid_index_at_v_location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:standard_name')
         
         nc_err = nf90_put_att(nc_id,ni_f_id,'long_name','x-dimension of the grid at f location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:long_name')
         
         nc_err = nf90_put_att(nc_id,ni_f_id,'standard_name','x_grid_index_at_f_location')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:standard_name')
         
         nc_err = nf90_put_att(nc_id,ni_id,'c_grid_axis_shift',0.)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:c_grid_axis_shift')

         nc_err = nf90_put_att(nc_id,ni_u_id,'c_grid_axis_shift',0.5)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI_U:c_grid_axis_shift')

         nc_err = nf90_put_att(nc_id,ni_v_id,'c_grid_axis_shift',0.)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI_V:c_grid_axis_shift')
         
         nc_err = nf90_put_att(nc_id,ni_f_id,'c_grid_axis_shift',0.5)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI_U:c_grid_axis_shift')

         nc_err = nf90_put_att(nc_id,nj_id,'c_grid_axis_shift',0.)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ:c_grid_axis_shift')

         nc_err = nf90_put_att(nc_id,nj_u_id,'c_grid_axis_shift',0.)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ_U:c_grid_axis_shift')

         nc_err = nf90_put_att(nc_id,nj_v_id,'c_grid_axis_shift',0.5)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ_V:c_grid_axis_shift')
        
         nc_err = nf90_put_att(nc_id,nj_f_id,'c_grid_axis_shift',0.5)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NJ_V:c_grid_axis_shift')
        
         !nc_err = nf90_put_att(nc_id, ni_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:units')

         !nc_err = nf90_put_att(nc_id, ni_u_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:units')

         !nc_err = nf90_put_att(nc_id, ni_v_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:units')

         !nc_err = nf90_put_att(nc_id, ni_f_id,'units', '1')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:units')

         nc_err = nf90_put_att(nc_id, ni_id,'axis','X')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:axis')
         
         nc_err = nf90_put_att(nc_id, ni_u_id,'axis','X')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:axis')
         
         nc_err = nf90_put_att(nc_id, ni_v_id,'axis','X')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:axis')
         
         nc_err = nf90_put_att(nc_id, ni_f_id,'axis','X')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','NI:axis')
         
         nc_err = nf90_put_att(nc_id,lat_id,'long_name','latitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:long_name')
         
         nc_err = nf90_put_att(nc_id,lat_id,'standard_name','latitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:standard_name')
         
         nc_err = nf90_put_att(nc_id, lat_id,'units', 'degrees_north')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:units')
         
         !nc_err = nf90_put_att(nc_id, lat_id,'axis','Y')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:axis')

         !nc_err = nf90_put_att(nc_id, lat_id,'content','yx')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:content')

         !nc_err = nf90_put_att(nc_id, lat_id,'coordinates','latitude longitude')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:coordinates')

         !nc_err = nf90_put_att(nc_id, lat_id,'associate','latitude longitude')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:associate')

         nc_err = nf90_put_att(nc_id,lat_id,'valid_min',REAL(-90.0,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',trim(ionc_nomlat))
         
         nc_err = nf90_put_att(nc_id,lat_id,'valid_max',REAL(90.0,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',trim(ionc_nomlat))
         
         nc_err = nf90_put_att(nc_id,lat_id,'_FillValue',REAL(1.7e+38,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue','LAT:_FillValue')

         nc_err = nf90_put_att(nc_id,lon_id,'long_name','longitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:long_name')
         
         nc_err = nf90_put_att(nc_id,lon_id,'standard_name','longitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:standard_name')
         
         nc_err = nf90_put_att(nc_id, lon_id,'units', 'degrees_east')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:units')
         
         !nc_err = nf90_put_att(nc_id, lon_id,'axis','X')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:axis')
        
         !nc_err = nf90_put_att(nc_id, lon_id,'content','yx')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:content')

         !nc_err = nf90_put_att(nc_id, lon_id,'coordinates','latitude longitude')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:coordinates')

         !nc_err = nf90_put_att(nc_id, lon_id,'associate','latitude longitude')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:associate')
         
         nc_err = nf90_put_att(nc_id,lon_id,'valid_min',REAL(-180.0,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',trim(ionc_nomlon))
         
         nc_err = nf90_put_att(nc_id,lon_id,'valid_max',REAL(180.0,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att_double:valid_max',trim(ionc_nomlon))
         
         nc_err = nf90_put_att(nc_id,lon_id,'_FillValue',REAL(1.7e+38,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue','LON:_FillValue')

         if (kmax .GT. 0) then
            
            !nc_err = nf90_put_att(nc_id,sig_id,'units','level')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:units')
            
            !nc_err = nf90_put_att(nc_id,sig_id,'standard_name','z_grid_index')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:standard_name')

            !nc_err = nf90_put_att(nc_id,sig_id,'content','z')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:content')

            !nc_err = nf90_put_att(nc_id,sig_id,'coordinates','level')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:coordinates')

            !nc_err = nf90_put_att(nc_id,sig_id,'associate','level')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:associate')
            
            nc_err = nf90_put_att(nc_id,sig_id, 'long_name', 'sigma level')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')

            nc_err = nf90_put_att(nc_id,sig_id,'axis','Z')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')

            nc_err = nf90_put_att(nc_id,sig_id,'c_grid_axis_shift',0.)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')

            !nc_err = nf90_put_att(nc_id,sig_w_id, 'long_name', 'sigma level at w location')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:long_name')

            !nc_err = nf90_put_att(nc_id,sig_w_id,'axis','Z')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')

            !nc_err = nf90_put_att(nc_id,sig_w_id,'c_grid_axis_shift',0.5)
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')

         endif
         
         if (dim_time .ne. -1) then
            nc_err = nf90_put_att(nc_id, time_id, 'long_name',ionc_longnamet)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:long_name')
            nc_err = nf90_put_att(nc_id, time_id, 'standard_name','time')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:standard_name')
            nc_err = nf90_put_att(nc_id, time_id, 'units',ionc_unitst)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:units')
            !nc_err = nf90_put_att(nc_id, time_id,'content','t')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:content')
            !nc_err = nf90_put_att(nc_id, time_id,'coordinates','time')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:coordinate')
            !nc_err = nf90_put_att(nc_id, time_id,'associate','time')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:associate')
            nc_err = nf90_put_att(nc_id, time_id,'axis','T')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:axis')
            nc_err = nf90_put_att(nc_id, time_id, 'time_origin',ionc_originet)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:time_origin')
         endif
         
! ******** On peut maintenant quitter le mode definition *****************
         nc_err = nf90_enddef(nc_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_enddef',' ')
         
! ******** Ajout des ni, nj ********************************************
         
         ini(:) = (/ (i,i=0,dim_lon-1) /)
         nc_err = nf90_put_var(nc_id,ni_id,ini)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlon))
         nc_err = nf90_put_var(nc_id,ni_v_id,ini)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlon_v))
         ini(:) = ini(:) + 0.5
         nc_err = nf90_put_var(nc_id,ni_u_id,ini)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlon_u))
         nc_err = nf90_put_var(nc_id,ni_f_id,ini)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlon_f))

         inj(:) = (/ (i,i=0,dim_lat-1) /)
         nc_err = nf90_put_var(nc_id,nj_id,inj)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlat))
         nc_err = nf90_put_var(nc_id,nj_u_id,inj)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlat_u))
         inj(:) = inj(:) + 0.5
         nc_err = nf90_put_var(nc_id,nj_v_id,inj)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlat_v))
         nc_err = nf90_put_var(nc_id,nj_f_id,inj)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_pnomlat_f))
         
! ******** Ajout des lon, lat ********************************************
         
         rlon(:,:) = lon2d(:,:)
         nc_err = nf90_put_var(nc_id,lon_id,rlon)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomlon))
         
         rlat(:,:) = lat2d(:,:)
         nc_err = nf90_put_var(nc_id,lat_id,rlat)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var_double',trim(ionc_nomlat))
         
! ******* Ajout de SIG(k)  ***************************************************
         
         if (kmax .GT. 0) then
           indice = 1
           if (present(zaxis_var)) then
             do ik = lkmi,lkma,pask
               rsig(indice) = zaxis_var(ik)
               indice = indice + 1
             enddo
           else
             do ik = lkmi,lkma,pask
               rsig(indice) = ik
               indice = indice + 1
             enddo
           endif
           nc_err = nf90_put_var(nc_id,sig_id,rsig)
           call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomz))
         endif
         
      endif
      
   else
      call ionc4_err(ionc_errcreation,ionc_rout, ' ' ,nom_fichier)
   endif
  
   deallocate(rlon,rlat)
   if(allocated(rsig)) deallocate(rsig)
 
   return
   end subroutine ionc4_createfile_2d_double
   
! ***************************************************************
! * SUBROUTINE ionc4_createfile_1d_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    20/02/01 par jfleroux                   *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)          *
!                       30/07/13 par A. Thevenin (CERFACS)      *
! *      *
! * Role : creer le fichier Netcdf qui contiendra les donnees      *
! *         Initialiser les donnees ne dependant pas du time :*
! *      - lon,lat*
! *      - SIG(k)*
! *         initialise les variables correspondants a la loupe*
! *         initialise les unites de temps et l'origine          *
! *        aux valeurs par defaut                               *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier : chemin complet du fichier a *
! *      creer*
! *      - dg   :               *
! *      - g0   :               *
! *      - dfi  :               *
! *      - fi0  :               *
! *      - limi : limite min de l'indice i de la loupe   *
! *      - lima : limite max de l'indice i de la loupe   *
! *      - pasi : pas de l'indice i de la loupe          *
! *      - ljmi : limite min de l'indice j de la loupe   *
! *      - ljma : limite max de l'indice j de la loupe*
! *      - pasj : pas de l'indice j de la loupe          *
! *      - lkmi : limite min de l'indice k de la loupe   *
! *      - lkma : limite max de l'indice k de la loupe   *
! *      - pask : pas de l'indice k de la loupe          *
! *      - kmax : borne max pour k                       *
! *      0 si 2D*
! *      - dim_time : nombre de pas de time que le *
! *           fichier contiendra si dim_time =0  *
! *           la dimension time est unlimited    *
! *           si dim_time = -1 pas de dim temps  *
! *      - SIG : niveaux           *
! ***************************************************************
   subroutine ionc4_createfile_1d_real(nom_fichier,dg,g0,dfi,fi0,  &
                                    limi,lima,pasi,ljmi,ljma,pasj, &
                                    lkmi,lkma,pask,kmax,dim_time,  &
                                    l_out_nc4par)
#ifdef MPI   
   include 'mpif.h'
#endif 
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier
   integer  :: dim_time
   integer  :: limi,lima,pasi,ljmi,ljma,pasj,lkmi,lkma,pask,kmax
   real(kind=4)     :: dg,g0,dfi,fi0
   logical,intent(in),optional :: l_out_nc4par
   
! ******** VARIABLES DE TRAVAIL **************
   
   logical :: l_useless
   integer :: dim_lon,dim_lat,dim_k
   
   integer :: nc_id,lat_id,lon_id,time_id,sig_id
   
! ** Declaration des IDS identifiant les dimensions des variables.
   integer :: dim_lon_id,dim_lat_id,dim_time_id,dim_k_id
   
! ** Declarations des tableaux qui contiendront les IDS des dimensions des variables
   integer,dimension(1) :: tab_dim_lon,tab_dim_lat,tab_dim_time,tab_dim_sigma
   integer  :: nc_err,ilon,ilat,ik,indice,cmode,cache_size
   real(kind=4),dimension(:),allocatable :: rlon,rlat,rsig
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_createfile_1d_real"
   
! *** Recuperation des dimensions
   dim_lon = (lima-limi)/pasi +1
   allocate(rlon(dim_lon))
   dim_lat = (ljma-ljmi)/pasj +1
   allocate(rlat(dim_lat))
   cache_size = 10*dim_lon*dim_lat 
   if (kmax .GT. 0) then
      dim_k   = (lkma-lkmi)/pask +1
      allocate(rsig(dim_k))
      cache_size = cache_size*dim_k
   endif

   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      
      ionc_nfich = ionc_nfich + 1
      if (ionc_nfich .gt. ionc_longtabfich) then
         call ionc4_err(ionc_errcreation,ionc_rout, ' ' ,nom_fichier)
         stop
      else
         ionc_nomfich(ionc_nfich) = nom_fichier
         cmode = nf90_share
#ifdef MPI
         if (present(l_out_nc4par)) then
            if (l_out_nc4par) then
               !print*,'Using parallel I/O features !!!' 
               cmode = ior(nf90_netcdf4,nf90_classic_model)
               cmode = ior(cmode,nf90_mpiio)
#ifdef key_oasis
               nc_err = nf90_create(nom_fichier, cmode, nc_id, cache_size=cache_size,comm=MPI_COMM_MARS, info=MPI_INFO_NULL)  
#else
               nc_err = nf90_create(nom_fichier, cmode, nc_id, cache_size=cache_size,comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)
#endif
            else
               nc_err = nf90_create(nom_fichier, cmode, nc_id)
            endif
         else
            nc_err = nf90_create(nom_fichier, cmode, nc_id)        
         endif
#else
         if (present(l_out_nc4par)) l_useless = l_out_nc4par  ! useless for portability only
         nc_err = nf90_create(nom_fichier, cmode, nc_id)
#endif
         
         call ionc4_err(nc_err,ionc_rout, 'nf90_create',nom_fichier)
         
         ionc_idfich(ionc_nfich) = nc_id
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_nomlon),dim_lon,dim_lon_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_nomlon)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_nomlat),dim_lat,dim_lat_id)
         call ionc4_err(nc_err,ionc_rout, 'nf90_def_dim',ionc_nomlat)
         
         
         if (kmax .GT. 0) then
            nc_err = nf90_def_dim(nc_id,trim(ionc_nomz),dim_k,dim_k_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_nomz)
         endif
         
         if (dim_time .ne. -1) then
            nc_err = nf90_def_dim(nc_id,trim(ionc_nomtime),dim_time,dim_time_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_nomtime)
         endif
         
         
! ******** On definit les variables  **********
         
         tab_dim_lat(1) = dim_lat_id
         nc_err = nf90_def_var(nc_id, trim(ionc_nomlat), nf90_real, tab_dim_lat, lat_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomlat))
         
         tab_dim_lon(1) = dim_lon_id
         nc_err = nf90_def_var(nc_id, trim(ionc_nomlon), nf90_real, tab_dim_lon, lon_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomlon))
         
         if (dim_time .ne. -1) then
            tab_dim_time(1) = dim_time_id
            nc_err = nf90_def_var(nc_id, trim(ionc_nomtime), nf90_real, tab_dim_time, time_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomtime))
#ifdef MPI
            if (present(l_out_nc4par)) then
              if (l_out_nc4par) then
                nc_err = nf90_var_par_access(nc_id,time_id,nf90_collective)
                call ionc4_err(nc_err,ionc_rout,'nf90_var_par_access',trim(ionc_nomtime))
              end if
            end if
#endif
         endif
         
         if (kmax .GT. 0) then
            tab_dim_sigma(1) = dim_k_id
            nc_err = nf90_def_var(nc_id, trim(ionc_nomz), nf90_real, tab_dim_sigma, sig_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomz))
         endif
         
! ******** On definit les attributs de la loupe  **********
         nc_err = nf90_put_att(nc_id, nf90_global, 'limi', limi)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','limi')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'lima', lima)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','lima')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'pasi', pasi)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','pasi')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'ljmi', ljmi)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','ljmi')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'ljma', ljma)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','ljma')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'pasj', pasj)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','pasi')
         
         if (kmax .GT. 0) then
            nc_err = nf90_put_att(nc_id, nf90_global, 'lkmi', lkmi)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','lkmi')
            
            nc_err = nf90_put_att(nc_id, nf90_global, 'lkma', lkma)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','lkma')
            
            nc_err = nf90_put_att(nc_id, nf90_global, 'pask', pask)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','pask')
            
         endif
         
! ******** On definit les attributs des variables  **********
         nc_err = nf90_put_att(nc_id,lat_id,'long_name','latitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:long_name')
         
         nc_err = nf90_put_att(nc_id,lat_id,'standard_name','latitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:standard_name')
         
         nc_err = nf90_put_att(nc_id, lat_id,'units', 'degrees_north')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:units')
         
         nc_err = nf90_put_att(nc_id, lat_id,'axis','Y')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:axis')
         
         nc_err = nf90_put_att(nc_id,lat_id,'valid_min',REAL(-90.0,4))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',trim(ionc_nomlat))
         
         nc_err = nf90_put_att(nc_id,lat_id,'valid_max',REAL(90.0,4))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',trim(ionc_nomlat))
         
         nc_err = nf90_put_att(nc_id,lat_id,'_FillValue',REAL(999.0,4))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue','LAT:_FillValue')

         nc_err = nf90_put_att(nc_id,lon_id,'long_name','longitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:long_name')
         
         nc_err = nf90_put_att(nc_id,lon_id,'standard_name','longitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:standard_name')
         
         nc_err = nf90_put_att(nc_id, lon_id,'units', 'degrees_east')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:units')
         
         nc_err = nf90_put_att(nc_id, lon_id,'axis','X')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:axis')
         
         nc_err = nf90_put_att(nc_id,lon_id,'valid_min',REAL(-180.0,4))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',trim(ionc_nomlon))
         
         nc_err = nf90_put_att(nc_id,lon_id,'valid_max',REAL(180.0,4))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',trim(ionc_nomlon))
         
         nc_err = nf90_put_att(nc_id,lon_id,'_FillValue',REAL(999.0,4))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue','LON:_FillValue')

         if (kmax .GT. 0) then
            
            !nc_err = nf90_put_att(nc_id,sig_id,'units','level')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:units')
            
            nc_err = nf90_put_att(nc_id,sig_id, 'long_name', 'sigma level')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:long_name')

            nc_err = nf90_put_att(nc_id,sig_id,'axis','Z')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')
            
            nc_err = nf90_put_att(nc_id,sig_id,'c_grid_axis_shift',0.)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')

            !nc_err = nf90_put_att(nc_id,sig_id,'standard_name','z_grid_index')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:standard_name')

            !nc_err = nf90_put_att(nc_id,sig_id,'content','z')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:content')

            !nc_err = nf90_put_att(nc_id,sig_id,'coordinates','level')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:coordinates')

            !nc_err = nf90_put_att(nc_id,sig_id,'associate','level')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:associate')
            
         endif
         
         if (dim_time .ne. -1) then
            nc_err = nf90_put_att(nc_id, time_id, 'long_name',ionc_longnamet)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:long_name')
            nc_err = nf90_put_att(nc_id, time_id, 'standard_name','time')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:standard_name')
            nc_err = nf90_put_att(nc_id, time_id, 'units',ionc_unitst)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:units')
            !nc_err = nf90_put_att(nc_id, time_id,'content','t')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:content')
            !nc_err = nf90_put_att(nc_id, time_id,'coordinates','time')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:coordinate')
            !nc_err = nf90_put_att(nc_id, time_id,'associate','time')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:associate')
            nc_err = nf90_put_att(nc_id, time_id,'axis','T')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:axis')
            nc_err = nf90_put_att(nc_id, time_id, 'time_origin',ionc_originet)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:time_origin')
         endif
         
! ******** On peut maintenant quitter le mode definition *****************
         
         nc_err = nf90_enddef(nc_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_enddef',' ')
         
! ******** Ajout des lon, lat ********************************************
         
         indice = 1
         do ilon = limi,lima,pasi
            rlon(indice) = real(ilon-1)*dg + g0
            indice = indice + 1
         enddo
         nc_err = nf90_put_var(nc_id,lon_id,rlon)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomlon))
         
         indice = 1
         do ilat = ljmi,ljma,pasj
            rlat(indice) = real(ilat-1)*dfi + fi0
            indice = indice + 1
         enddo
         nc_err = nf90_put_var(nc_id,lat_id,rlat)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomlat))
         
! ******* Ajout de SIG(k)  ***************************************************
         
         if (kmax .GT. 0) then
            indice = 1
            do ik = lkmi,lkma,pask
               !rsig(indice) = sig(ik)
               rsig(indice) = ik
               indice = indice + 1
            enddo
            nc_err = nf90_put_var(nc_id,sig_id,rsig)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomz))
         endif
         
      endif
      
   else
      call ionc4_err(ionc_errcreation,ionc_rout, ' ' ,nom_fichier)
   endif
  
   deallocate(rlon,rlat)
   if(allocated(rsig)) deallocate(rsig)
 
   return
   end subroutine ionc4_createfile_1d_real

! ***************************************************************
! * SUBROUTINE ionc4_createfile_1d_double                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    20/02/01 par jfleroux                   *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)          *
!                       30/07/13 par A. Thevenin (CERFACS)      *
! *      *
! * Role : creer le fichier Netcdf qui contiendra les donnees      *
! *         Initialiser les donnees ne dependant pas du time :*
! *      - lon,lat*
! *      - SIG(k)*
! *         initialise les variables correspondants a la loupe*
! *         initialise les unites de temps et l'origine          *
! *        aux valeurs par defaut                               *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier : chemin complet du fichier a *
! *      creer*
! *      - dg   :               *
! *      - g0   :               *
! *      - dfi  :               *
! *      - fi0  :               *
! *      - limi : limite min de l'indice i de la loupe   *
! *      - lima : limite max de l'indice i de la loupe   *
! *      - pasi : pas de l'indice i de la loupe          *
! *      - ljmi : limite min de l'indice j de la loupe   *
! *      - ljma : limite max de l'indice j de la loupe*
! *      - pasj : pas de l'indice j de la loupe          *
! *      - lkmi : limite min de l'indice k de la loupe   *
! *      - lkma : limite max de l'indice k de la loupe   *
! *      - pask : pas de l'indice k de la loupe          *
! *      - kmax : borne max pour k                       *
! *      0 si 2D*
! *      - dim_time : nombre de pas de time que le *
! *           fichier contiendra si dim_time =0  *
! *           la dimension time est unlimited    *
! *           si dim_time = -1 pas de dim temps  *
! *      - SIG : niveaux*
! ***************************************************************
   subroutine ionc4_createfile_1d_double(nom_fichier,dg,g0,dfi,fi0,  &
                                      limi,lima,pasi,ljmi,ljma,pasj, &
                                      lkmi,lkma,pask,kmax,dim_time,  &
                                      l_out_nc4par)
#ifdef MPI 
   include 'mpif.h'
#endif

! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier
   integer :: dim_time
   integer :: limi,lima,pasi,ljmi,ljma,pasj,lkmi,lkma,pask,kmax
   real(kind=8) :: dg,g0,dfi,fi0
   logical,intent(in),optional :: l_out_nc4par
   
! ******** VARIABLES DE TRAVAIL **************
   
   logical :: l_useless
   integer :: dim_lon,dim_lat,dim_k
   
   integer :: nc_id,lat_id,lon_id,time_id,sig_id
   
! ** Declaration des IDS identifiant les dimensions des variables.
   integer :: dim_lon_id,dim_lat_id,dim_time_id,dim_k_id
   
! ** Declarations des tableaux qui contiendront les IDS des dimensions des variables
   integer,dimension(1) :: tab_dim_lon,tab_dim_lat,tab_dim_time,tab_dim_sigma
   integer  :: nc_err,ilon,ilat,ik,indice,cmode,cache_size
   real(kind=8),dimension(:),allocatable :: rlon,rlat
   real(kind=4),dimension(:),allocatable :: rsig
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   
   ionc_rout = "ionc4_createfile_1d_double"
   
! *** Recuperation des dimensions
   dim_lon = (lima-limi)/pasi +1
   allocate(rlon(dim_lon))
   dim_lat = (ljma-ljmi)/pasj +1
   allocate(rlat(dim_lat))
   cache_size = 10*dim_lon*dim_lat
    
   if (kmax .GT. 0) then
      dim_k   = (lkma-lkmi)/pask +1
      allocate(rsig(dim_k))
      cache_size = cache_size*dim_k
   end if
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      
      ionc_nfich = ionc_nfich + 1
      if (ionc_nfich .gt. ionc_longtabfich) then
         call ionc4_err(ionc_errcreation,ionc_rout, ' ' ,nom_fichier)
         stop
      else
         ionc_nomfich(ionc_nfich) = nom_fichier
         cmode = nf90_share
#ifdef MPI
         if (present(l_out_nc4par)) then
            if (l_out_nc4par) then
               !print*,'Using parallel I/O features !!!' 
               cmode = ior(NF90_NETCDF4,NF90_CLASSIC_MODEL)
               cmode = ior(cmode,NF90_MPIIO)
#ifdef key_oasis
               nc_err = nf90_create(nom_fichier, cmode, nc_id, cache_size = cache_size,comm=MPI_COMM_MARS, info=MPI_INFO_NULL)
#else
               nc_err = nf90_create(nom_fichier, cmode, nc_id, cache_size = cache_size,comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)
#endif
            else
               nc_err = nf90_create(nom_fichier, cmode, nc_id)
            endif
         else
            nc_err = nf90_create(nom_fichier, cmode, nc_id)
         endif
#else
         if (present(l_out_nc4par)) l_useless = l_out_nc4par  ! useless for portability only
         nc_err = nf90_create(nom_fichier, cmode, nc_id)
#endif
         call ionc4_err(nc_err,ionc_rout, 'nf90_create',nom_fichier)
         
         ionc_idfich(ionc_nfich) = nc_id
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_nomlon),dim_lon,dim_lon_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_nomlon)
         
         nc_err = nf90_def_dim(nc_id,trim(ionc_nomlat),dim_lat,dim_lat_id)
         call ionc4_err(nc_err,ionc_rout, 'nf90_def_dim',ionc_nomlat)
         
         
         if (kmax .GT. 0) then
            nc_err = nf90_def_dim(nc_id,trim(ionc_nomz),dim_k,dim_k_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_nomz)
         end if
         
         if (dim_time .ne. -1) then
            nc_err = nf90_def_dim(nc_id,trim(ionc_nomtime),dim_time,dim_time_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',ionc_nomtime)
         endif
         
! ******** On definit les variables  **********
         
         tab_dim_lat(1) = dim_lat_id
         nc_err = nf90_def_var(nc_id, trim(ionc_nomlat), nf90_double, tab_dim_lat, lat_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomlat))
         
         tab_dim_lon(1) = dim_lon_id
         nc_err = nf90_def_var(nc_id, trim(ionc_nomlon), nf90_double, tab_dim_lon, lon_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomlon))
         
         if (dim_time .ne. -1) then
            tab_dim_time(1) = dim_time_id
            nc_err = nf90_def_var(nc_id, trim(ionc_nomtime), nf90_double, tab_dim_time, time_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomtime))
#ifdef MPI
            if (present(l_out_nc4par)) then
              if (l_out_nc4par) then
                nc_err = nf90_var_par_access(nc_id,time_id,nf90_collective)
                call ionc4_err(nc_err,ionc_rout,'nf90_var_par_access',trim(ionc_nomtime))
              end if
            end if
#endif
         endif
         
         if (kmax .GT. 0) then
            tab_dim_sigma(1) = dim_k_id
            nc_err = nf90_def_var(nc_id, trim(ionc_nomz), nf90_real, tab_dim_sigma, sig_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomz))
         endif
         
! ******** On definit les attributs de la loupe  **********
         nc_err = nf90_put_att(nc_id, nf90_global, 'limi', limi)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','limi')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'lima', lima)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','lima')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'pasi', pasi)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','pasi')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'ljmi', ljmi)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','ljmi')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'ljma', ljma)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','ljma')
         
         nc_err = nf90_put_att(nc_id, nf90_global, 'pasj', pasj)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','pasi')
         
         if (kmax .GT. 0) then
            nc_err = nf90_put_att(nc_id, nf90_global, 'lkmi', lkmi)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','lkmi')
            
            nc_err = nf90_put_att(nc_id, nf90_global, 'lkma', lkma)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','lkma')
            
            nc_err = nf90_put_att(nc_id, nf90_global, 'pask', pask)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','pask')
            
         endif
         
! ******** On definit les attributs des variables  **********
         nc_err = nf90_put_att(nc_id,lat_id,'long_name','latitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:long_name')
         
         nc_err = nf90_put_att(nc_id,lat_id,'standard_name','latitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:standard_name')
         
         nc_err = nf90_put_att(nc_id, lat_id,'units', 'degrees_north')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:units')
         
         nc_err = nf90_put_att(nc_id, lat_id,'axis','Y')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:axis')
         
         nc_err = nf90_put_att(nc_id,lat_id,'valid_min',REAL(-90.0,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',trim(ionc_nomlat))
         
         nc_err = nf90_put_att(nc_id,lat_id,'valid_max',REAL(90.0,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',trim(ionc_nomlat))
         
         nc_err = nf90_put_att(nc_id,lat_id,'_FillValue',REAL(1.7e+38,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue','LAT:_FillValue')

         nc_err = nf90_put_att(nc_id,lon_id,'long_name','longitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:long_name')
         
         nc_err = nf90_put_att(nc_id,lon_id,'standard_name','longitude')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:standard_name')
         
         nc_err = nf90_put_att(nc_id, lon_id,'units', 'degrees_east')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:units')
         
         nc_err = nf90_put_att(nc_id, lon_id,'axis','X')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:axis')
         
         nc_err = nf90_put_att(nc_id,lon_id,'valid_min',REAL(-180.0,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',trim(ionc_nomlon))
         
         nc_err = nf90_put_att(nc_id,lon_id,'valid_max',REAL(180.0,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att_double:valid_max',trim(ionc_nomlon))
         
         nc_err = nf90_put_att(nc_id,lon_id,'_FillValue',REAL(1.7e+38,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue','LON:_FillValue')

         if (kmax .GT. 0) then
            
            !nc_err = nf90_put_att(nc_id,sig_id,'units','level')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:units')
            
            nc_err = nf90_put_att(nc_id,sig_id, 'long_name', 'sigma level')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:long_name')
            
            nc_err = nf90_put_att(nc_id,sig_id,'standard_name','z_grid_index')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:standard_name')

            !nc_err = nf90_put_att(nc_id,sig_id,'content','z')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:content')

            !nc_err = nf90_put_att(nc_id,sig_id,'coordinates','level')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:coordinates')

            !nc_err = nf90_put_att(nc_id,sig_id,'associate','level')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:associate')

            nc_err = nf90_put_att(nc_id,sig_id,'axis','Z')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')
            
         endif
         
         if (dim_time .ne. -1) then
            nc_err = nf90_put_att(nc_id, time_id, 'long_name',ionc_longnamet)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:long_name')
            nc_err = nf90_put_att(nc_id, time_id, 'standard_name','time')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:standard_name')
            nc_err = nf90_put_att(nc_id, time_id, 'units',ionc_unitst)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:units')
            !nc_err = nf90_put_att(nc_id, time_id,'content','t')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:content')
            !nc_err = nf90_put_att(nc_id, time_id,'coordinates','time')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:coordinate')
            !nc_err = nf90_put_att(nc_id, time_id,'associate','time')
            !call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:associate')
            nc_err = nf90_put_att(nc_id, time_id,'axis','T')
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:axis')
            nc_err = nf90_put_att(nc_id, time_id, 'time_origin',ionc_originet)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:time_origin')
         endif
         
! ******** On peut maintenant quitter le mode definition *****************
         nc_err = nf90_enddef(nc_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_enddef',' ')
         
! ******** Ajout des lon, lat ********************************************
         
         indice = 1
         do ilon = limi,lima,pasi
            rlon(indice) = REAL(ilon-1,8)*dg+g0
            indice = indice + 1
         enddo
         nc_err = nf90_put_var(nc_id,lon_id,rlon)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomlon))
         
         indice = 1
         do ilat = ljmi,ljma,pasj
            rlat(indice) = REAL(ilat-1,8)*dfi+fi0
            indice = indice + 1
         enddo
         nc_err = nf90_put_var(nc_id,lat_id,rlat)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var_double',trim(ionc_nomlat))
         
! ******* Ajout de SIG(k)  ***************************************************
         
         if (kmax .GT. 0) then
            indice = 1
            do ik = lkmi,lkma,pask
               !rsig(indice) = sig(ik)
               rsig(indice) = ik
               indice = indice + 1
            enddo
            nc_err = nf90_put_var(nc_id,sig_id,rsig)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomz))
         endif
         
      endif
      
   else
      call ionc4_err(ionc_errcreation,ionc_rout, ' ' ,nom_fichier)
   endif
  
   deallocate(rlon,rlat)
   if(allocated(rsig)) deallocate(rsig)
 
   return
   end subroutine ionc4_createfile_1d_double
   
! ***************************************************************
! * subroutine ionc4_corres                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    16/02/01                                *
! * derniere modif :    20/02/01 par jfleroux                   *
! *      *
! * Role : rend le nc_id du fichier       *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier : chemin complet du fichier a *
! *      creer*
! *       sortie:- nc_id (0 -> pas de correspondance)*
! ***************************************************************
   subroutine ionc4_corres(nom_fichier,nc_id)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   integer :: nc_id
   
! ******** VARIABLES DE TRAVAIL **************
   integer  :: indfich
   
! ******** FIN DES DECLARATIONS **************
   
! ******test ionc_init
   
   if (ionc_init_ok .NE. 1) then
      write(*,*) 'ionc_init_ok= ',ionc_init_ok
      write(*,*) 'LIBIONC : ERREUR , initialisation avec ionc_init '
      write(*,*) '  obligatoire en debut de programme'
      stop
   endif
   
   
   nc_id = 0
   
   if (ionc_nfich .GT. 0) then
      do indfich = 1,ionc_nfich
         if (trim(ionc_nomfich(indfich)) .EQ. trim(nom_fichier)) then
            nc_id = ionc_idfich(indfich)
            EXIT
         endif
      enddo
   endif
   
   
   return
   end subroutine ionc4_corres
   
! ***************************************************************
! * subroutine ionc4_err                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    20/02/01 par jfleroux                   *
! *      *
! * Role : message d'erreur lie au acces netcdf      *
! *      *
! * Parametres :      *
! *       entree:- status    :  numero netcdf de l'erreur    *
! *       - croutine  :  nom de la routine        *
! *       - cfonction :  nom de la fonction netcdf        *
! *       - cvar      :  nom de la variable        *
! ***************************************************************
   subroutine ionc4_err(status,croutine,cfonction,cvar)
   
#ifdef MPI   
   include 'mpif.h'
#endif

! ******** PARAMETRES DE LA SUBROUTINE *******
   integer :: status,IERR_MPI
   character(len=*) :: croutine,cfonction,cvar
   
! ******** VARIABLES DE TRAVAIL **************
   
! ******** FIN DES DECLARATIONS **************
   
   if(status .ne. nf90_noerr) then
      write(*,*) 'Erreur : routine  = ' // trim(croutine)  // new_line('') // &
                '          fonction = ' // trim(cfonction) // new_line('') // &
                '          variable = ' //trim(cvar)
      if(status .gt. ionc_errcreation) then
         write(*,'(A,I5,A)') 'ERROR', status , ' = ' // trim(nf90_strerror(status))
      endif
#ifdef MPI
      call MPI_FINALIZE(IERR_MPI)
#endif
      stop
   endif
   end subroutine ionc4_err
   
   
! ***************************************************************
! * SUBROUTINE ionc4_createfile_sta_double                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/04/02                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)  *
! *      *
! * Role : creer le fichier Netcdf qui contiendra les donnees   *
! *        non maillees (stationnelles)                         *
! *         Initialiser les donnees ne dependant pas du time :*
! *      - sta*
! *      - SIG(k)*
! *         initialise les unites de temps et l'origine          *
! *        aux valeurs par defaut                               *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier : chemin complet du fichier a *
! *      creer*
! *      - lkmi : limite min de l'indice k de la loupe   *
! *      - lkma : limite max de l'indice k de la loupe   *
! *      - pask : pas de l'indice k de la loupe          *
! *      - kmax : borne max pour k                       *
! *      0 si 2D*
! *      - dim_sta : nombre de stations                  *
! *      - dim_time : nombre de pas de time que le *
! *           fichier contiendra si dim_time =0  *
! *           la dimension time est unlimited    *
! *           si dim_time = -1 pas de dim temps  *
! *      - SIG : niveaux*
! ***************************************************************
   subroutine ionc4_createfile_sta_double(nom_fichier,lkmi,lkma,pask,kmax, &
                                          dim_sta,dim_time,lon,lat)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier
   integer :: lkmi,lkma,pask,dim_time,dim_sta,kmax
   real(kind=8),dimension(dim_sta) :: lon,lat
   real(kind=4),dimension(:),allocatable :: rsig
   
! ******** VARIABLES DE TRAVAIL **************
   integer  :: dim_k,nc_id,sta_id,lat_id,lon_id,time_id,sig_id
   
! ** Declaration des IDS identifiant les dimensions des variables.
   integer  :: dim_sta_id,dim_time_id,dim_k_id
   
! ** Declarations des tableaux qui contiendront les IDS des dimensions des variables
   integer,dimension(1) :: tab_dim_sta,tab_dim_time,tab_dim_sigma
   integer  :: nc_err,ista,ik,indice
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_createfile_sta_double"
   
! *** Recuperation des dimensions
   
   if (kmax .GT. 0) then
      dim_k   = (lkma-lkmi)/pask +1
      allocate(rsig(dim_k))
   end if
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      
      ionc_nfich = ionc_nfich + 1
      if (ionc_nfich .gt. ionc_longtabfich) then
         call ionc4_err(ionc_errcreation,ionc_rout, ' ' ,nom_fichier)
         stop
      else
         ionc_nomfich(ionc_nfich) = nom_fichier
         
         nc_err = nf90_create(nom_fichier, nf90_share, nc_id)
         call ionc4_err(nc_err,ionc_rout, 'nf90_create',nom_fichier)
         
         ionc_idfich(ionc_nfich) = nc_id
         
         nc_err = nf90_def_dim(nc_id,'sta',dim_sta,dim_sta_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_dim','sta')
         
         if (kmax .GT. 0) then
            nc_err = nf90_def_dim(nc_id,trim(ionc_nomz),dim_k,dim_k_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',trim(ionc_nomz) )
         end if
         
         if (dim_time .ne. -1) then
            nc_err = nf90_def_dim(nc_id,trim(ionc_nomtime),dim_time,dim_time_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',trim(ionc_nomtime))
         endif
         
! ******** On definit les variables  **********
         
         tab_dim_sta(1) = dim_sta_id
         nc_err = nf90_def_var(nc_id, 'sta', nf90_int,tab_dim_sta, sta_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var','sta')
         
         if (dim_time .ne. -1) then
            tab_dim_time(1) = dim_time_id
            nc_err = nf90_def_var(nc_id, trim(ionc_nomtime),nf90_double,tab_dim_time, time_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomtime))
         endif
         
         if (kmax .GT. 0) then
            tab_dim_sigma(1) = dim_k_id
            nc_err = nf90_def_var(nc_id, trim(ionc_nomz),nf90_real, tab_dim_sigma, sig_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomz))
         end if
         
         nc_err = nf90_def_var(nc_id, trim(ionc_nomlat),nf90_double, tab_dim_sta, lat_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomlat))
         
         nc_err = nf90_def_var(nc_id, trim(ionc_nomlon),nf90_double, tab_dim_sta, lon_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomlon))
         
      end if
      
! ******** On definit les attributs des variables  **********
      nc_err = nf90_put_att(nc_id, lat_id,'long_name','latitude')
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:long_name')
      
      nc_err = nf90_put_att(nc_id, lat_id,'units', 'degrees_north')
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:units')
      
      nc_err = nf90_put_att(nc_id, lon_id,'long_name', 'longitude')
      call ionc4_err(nc_err,ionc_rout, 'nf90_put_att','LON:long_name')
      
      nc_err = nf90_put_att(nc_id, lon_id,'units', 'degrees_east')
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:units')
      
      nc_err = nf90_put_att(nc_id, sta_id,'long_name', 'STATION')
      call ionc4_err(nc_err,ionc_rout, 'nf90_put_att','STATION:long_name')
      
      if (kmax .GT. 0) then
         !nc_err = nf90_put_att(nc_id, sig_id,'units', 'level')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:units')
         
         nc_err = nf90_put_att(nc_id, sig_id, 'long_name ', 'sigma level')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:long_name')
         
         !nc_err = nf90_put_att(nc_id, sig_id, 'standard_name','z_grid_index')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:standard_name')
         
         nc_err = nf90_put_att(nc_id, sig_id,'axis', 'Z')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')

         nc_err = nf90_put_att(nc_id,sig_id,'c_grid_axis_shift',0.)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')
         
      end if
      
      if (dim_time .ne. -1) then
         nc_err = nf90_put_att(nc_id, time_id, 'long_name',ionc_longnamet)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:long_name')
         
         nc_err = nf90_put_att(nc_id, time_id, 'units', ionc_unitst)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:units')
         
         nc_err = nf90_put_att(nc_id, time_id, 'time_origin', ionc_originet)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:time_origin')
      endif
      
! ******** On peut maintenant quitter le mode definition *****************
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',' ')
      
! ******** Ajout des stations avec leurs coordonnees**********************
      
      indice = 1
      do ista = 1,dim_sta
         nc_err = nf90_put_var(nc_id,sta_id,indice,start=(/indice/))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var','STA')
         nc_err = nf90_put_var(nc_id,lat_id,lat(indice),start=(/indice/))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomlat))
         nc_err = nf90_put_var(nc_id,lon_id,lon(indice),start=(/indice/))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomlon))
         indice = indice + 1
      enddo
      
! ******* Ajout de SIG(k)  ***************************************************
      
      if (kmax .GT. 0) then
        indice = 1
        do ik = lkmi,lkma,pask
          !rsig(indice) = sig(ik)
          rsig(indice) = ik
          indice = indice + 1
        enddo
        nc_err = nf90_put_var(nc_id,sig_id,rsig)
        call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomz))
      endif
      
   else
      call ionc4_err(ionc_errcreation,ionc_rout, ' ' ,nom_fichier)
   endif
   
   return
   end subroutine ionc4_createfile_sta_double
   
! ***************************************************************
! * SUBROUTINE ionc4_createfile_sta_real                 *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/04/02                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)  *
! *      *
! * Role : creer le fichier Netcdf qui contiendra les donnees   *
! *        non maillees (stationnelles)                         *
! *         Initialiser les donnees ne dependant pas du time :*
! *      - sta*
! *      - SIG(k)*
! *         initialise les unites de temps et l'origine          *
! *        aux valeurs par defaut                               *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier : chemin complet du fichier a *
! *      creer*
! *      - lkmi : limite min de l'indice k de la loupe   *
! *      - lkma : limite max de l'indice k de la loupe   *
! *      - pask : pas de l'indice k de la loupe          *
! *      - kmax : borne max pour k                       *
! *      0 si 2D*
! *      - dim_sta : nombre de stations                  *
! *      - dim_time : nombre de pas de time que le *
! *           fichier contiendra si dim_time =0  *
! *           la dimension time est unlimited    *
! *           si dim_time = -1 pas de dim temps  *
! *      - SIG : niveaux*
! ***************************************************************
   subroutine ionc4_createfile_sta_real(nom_fichier,lkmi,lkma,pask,kmax, &
   dim_sta,dim_time,lon,lat)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier
   integer :: lkmi,lkma,pask,dim_time,dim_sta,kmax
   real(kind=4),dimension(dim_sta) :: lon,lat
   real(kind=4),dimension(:),allocatable :: rsig
   
! ******** VARIABLES DE TRAVAIL **************
   integer  :: dim_k,nc_id,sta_id,lat_id,lon_id,time_id,sig_id
   
! ** Declaration des IDS identifiant les dimensions des variables.
   integer  :: dim_sta_id,dim_time_id,dim_k_id
   
! ** Declarations des tableaux qui contiendront les IDS des dimensions des variables
   integer,dimension(1) :: tab_dim_sta,tab_dim_time,tab_dim_sigma
   integer  :: nc_err,ista,ik,indice
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_createfile_sta_real"
   
! *** Recuperation des dimensions
   
   if (kmax .GT. 0) then
      dim_k   = (lkma-lkmi)/pask +1
      allocate(rsig(dim_k))
   end if
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      
      ionc_nfich = ionc_nfich + 1
      if (ionc_nfich .gt. ionc_longtabfich) then
         call ionc4_err(ionc_errcreation,ionc_rout, ' ' ,nom_fichier)
         stop
      else
         ionc_nomfich(ionc_nfich) = nom_fichier
         
         nc_err = nf90_create(nom_fichier, nf90_share, nc_id)
         call ionc4_err(nc_err,ionc_rout, 'nf90_create',nom_fichier)
         
         ionc_idfich(ionc_nfich) = nc_id
         
         nc_err = nf90_def_dim(nc_id,'sta',dim_sta,dim_sta_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_dim','sta')
         
         if (kmax .GT. 0) then
            nc_err = nf90_def_dim(nc_id,trim(ionc_nomz),dim_k,dim_k_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',trim(ionc_nomz) )
         end if
         
         if (dim_time .ne. -1) then
            nc_err = nf90_def_dim(nc_id,trim(ionc_nomtime),dim_time,dim_time_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',trim(ionc_nomtime))
         endif
         
! ******** On definit les variables  **********
         
         tab_dim_sta(1) = dim_sta_id
         nc_err = nf90_def_var(nc_id, 'sta', nf90_int, tab_dim_sta, sta_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var','sta')
         
         if (dim_time .ne. -1) then
            tab_dim_time(1) = dim_time_id
            nc_err = nf90_def_var(nc_id, trim(ionc_nomtime),nf90_double,  tab_dim_time, time_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomtime))
         endif
         
         if (kmax .GT. 0) then
            tab_dim_sigma(1) = dim_k_id
            nc_err = nf90_def_var(nc_id, trim(ionc_nomz),nf90_real,  tab_dim_sigma, sig_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomz))
         end if
         
         nc_err = nf90_def_var(nc_id, trim(ionc_nomlat),nf90_real,  tab_dim_sta, lat_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomlat))
         
         nc_err = nf90_def_var(nc_id, trim(ionc_nomlon),nf90_real,  tab_dim_sta, lon_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomlon))
         
      end if
      
! ******** On definit les attributs des variables  **********
      nc_err = nf90_put_att(nc_id, lat_id,'long_name','latitude')
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:long_name')
      
      nc_err = nf90_put_att(nc_id, lat_id,'units', 'degrees_north')
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LAT:units')
      
      nc_err = nf90_put_att(nc_id, lon_id,'long_name', 'longitude')
      call ionc4_err(nc_err,ionc_rout, 'nf90_put_att','LON:long_name')
      
      nc_err = nf90_put_att(nc_id, lon_id,'units', 'degrees_east')
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att','LON:units')
      
      nc_err = nf90_put_att(nc_id, sta_id,'long_name', 'STATION')
      call ionc4_err(nc_err,ionc_rout, 'nf90_put_att','STATION:long_name')

      if (kmax .GT. 0) then
         !nc_err = nf90_put_att(nc_id, sig_id,'units', 'level')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:units')
         
         nc_err = nf90_put_att(nc_id, sig_id, 'long_name ', 'sigma level')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:long_name')
         
         !nc_err = nf90_put_att(nc_id, sig_id, 'standard_name ', 'z_grid_index')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:standard_name')
         
         nc_err = nf90_put_att(nc_id, sig_id,'axis', 'Z')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')

         nc_err = nf90_put_att(nc_id,sig_id,'c_grid_axis_shift',0.)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')

      end if
      
      if (dim_time .ne. -1) then
         nc_err = nf90_put_att(nc_id, time_id, 'long_name',ionc_longnamet)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:long_name')
         
         nc_err = nf90_put_att(nc_id, time_id, 'units', ionc_unitst)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:units')
         
         nc_err = nf90_put_att(nc_id, time_id, 'time_origin', ionc_originet)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:time_origin')
      endif
      
! ******** On peut maintenant quitter le mode definition *****************
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',' ')
      
! ******** Ajout des stations avec leurs coordonnees**********************
      
      indice = 1
      do ista = 1,dim_sta
         nc_err = nf90_put_var(nc_id,sta_id,indice,start=(/indice/))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var','STA')
         nc_err = nf90_put_var(nc_id,lat_id,lat(indice),start=(/indice/))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomlat))
         nc_err = nf90_put_var(nc_id,lon_id,lon(indice),start=(/indice/))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomlon))
         indice = indice + 1
      enddo
      
! ******* Ajout de SIG(k)  ***************************************************
      
      if (kmax .GT. 0) then
        indice = 1
        do ik = lkmi,lkma,pask
          !rsig(indice) = sig(ik)
          rsig(indice) = ik
          indice = indice + 1
        enddo
        nc_err = nf90_put_var(nc_id,sig_id,rsig)
        call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomz))
      endif
      
   else
      call ionc4_err(ionc_errcreation,ionc_rout, ' ' ,nom_fichier)
   endif
   
   return
   end subroutine ionc4_createfile_sta_real

! ****************************************************************
! * subroutine ionc4_createvar_real                              *
! *                                                              *
! * auteur         :    pgarreau,jfleroux                        *
! * org            :    IFREMER                                  *
! * date creation  :    xx/xx/xx                                 *
! * derniere modif :    21/12/10 par rramel                      *
! *                     genericite et nf90                       *
! * Role : creer une variable reelle simple precision            *
! *                                                              *
! * Parametres :                                                 *
! *       entree:- nom_fichier :  nom du fichier                 *
! *       - var_name    :  nom de la variable                    *
! *       - units       :  unite                                 *
! *       - long_name   :  nom long de la variable               *
! *       - dims        :  dimensions de la variable (optionnel) *
! *                        n'importe quelle combinaison de 'xyzt'*
! *                        (eg :'xt', 'xy', 'xyt', 'z', ... )    *
! *                        Si absent creation d'un scalaire      *
! *       - standard_name: nom standard (optionnel)              *
! *       - valid_min    : valeur minimale acceptable (optionnel)*
! *       - valid_max    : valeur maximale acceptable (optionnel)*
! *       - fill_value   : valeur manquante*
! *       - l_pack       : logique pour compression en short   *
! *       - l_out_nc4par     : logique pour ecriture // netcdf4  *
! *
! ****************************************************************
   subroutine ionc4_createvar_real(nom_fichier,var_name,units,&
                               long_name,standard_name,  &
                               location_in_cell,         &
                               l_auxcoordinate,          &
                               valid_min,valid_max,fill_value,&
                               dims,l_pack,l_out_nc4par)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),intent(in) :: nom_fichier,var_name,units,long_name
   character(len=*),intent(in),optional :: dims,standard_name,location_in_cell
   real(kind=4),intent(in),optional :: valid_min,valid_max
   real(kind=4),intent(in) :: fill_value
   logical,intent(in),optional :: l_out_nc4par,l_pack,l_auxcoordinate
   
! ******** VARIABLES DE TRAVAIL **************
   logical :: l_useless
   integer :: nc_id, dim_id,indice,dimlen
   integer :: var_id,nc_err,ndim_var,i,j,nf_type
   integer,DIMENSION(:),ALLOCATABLE :: tab_dim_var,chunksizes
   logical :: l_dim,packing
   character(len=30),dimension(4,2) :: tabdim
   character(len=30),dimension(4) :: ionc_nom
   character(len=40)              :: coordinate
   !character(len=40)              :: associate, content
   character(len=1)               :: location
   real(kind=4) :: scale_factor,add_offset

! ******** FONCTIONS **************
   integer trim

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_createvar_real"

#ifdef key_gfortran
   tabdim(:,1) = (/'x','y','z','t'/)
   tabdim(1,2) = ionc_nomlon
   tabdim(2,2) = ionc_nomlat
!cval because size of ionc_nomz different from ionc_nomlat and ionc_nomlon ? see comionc4
   tabdim(3,2) = ionc_nomz
   tabdim(4,2) = ionc_nomtime
   ionc_nom(1) = ionc_pnomlon
   ionc_nom(2) = ionc_pnomlat
   ionc_nom(3) = ionc_nomz
   ionc_nom(4) = ionc_nomtime
#else
   tabdim(:,1) = (/'x','y','z','t'/)
   tabdim(:,2) = (/ionc_nomlon,ionc_nomlat,ionc_nomz,ionc_nomtime/)
   ionc_nom = (/ionc_pnomlon,ionc_pnomlat,ionc_nomz,ionc_nomtime/)
#endif
   if (present(location_in_cell)) then
     location = TRIM(ADJUSTL(ADJUSTR(location_in_cell)))
     if (location=='u' .or. location=='v' .or. location=='f') then
       tabdim(1,2) = TRIM(tabdim(1,2))//'_'//location
       tabdim(2,2) = TRIM(tabdim(2,2))//'_'//location
#ifdef key_temp_virtualsed_out
     else if (location == 'w' .or. location == 's') then
#else
     else if (location == 'w') then
#endif
       tabdim(3,2) = TRIM(tabdim(3,2))//'_'//location
     end if
     if (location=='u') ionc_nom(1) = ionc_pnomlon_u
     if (location=='u') ionc_nom(2) = ionc_pnomlat_u
     if (location=='v') ionc_nom(1) = ionc_pnomlon_v
     if (location=='v') ionc_nom(2) = ionc_pnomlat_v
     if (location=='f') ionc_nom(1) = ionc_pnomlon_f
     if (location=='f') ionc_nom(2) = ionc_pnomlat_f
     if (location=='w') ionc_nom(3) = ionc_pnomz_w
#ifdef key_temp_virtualsed_out
     if (location=='s') ionc_nom(3) = ionc_pnomz_s
#endif
   end if
   

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else

      nf_type = nf90_real
      packing = .false.
! si on utilise la compression avec scale_fator et add_offset
      if (present(l_pack)) then
         if( l_pack             .and. &
             present(valid_min) .and. &
             present(valid_max)        ) packing = .true.
      endif

      if ( packing ) nf_type = nf90_short

      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      if(present(dims)) then
         ndim_var = len(dims)
         allocate(tab_dim_var(ndim_var))
         allocate(chunksizes(ndim_var))
         !associate=''
         coordinate=''
         !content=''
         indice = 0
         do j =1,size(tabdim,1)
            l_dim = .false.
            do i=1,ndim_var
               if (dims(i:i) == tabdim(j,1)) then
                  l_dim = .true.
                   exit
               endif
            enddo
            if (l_dim) then
               indice = indice + 1
               nc_err =  nf90_inq_dimid(nc_id,trim(ionc_nom(j)),dim_id)
               if (j <= 2) coordinate = TRIM(tabdim(j,2))//' '//TRIM(coordinate)
               !associate = TRIM(tabdim(j,2))//' '//TRIM(associate)
               !content = TRIM(tabdim(j,1))//TRIM(content)
               call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nom(j)))
               tab_dim_var(indice) = dim_id
            endif
         enddo
      endif

      if (present(dims)) then
#ifdef MPI
         if (present(l_out_nc4par)) then
            if (l_out_nc4par) then
               do i = 1,ndim_var
                  nc_err = nf90_inquire_dimension(nc_id, tab_dim_var(i),len=dimlen)
                  call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(var_name))
                  chunksizes(i)=max(dimlen,1)
               enddo
               nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id,chunksizes=chunksizes)
               nc_err = nf90_var_par_access(nc_id,var_id,nf90_collective)
            else
               nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
            endif
         else
            nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
         endif
#else
         if (present(l_out_nc4par)) l_useless = l_out_nc4par  ! useless for portability only
         nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
#endif
      else
         nc_err = nf90_def_var(nc_id, var_name, nf_type,var_id)
      endif

      call ionc4_err(nc_err,ionc_rout,'nf90_def_var',var_name)

! ******** On definit maintenant les attributs de cette variable *************

      nc_err = nf90_put_att(nc_id, var_id, 'long_name', long_name)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:long_name',var_name)
      if (present(standard_name)) then
         nc_err = nf90_put_att(nc_id, var_id, 'standard_name',standard_name)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:standard_name',var_name)
      endif
      nc_err = nf90_put_att(nc_id, var_id, 'units',units)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:units',var_name)
      
      if (present(dims)) then
        !nc_err = nf90_put_att(nc_id, var_id, 'associate', associate)
        !call ionc4_err(nc_err,ionc_rout,'nf90_put_att:associate',var_name)
        if(.not. present(l_auxcoordinate) ) then
          if (len(coordinate)>1) then
            nc_err = nf90_put_att(nc_id, var_id, 'coordinates', coordinate)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:coordinates',var_name)
          endif
        endif
        !nc_err = nf90_put_att(nc_id, var_id, 'content', content)
        !    call ionc4_err(nc_err,ionc_rout,'nf90_put_att:content',var_name)
      endif
! ****** Calcul du scale factor et de l'add_offset *********************
      if ( packing ) then
         scale_factor = (valid_max - valid_min) / (2**nb_bits - 2)
         add_offset   = (valid_max + valid_min) / 2
         nc_err = nf90_put_att(nc_id, var_id, 'scale_factor',scale_factor)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:scale_factor',var_name)
         nc_err = nf90_put_att(nc_id, var_id, 'add_offset',add_offset)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:add_offset',var_name)
      endif

      if (present(valid_min)) then
         if ( packing ) then
            nc_err = nf90_put_att(nc_id,var_id,'valid_min', nint((valid_min-add_offset)/scale_factor,2))
         else
            nc_err = nf90_put_att(nc_id,var_id,'valid_min', valid_min)
         endif
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
      endif

      if (present(valid_max)) then
         if ( packing ) then
            nc_err = nf90_put_att(nc_id,var_id,'valid_max', nint((valid_max-add_offset)/scale_factor,2))
         else
            nc_err = nf90_put_att(nc_id,var_id,'valid_max', valid_max)
         endif
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
      endif

      if ( packing ) then
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue', -HUGE(1_2)-1_2)
      else
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue', fill_value)
      endif
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)

      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif

   if (allocated(tab_dim_var)) deallocate(tab_dim_var)
   if (allocated(chunksizes)) deallocate(chunksizes)
   return
   end subroutine ionc4_createvar_real

! ****************************************************************
! * subroutine ionc4_createvar_double                            *
! *                                                              *
! * auteur         :    pgarreau,jfleroux                        *
! * org            :    IFREMER                                  *
! * date creation  :    xx/xx/xx                                 *
! * derniere modif :    21/12/10 par rramel                      *
! *                     genericite et nf90                       *
! * Role : creer une variable double                             *
! *                                                              *
! * Parametres :                                                 *
! *       entree:- nom_fichier :  nom du fichier                 *
! *       - var_name    :  nom de la variable                    *
! *       - units       :  unite                                 *
! *       - long_name   :  nom long de la variable               *
! *       - scale_factor:  coeff de multiplication               *
! *       - add_offset  :  offset                                *
! *       - dims        :  dimensions de la variable (optionnel) *
! *                        n'importe quelle combinaison de 'xyzt'*
! *                        (eg :'xt', 'xy', 'xyt', 'z', ... )    *
! *                        Si absent creation d'un scalaire      *
! *       - standard_name: nom standard (optionnel)              *
! *       - valid_min    : valeur minimale acceptable (optionnel)*
! *       - valid_max    : valeur maximale acceptable (optionnel)*
! *       - fill_value   : valeur manquante *
! *       - l_pack       : logique pour compression en short   *
! *       - l_out_nc4par     : logique pour ecriture // netcdf4  *
! *
! ****************************************************************
   subroutine ionc4_createvar_double(nom_fichier,var_name,units,&
                               long_name,standard_name,  &
                               location_in_cell,         &
                               l_auxcoordinate,          &
                               valid_min,valid_max,fill_value,&
                               dims,l_pack,l_out_nc4par)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),intent(in) :: nom_fichier,var_name,units,long_name
   character(len=*),intent(in),optional :: dims,standard_name,location_in_cell
   real(kind=8),intent(in),optional :: valid_min,valid_max
   real(kind=8),intent(in) :: fill_value
   logical,intent(in),optional :: l_out_nc4par,l_pack,l_auxcoordinate
   
! ******** VARIABLES DE TRAVAIL **************
   logical :: l_useless
   integer :: nc_id, dim_id,indice,dimlen
   integer :: var_id,nc_err,ndim_var,i,j,nf_type
   integer,DIMENSION(:),ALLOCATABLE :: tab_dim_var,chunksizes
   logical :: l_dim,packing
   character(len=30),dimension(4,2) :: tabdim
   character(len=20),dimension(4) :: ionc_nom
   !character(len=30)              :: associate,content
   character(len=30)              :: coordinate
   character(len=1)               :: location
   real(kind=8) :: scale_factor,add_offset

! ******** FONCTIONS **************
   integer trim

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_createvar_double"

#ifdef key_gfortran
   tabdim(:,1) = (/'x','y','z','t'/)
   tabdim(1,2) = ionc_nomlon
   tabdim(2,2) = ionc_nomlat
!cval because size of ionc_nomz different from ionc_nomlat and ionc_nomlon ? see comionc4
   tabdim(3,2) = ionc_nomz
   tabdim(4,2) = ionc_nomtime
   ionc_nom(1) = ionc_pnomlon
   ionc_nom(2) = ionc_pnomlat
   ionc_nom(3) = ionc_nomz
   ionc_nom(4) = ionc_nomtime
#else
   tabdim(:,1) = (/'x','y','z','t'/)
   tabdim(:,2) = (/ionc_nomlon,ionc_nomlat,ionc_nomz,ionc_nomtime/)
   ionc_nom = (/ionc_pnomlon,ionc_pnomlat,ionc_nomz,ionc_nomtime/)
#endif
   if (present(location_in_cell)) then
     location = TRIM(ADJUSTL(ADJUSTR(location_in_cell)))
     if (location=='u' .or. location=='v' .or. location=='f') then
       tabdim(1,2) = TRIM(tabdim(1,2))//'_'//location
       tabdim(2,2) = TRIM(tabdim(2,2))//'_'//location
#ifdef key_temp_virtualsed_out
     else if (location == 'w' .or. location == 's') then
#else
     else if (location == 'w') then
#endif
       tabdim(3,2) = TRIM(tabdim(3,2))//'_'//location
     end if
     if (location=='u') ionc_nom(1) = ionc_pnomlon_u
     if (location=='u') ionc_nom(2) = ionc_pnomlat_u
     if (location=='v') ionc_nom(1) = ionc_pnomlon_v
     if (location=='v') ionc_nom(2) = ionc_pnomlat_v
     if (location=='f') ionc_nom(1) = ionc_pnomlon_f
     if (location=='f') ionc_nom(2) = ionc_pnomlat_f
     if (location=='w') ionc_nom(3) = ionc_pnomz_w
#ifdef key_temp_virtualsed_out
     if (location=='s') ionc_nom(3) = ionc_pnomz_s
#endif
   end if
   

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nf_type = nf90_double
      packing = .false.
! si on utilise la compression avec scale_fator et add_offset
      if (present(l_pack)) then
         if( l_pack             .and. &
             present(valid_min) .and. &
             present(valid_max)        ) packing = .true.
      endif

      if ( packing ) nf_type = nf90_short

      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      if(present(dims)) then
         ndim_var = len(dims)
         allocate(tab_dim_var(ndim_var))
         allocate(chunksizes(ndim_var))
         !associate=''
         coordinate=''
         !content=''
         indice = 0
         do j =1,size(tabdim,1)
            l_dim = .false.
            do i=1,ndim_var
               if (dims(i:i) == tabdim(j,1)) then
                  l_dim = .true.
                   exit
               endif
            enddo
            if (l_dim) then
               indice = indice + 1
               if (j <= 2) coordinate = TRIM(tabdim(j,2))//' '//TRIM(coordinate)
               !associate = TRIM(tabdim(j,2))//' '//TRIM(associate)
               !content = TRIM(tabdim(j,1))//TRIM(content)      
               nc_err =  nf90_inq_dimid(nc_id,trim(ionc_nom(j)),dim_id)
               call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nom(j)))
               tab_dim_var(indice) = dim_id
            endif
         enddo
      endif

      if (present(dims)) then
#ifdef MPI
         if (present(l_out_nc4par)) then
            if (l_out_nc4par) then
               do i = 1,ndim_var
                  nc_err = nf90_inquire_dimension(nc_id, tab_dim_var(i),len=dimlen)
                  call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(var_name))
                  chunksizes(i)=max(dimlen,1)
               enddo
               nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id,chunksizes=chunksizes)
               nc_err = nf90_var_par_access(nc_id,var_id,nf90_collective)
            else
               nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
            endif
         else
            nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
         endif
#else
         if (present(l_out_nc4par)) l_useless = l_out_nc4par  ! useless for portability only
         nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
#endif
      else
         nc_err = nf90_def_var(nc_id, var_name, nf_type,var_id)
      endif

      call ionc4_err(nc_err,ionc_rout,'nf90_def_var',var_name)

! ******** On definit maintenant les attributs de cette variable *************

      nc_err = nf90_put_att(nc_id, var_id, 'long_name', long_name)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:long_name',var_name)
      if (present(standard_name)) then
         nc_err = nf90_put_att(nc_id, var_id, 'standard_name',standard_name)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:standard_name',var_name)
      endif
      nc_err = nf90_put_att(nc_id, var_id, 'units',units)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:units',var_name)

      if (present(dims)) then
        !nc_err = nf90_put_att(nc_id, var_id, 'associate', associate)
        !call ionc4_err(nc_err,ionc_rout,'nf90_put_att:associate',var_name)
        if(.not. present(l_auxcoordinate) ) then
          if (len(coordinate)>1) then
            nc_err = nf90_put_att(nc_id, var_id, 'coordinates', coordinate)
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:coordinates',var_name)
          endif
        endif
         !nc_err = nf90_put_att(nc_id, var_id, 'content', content)
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att:content',var_name)
      endif
! ****** Calcul du scale factor et de l'add_offset *********************
      if ( packing ) then
         scale_factor = (valid_max - valid_min) / (2**nb_bits - 2)
         add_offset   = (valid_max + valid_min) / 2
         nc_err = nf90_put_att(nc_id, var_id, 'scale_factor',scale_factor)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:scale_factor',var_name)
         nc_err = nf90_put_att(nc_id, var_id, 'add_offset',add_offset)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:add_offset',var_name)
      endif

      if (present(valid_min)) then
         if ( packing ) then
            nc_err = nf90_put_att(nc_id,var_id,'valid_min', nint((valid_min-add_offset)/scale_factor,2))
         else
            nc_err = nf90_put_att(nc_id,var_id,'valid_min', valid_min)
         endif
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
      endif

      if (present(valid_max)) then
         if ( packing ) then
            nc_err = nf90_put_att(nc_id,var_id,'valid_max', nint((valid_max-add_offset)/scale_factor,2))
         else
            nc_err = nf90_put_att(nc_id,var_id,'valid_max', valid_max)
         endif
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
      endif

      if ( packing ) then
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue', -HUGE(1_2)-1_2)
      else
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue', fill_value)
      endif
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)

      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif

   if (allocated(tab_dim_var)) deallocate(tab_dim_var)
   if (allocated(chunksizes)) deallocate(chunksizes)
   return
   end subroutine ionc4_createvar_double

! ****************************************************************
! * subroutine ionc4_createvar_int                            *
! *                                                              *
! * auteur         :    pgarreau,jfleroux                        *
! * org            :    IFREMER                                  *
! * date creation  :    xx/xx/xx                                 *
! * derniere modif :    21/12/10 par rramel                      *
! *                     genericite et nf90                       *
! * Role : creer une variable entiere !!                      *
! *                                                              *
! * Parametres :                                                 *
! *       entree:- nom_fichier :  nom du fichier                 *
! *       - var_name    :  nom de la variable                    *
! *       - units       :  unite                                 *
! *       - long_name   :  nom long de la variable               *
! *       - scale_factor:  coeff de multiplication               *
! *       - add_offset  :  offset                                *
! *       - dims        :  dimensions de la variable (optionnel) *
! *                        n'importe quelle combinaison de 'xyzt'*
! *                        (eg :'xt', 'xy', 'xyt', 'z', ... )    *
! *                        Si absent creation d'un scalaire      *
! *       - standard_name: nom standard (optionnel)              *
! *       - valid_min    : valeur minimale acceptable (optionnel)*
! *       - valid_max    : valeur maximale acceptable (optionnel)*
! *       - fill_value   : valeur manquante *
! *       - l_pack       : logique pour compression en short   *
! *       - l_out_nc4par     : logique pour ecriture // netcdf4  *
! *
! ****************************************************************
   subroutine ionc4_createvar_int(nom_fichier,var_name,units,&
                               long_name,standard_name,  &
                               location_in_cell,         &
                               l_auxcoordinate,          &
                               valid_min,valid_max,fill_value,&
                               dims,l_out_nc4par)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),intent(in) :: nom_fichier,var_name,units,long_name
   character(len=*),intent(in),optional :: dims,standard_name,location_in_cell
   integer,intent(in),optional :: valid_min,valid_max
   integer,intent(in) :: fill_value
   logical,intent(in),optional :: l_out_nc4par,l_auxcoordinate
   
! ******** VARIABLES DE TRAVAIL **************
   logical :: l_useless
   integer :: nc_id, dim_id,indice,dimlen
   integer :: var_id,nc_err,ndim_var,i,j,nf_type
   integer,DIMENSION(:),ALLOCATABLE :: tab_dim_var,chunksizes
   logical :: l_dim
   character(len=30),dimension(4,2) :: tabdim
   character(len=20),dimension(4) :: ionc_nom
   !character(len=30)              :: associate,content
   character(len=30)              :: coordinate
   character(len=1)               :: location

! ******** FONCTIONS **************
   integer trim

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_createvar_int"

#ifdef key_gfortran
   tabdim(:,1) = (/'x','y','z','t'/)
   tabdim(1,2) = ionc_nomlon
   tabdim(2,2) = ionc_nomlat
!cval because size of ionc_nomz different from ionc_nomlat and ionc_nomlon ? see comionc4
   tabdim(3,2) = ionc_nomz
   tabdim(4,2) = ionc_nomtime
   ionc_nom(1) = ionc_pnomlon
   ionc_nom(2) = ionc_pnomlat
   ionc_nom(3) = ionc_nomz
   ionc_nom(4) = ionc_nomtime
#else
   tabdim(:,1) = (/'x','y','z','t'/)
   tabdim(:,2) = (/ionc_nomlon,ionc_nomlat,ionc_nomz,ionc_nomtime/)
   ionc_nom = (/ionc_pnomlon,ionc_pnomlat,ionc_nomz,ionc_nomtime/)
#endif
   if (present(location_in_cell)) then
     location = TRIM(ADJUSTL(ADJUSTR(location_in_cell)))
     if (location=='u' .or. location=='v' .or. location=='f') then
       tabdim(1,2) = TRIM(tabdim(1,2))//'_'//location
       tabdim(2,2) = TRIM(tabdim(2,2))//'_'//location
#ifdef key_temp_virtualsed_out
     else if (location == 'w' .or. location == 's') then
#else
     else if (location == 'w') then
#endif
       tabdim(3,2) = TRIM(tabdim(3,2))//'_'//location
     end if
     if (location=='u') ionc_nom(1) = ionc_pnomlon_u
     if (location=='u') ionc_nom(2) = ionc_pnomlat_u
     if (location=='v') ionc_nom(1) = ionc_pnomlon_v
     if (location=='v') ionc_nom(2) = ionc_pnomlat_v
     if (location=='f') ionc_nom(1) = ionc_pnomlon_f
     if (location=='f') ionc_nom(2) = ionc_pnomlat_f
     if (location=='w') ionc_nom(3) = ionc_pnomz_w
#ifdef key_temp_virtualsed_out
     if (location=='s') ionc_nom(3) = ionc_pnomz_s
#endif
   end if
   

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nf_type = nf90_int
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      if(present(dims)) then
         ndim_var = len(dims)
         allocate(tab_dim_var(ndim_var))
         allocate(chunksizes(ndim_var))
         !associate=''
         coordinate=''
         !content=''

         indice = 0
         do j =1,size(tabdim,1)
            l_dim = .false.
            do i=1,ndim_var
               if (dims(i:i) == tabdim(j,1)) then
                  l_dim = .true.
                   exit
               endif
            enddo
            if (l_dim) then
               indice = indice + 1
               if (j <= 2) coordinate = TRIM(tabdim(j,2))//' '//TRIM(coordinate)
               !associate = TRIM(tabdim(j,2))//' '//TRIM(associate)
               !content = TRIM(tabdim(j,1))//TRIM(content)      
               nc_err =  nf90_inq_dimid(nc_id,trim(ionc_nom(j)),dim_id)
               call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nom(j)))
               tab_dim_var(indice) = dim_id
            endif
         enddo
      endif

      if (present(dims)) then
#ifdef MPI
         if (present(l_out_nc4par)) then
            if (l_out_nc4par) then
               do i = 1,ndim_var
                  nc_err = nf90_inquire_dimension(nc_id, tab_dim_var(i),len=dimlen)
                  call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(var_name))
                  chunksizes(i)=max(dimlen,1)
               enddo
               nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id,chunksizes=chunksizes)
               nc_err = nf90_var_par_access(nc_id,var_id,nf90_collective)
            else
               nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
            endif
         else
            nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
         endif
#else
         if (present(l_out_nc4par)) l_useless = l_out_nc4par  ! useless for portability only
         nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
#endif
      else
         nc_err = nf90_def_var(nc_id, var_name, nf_type,var_id)
      endif

      call ionc4_err(nc_err,ionc_rout,'nf90_def_var',var_name)

! ******** On definit maintenant les attributs de cette variable *************

      nc_err = nf90_put_att(nc_id, var_id, 'long_name', long_name)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:long_name',var_name)
      if (present(standard_name)) then
         nc_err = nf90_put_att(nc_id, var_id, 'standard_name',standard_name)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:standard_name',var_name)
      endif
      nc_err = nf90_put_att(nc_id, var_id, 'units',units)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:units',var_name)

      if (present(dims)) then
         !nc_err = nf90_put_att(nc_id, var_id, 'associate', associate)
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att:associate',var_name)
         if(.not. present(l_auxcoordinate) ) then
           if (len(coordinate)>1) then
             nc_err = nf90_put_att(nc_id, var_id, 'coordinates', coordinate)
             call ionc4_err(nc_err,ionc_rout,'nf90_put_att:coordinates',var_name)
           endif
         endif
         !nc_err = nf90_put_att(nc_id, var_id, 'content', content)
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att:content',var_name)
      endif
      
      if (present(valid_min)) then  
         nc_err = nf90_put_att(nc_id,var_id,'valid_min', valid_min)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
      endif

      if (present(valid_max)) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_max', valid_max)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
      endif

      nc_err = nf90_put_att(nc_id,var_id,'_FillValue', fill_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)

      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif

   if (allocated(tab_dim_var)) deallocate(tab_dim_var)
   if (allocated(chunksizes)) deallocate(chunksizes)
   return
   end subroutine ionc4_createvar_int
   
! ***************************************************************
! * subroutine ionc4_write_zxyt_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 4D = lon, lat, z, temps    *
! *        a partir d'une variable reelle 3D = z , lon, lat     *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - kmax        :  borne max dim k*
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! *       - l_out_nc4par     : logique I/O parallel
! ***************************************************************
   subroutine ionc4_write_zxyt_real(nom_fichier,var_name,var,imin,imax,jmin,jmax,kmin,kmax,nrec,fill_value,l_out_nc4par,sed)

! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,kmin,kmax,nrec
   real(kind=4),dimension(kmin:kmax,imin:imax,jmin:jmax) :: var
   real(kind=4),optional,intent(in) :: fill_value
   logical,optional,intent(in) :: l_out_nc4par,sed
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id, nc_err,dim_time_id,dim_time
   real(kind=4) :: add_offset,scale_factor
   integer :: dimi,dimj,dimk,deck
   integer :: i,j,k,indi,indj,indk
   integer,DIMENSION(4) :: start,count
   real(kind=4),dimension(:,:,:),allocatable :: bid
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_zxyt_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      call ionc4_readbounds(nc_id)
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      if (nrec .eq. 0) then
         nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
         call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         dim_time = nrec
      endif

      start = (/1,1,1,dim_time/)
! * on traite le cas des I/O paralleles en localisant le sous-domaine MPI
! * dans le fichier de sortie global     
      if (present(l_out_nc4par)) then
         if (l_out_nc4par) then
            call ionc4_locbounds(nc_id,imin,imax,jmin,jmax)
            start = (/ int((ionc_limi-ionc_global_imin)/ionc_pasi) + 1 , &
                       int((ionc_ljmi-ionc_global_jmin)/ionc_pasj) + 1 , &
                       1,dim_time /)
         endif
      endif

      dimi = (ionc_lima - ionc_limi)/ ionc_pasi + 1
      dimj = (ionc_ljma - ionc_ljmi)/ ionc_pasj + 1
      dimk= (ionc_lkma - ionc_lkmi)/ ionc_pask + 1


      if (present(sed)) then
        dimk=kmax
        ionc_lkmi=1
        ionc_lkma=kmax
        ionc_pask=1
      end if       

      deck=0 
      if (kmin==0 .AND. ionc_lkmi==1 .AND. ionc_pask==1) then
          dimk = dimk + 1
          deck = -1
      end if

      allocate(bid(dimi,dimj,dimk))
      count = (/dimi,dimj,dimk,1/)

      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
 
      if (scale_factor /= 1.e0 .and. present(fill_value)) then
         indk = 1
         do k=ionc_lkmi+deck,ionc_lkma, ionc_pask
            indj = 1
            do j=ionc_ljmi,ionc_ljma, ionc_pasj
               indi = 1
               do i=ionc_limi,ionc_lima, ionc_pasi
                  if (var(i,j,k) /= fill_value) then
                     bid(indi,indj,indk)=real(nint((var(i,j,k)- add_offset)/scale_factor))
                  else
                     ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
                     bid(indi,indj,indk)=-HUGE(1_2)-1.0_4
                  endif
                  indi = indi + 1
               enddo
               indj = indj + 1
            enddo
            indk = indk + 1
         enddo
      else
         indk = 1
         do k=ionc_lkmi+deck,ionc_lkma, ionc_pask
            indj = 1
            do j=ionc_ljmi,ionc_ljma, ionc_pasj
               indi = 1
               do i=ionc_limi,ionc_lima, ionc_pasi
                  bid(indi,indj,indk)=var(i,j,k)
                  indi = indi + 1
               enddo
               indj = indj + 1
            enddo
            indk = indk + 1
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      
      deallocate(bid)
   endif
   
   return
   end subroutine ionc4_write_zxyt_real
   
!***************************************************************
! * subroutine ionc4_write_zxyt_double                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 4D = lon, lat, z, temps    *
! *        a partir d'une variable reelle 3D = z , lon, lat     *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - kmax        :  borne max dim k*
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! ***************************************************************
   subroutine ionc4_write_zxyt_double(nom_fichier,var_name,var,imin,imax,jmin,jmax,kmin,kmax,nrec,fill_value,l_out_nc4par,sed)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,kmin,kmax,nrec
   real(kind=8) ,dimension(kmin:kmax,imin:imax,jmin:jmax) :: var
   real(kind=8),intent(in),optional :: fill_value
   logical,optional,intent(in) :: l_out_nc4par,sed
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id, nc_err,dim_time_id,dim_time
   real(kind=8) :: add_offset,scale_factor
   integer :: dimi,dimj,dimk
   integer :: i,j,k,indi,indj,indk,deck
   integer,DIMENSION(4) :: start,count
   real(kind=8),dimension(:,:,:),allocatable :: bid
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_zxyt_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      call ionc4_readbounds(nc_id)
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      if (nrec .eq. 0) then
         nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
         call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         dim_time = nrec
      endif

      start = (/1,1,1,dim_time/)
! * on traite le cas des I/O paralleles en localisant le sous-domaine MPI
! * dans le fichier de sortie global
      if (present(l_out_nc4par)) then
         if (l_out_nc4par) then
            call ionc4_locbounds(nc_id,imin,imax,jmin,jmax)
            start = (/ int((ionc_limi-ionc_global_imin)/ionc_pasi) + 1 , &
                       int((ionc_ljmi-ionc_global_jmin)/ionc_pasj) + 1 , &
                       1,dim_time /)
         endif
      endif
      
      dimi = (ionc_lima - ionc_limi)/ ionc_pasi + 1
      dimj = (ionc_ljma - ionc_ljmi)/ ionc_pasj + 1
      dimk= (ionc_lkma - ionc_lkmi)/ ionc_pask + 1

      if (present(sed)) then
        dimk=kmax
        ionc_lkmi=1
        ionc_lkma=kmax
        ionc_pask=1
      end if       

      deck=0 
      if (kmin==0 .AND. ionc_lkmi==1 .AND. ionc_pask==1) then
          dimk = dimk + 1
          deck = -1
      end if

      allocate(bid(dimi,dimj,dimk))
      count = (/dimi,dimj,dimk,1/)
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.d0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.d0
      endif
      
      if (scale_factor /= 1.d0 .and. present(fill_value)) then   
         indk = 1
         do k=ionc_lkmi+deck,ionc_lkma, ionc_pask
            indi = 1
            do i=ionc_limi,ionc_lima, ionc_pasi
               indj = 1
               do j=ionc_ljmi,ionc_ljma, ionc_pasj
                  if (var(i,j,k) /= fill_value) then
                     bid(indi,indj,indk)=REAL(nint((var(i,j,k)- add_offset)/scale_factor),8)
                  else
                     ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
                     bid(indi,indj,indk)=-HUGE(1_2)-1.0_8
                  endif
                  indj = indj + 1
               enddo
               indi = indi + 1
            enddo
            indk = indk + 1
         enddo
      else
         indk = 1
         do k=ionc_lkmi+deck,ionc_lkma, ionc_pask
            indi = 1
            do i=ionc_limi,ionc_lima, ionc_pasi
               indj = 1
               do j=ionc_ljmi,ionc_ljma, ionc_pasj
                  bid(indi,indj,indk)=var(i,j,k)
                  indj = indj + 1
               enddo
               indi = indi + 1
            enddo
            indk = indk + 1
         enddo
      endif      
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid)
      
   endif
   
   return
   end subroutine ionc4_write_zxyt_double
   
! ***************************************************************
! * subroutine ionc4_readbounds                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    22/02/01                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : lit les parametres de loupe dans le fichier netcdf   *
! *      *
! * Parametres :      *
! *       entree:- nc_id :  id du fichier        *
! ***************************************************************
   subroutine ionc4_readbounds(nc_id)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   integer :: nc_id, dim_lon_id,dim_lat_id,dim_z_id
   
! ******** VARIABLES DE TRAVAIL **************
  
   CHARACTER(LEN=30) :: ionc_rout_loc
   integer :: nc_err, nc_err1, nc_err2
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout_loc = "ionc4_readbounds"
   
   nc_err = nf90_get_att(nc_id, nf90_global,'limi',ionc_limi)
   if(nc_err .ne. nf90_noerr) ionc_limi = 1
   
   nc_err = nf90_get_att(nc_id, nf90_global,'lima',ionc_lima)
   if(nc_err .ne. nf90_noerr) then
      nc_err1 = nf90_inq_dimid(nc_id,trim(ionc_pnomlon),dim_lon_id)
      call ionc4_err(nc_err1,ionc_rout_loc, 'nf90_inq_dimid',trim(ionc_pnomlon))
      nc_err1 = nf90_inquire_dimension(nc_id,dim_lon_id,len=ionc_lima)
      call ionc4_err(nc_err1,ionc_rout_loc,'nf90_inquire_dimension',trim(ionc_pnomlon))
   end if
   
   nc_err = nf90_get_att(nc_id, nf90_global,'pasi',ionc_pasi)
   if(nc_err .ne. nf90_noerr) ionc_pasi = 1
   
   nc_err = nf90_get_att(nc_id, nf90_global,'ljmi',ionc_ljmi)
   if(nc_err .ne. nf90_noerr) ionc_ljmi = 1
   
   nc_err = nf90_get_att(nc_id, nf90_global,'ljma',ionc_ljma)
   if(nc_err .ne. nf90_noerr) then
      nc_err1 = nf90_inq_dimid(nc_id,trim(ionc_pnomlat),dim_lat_id)
      call ionc4_err(nc_err1,ionc_rout_loc, 'nf90_inq_dimid',trim(ionc_pnomlat))
      nc_err1 = nf90_inquire_dimension(nc_id,dim_lat_id,len=ionc_ljma)
      call ionc4_err(nc_err1,ionc_rout_loc,'nf90_inquire_dimension',trim(ionc_pnomlat))
   end if
   
   nc_err = nf90_get_att(nc_id, nf90_global,'pasj',ionc_pasj)
   if(nc_err .ne. nf90_noerr) ionc_pasj = 1
   
   nc_err = nf90_get_att(nc_id, nf90_global,'lkmi',ionc_lkmi)
   if(nc_err .ne. nf90_noerr) ionc_lkmi = 1
   
   nc_err = nf90_get_att(nc_id, nf90_global,'lkma',ionc_lkma)
   if(nc_err .ne. nf90_noerr) then
      nc_err1 = nf90_inq_dimid(nc_id,trim(ionc_nomz),dim_z_id)
      if(nc_err1 .ne. nf90_noerr) then
         ionc_lkma = 0
      else
         nc_err2 = nf90_inquire_dimension(nc_id,dim_z_id,len=ionc_lkma)
         call ionc4_err(nc_err,ionc_rout_loc,'nf90_inquire_dimension',trim(ionc_nomz))
      end if
   end if
   
   nc_err = nf90_get_att(nc_id, nf90_global,'pask',ionc_pask)
   if(nc_err .ne. nf90_noerr) then
      ionc_pask = 1
   end if
   
   return
   end subroutine ionc4_readbounds

! ***************************************************************
! * subroutine ionc4_locbounds                   *
! *      *
! * auteur         :    rramel                                *
! * org            :    ALYOTECH                                 *
! * date creation  :    18/01/11                                *
! * derniere modif :    18/01/11 par rramel (Alyotech)                   *
! *      *
! * Role : localise un sous-domaine MPI dans un fichier global  *
! *      * modification de ionc_limi & ionc_lima
! * Parametres :      *
! *       entree:
! *       - nom_fichier :  nom du fichier        *
! *       - temps       :  valeur du temps        *
! *       - nrec        :  numero d'enregistrement*
! ***************************************************************
   subroutine ionc4_locbounds(nc_id,imin,imax,jmin,jmax)
  
   
! ******** PARAMETRES DE LA SUBROUTINE *******

   integer :: nc_id,imin,imax,jmin,jmax

! ******** VARIABLES DE TRAVAIL **************

   CHARACTER(LEN=30) :: ionc_rout_loc
   integer :: nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout_loc = "ionc4_locbounds"

   nc_err = nf90_get_att(nc_id,nf90_global,'global_imin',ionc_global_imin)
   call ionc4_err(nc_err,ionc_rout_loc,'nf90_get_att','global_imin')
   nc_err = nf90_get_att(nc_id,nf90_global,'global_imax',ionc_global_imax)
   call ionc4_err(nc_err,ionc_rout_loc,'nf90_get_att','global_imax')
   nc_err = nf90_get_att(nc_id,nf90_global,'global_jmin',ionc_global_jmin)
   call ionc4_err(nc_err,ionc_rout_loc,'nf90_get_att','global_jmin')
   nc_err = nf90_get_att(nc_id,nf90_global,'global_jmax',ionc_global_jmax)
   call ionc4_err(nc_err,ionc_rout_loc,'nf90_get_att','global_jmax')
   nc_err = nf90_get_att(nc_id,nf90_global,'pasi',ionc_pasi)
   call ionc4_err(nc_err,ionc_rout_loc,'nf90_get_att','pasi')
   nc_err = nf90_get_att(nc_id,nf90_global,'pasj',ionc_pasj)
   call ionc4_err(nc_err,ionc_rout_loc,'nf90_get_att','pasj')

   ionc_limi = max(imin,ionc_global_imin)
   do while (  real(ionc_limi-ionc_global_imin)/ionc_pasi &
            /= real((ionc_limi-ionc_global_imin)/ionc_pasi) )
      ionc_limi = ionc_limi + 1
   enddo
   ionc_lima = min(imax,ionc_global_imax)
   do while (  real(ionc_lima-ionc_global_imin)/ionc_pasi &
            /= real((ionc_lima-ionc_global_imin)/ionc_pasi) )
      ionc_lima = ionc_lima - 1
   enddo

   ionc_ljmi = max(jmin,ionc_global_jmin)
   do while (  real(ionc_ljmi-ionc_global_jmin)/ionc_pasj &
            /= real((ionc_ljmi-ionc_global_jmin)/ionc_pasj) )
      ionc_ljmi = ionc_ljmi + 1
   enddo
   ionc_ljma = min(jmax,ionc_global_jmax)
   do while (  real(ionc_ljma-ionc_global_jmin)/ionc_pasj &
            /= real((ionc_ljma-ionc_global_jmin)/ionc_pasj) )
      ionc_ljma = ionc_ljma - 1
   enddo

   return

   end subroutine ionc4_locbounds

! ***************************************************************
! * subroutine ionc4_write_time                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    26/02/01                                *
! * derniere modif :    26/02/01 par jfleroux                   *
! *      *
! * Role : ecrit le temps                      *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - temps       :  valeur du temps        *
! *       - nrec        :  numero d'enregistrement*
! ***************************************************************
   subroutine ionc4_write_time(nom_fichier,nrec,temps)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier
   real(kind=8) :: temps
   integer :: nrec
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id,dim_time_id,dim_time,nc_err
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_time"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      if (nrec .eq. 0) then
         nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
         call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
         dim_time = dim_time + 1
      else
         dim_time = nrec
      endif
      
      nc_err = nf90_inq_varid(nc_id,trim(ionc_nomtime),var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',trim(ionc_nomtime))
      nc_err = nf90_put_var(nc_id,var_id,temps,start=(/dim_time/))
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomtime))

   endif
   
   return
   end subroutine ionc4_write_time
   
! ***************************************************************
! * subroutine ionc4_write_xyt_double                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 3D = lon, lat, temps    *
! *        a partir d'une variable reelle 2D =  lon, lat     *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! *       - l_out_nc4par     : logique pour I/O parallel
! ***************************************************************
   subroutine ionc4_write_xyt_double(nom_fichier,var_name,var,imin,imax,jmin,jmax,nrec,fill_value,l_out_nc4par)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,nrec
   real(kind=8) ,dimension(imin:imax,jmin:jmax) :: var
   real(kind=8),optional,intent(in) :: fill_value
   logical,optional,intent(in) :: l_out_nc4par 
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id, nc_err,dim_time_id,dim_time
   integer :: i,j,indi,indj,dimi,dimj
   real(kind=8) :: add_offset,scale_factor
   integer,DIMENSION(3) :: start,count
   real(kind=8),dimension(:,:),allocatable :: bid
!
!******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_xyt_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      call ionc4_readbounds(nc_id)
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      if (nrec .eq. 0) then
         nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
         call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         dim_time = nrec
      endif
      
      start  = (/1,1,dim_time/)
! * on traite le cas des I/O paralleles en localisant le sous-domaine MPI
! * dans le fichier de sortie global
      if (present(l_out_nc4par)) then
         if (l_out_nc4par) then
            call ionc4_locbounds(nc_id,imin,imax,jmin,jmax)
            start = (/ int((ionc_limi-ionc_global_imin)/ionc_pasi) + 1 , &
                       int((ionc_ljmi-ionc_global_jmin)/ionc_pasj) + 1 , &
                       dim_time /)
         endif
      endif

      dimi = (ionc_lima - ionc_limi)/ ionc_pasi + 1
      dimj = (ionc_ljma - ionc_ljmi)/ ionc_pasj + 1

      allocate(bid(dimi,dimj))
      count = (/dimi,dimj,1/)
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.d0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.d0
      endif
      
      if (scale_factor /= 1.d0 .and. present(fill_value)) then     
         indi = 1
         do i=ionc_limi,ionc_lima, ionc_pasi
            indj = 1
            do j=ionc_ljmi,ionc_ljma, ionc_pasj
               if (var(i,j) /= fill_value) then
                  bid(indi,indj)=REAL(nint((var(i,j)- add_offset)/scale_factor),8)
               else
                  ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
                  bid(indi,indj)=-HUGE(1_2)-1.0_8
               endif
               indj = indj + 1
            enddo
            indi = indi + 1
         enddo
      else
         indi = 1
         do i=ionc_limi,ionc_lima, ionc_pasi
            indj = 1
            do j=ionc_ljmi,ionc_ljma, ionc_pasj
               bid(indi,indj)=var(i,j)
               indj = indj + 1
            enddo
            indi = indi + 1
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid)
   endif
   
   return
   end subroutine ionc4_write_xyt_double
   
! ***************************************************************
! * subroutine ionc4_write_xyt_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 3D = lon, lat, temps    *
! *        a partir d'une variable reelle 2D = lon, lat     *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! *       - l_out_nc4par     :  logique pour I/O parallel
! ***************************************************************
   subroutine ionc4_write_xyt_real(nom_fichier,var_name,var,imin,imax,jmin,jmax,nrec,fill_value,l_out_nc4par)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,nrec
   real(kind=4),dimension(imin:imax,jmin:jmax) :: var
   real(kind=4),optional,intent(in) :: fill_value
   logical,optional,intent(in) :: l_out_nc4par 
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id, nc_err,dim_time_id,dim_time
   integer :: i,j,indi,indj,dimi,dimj
   real(kind=4) :: add_offset,scale_factor
   integer,DIMENSION(3) :: start,count
   real(kind=4),dimension(:,:),allocatable :: bid
!
!******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_xyt_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      call ionc4_readbounds(nc_id)
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      if (nrec .eq. 0) then
         nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
         call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         dim_time = nrec
      endif

      start  = (/1,1,dim_time/)
! * on traite le cas des I/O paralleles en localisant le sous-domaine MPI
! * dans le fichier de sortie global
      if (present(l_out_nc4par)) then
         if (l_out_nc4par) then
            call ionc4_locbounds(nc_id,imin,imax,jmin,jmax)
            start = (/ int((ionc_limi-ionc_global_imin)/ionc_pasi) + 1 , &
                       int((ionc_ljmi-ionc_global_jmin)/ionc_pasj) + 1 , &
                       dim_time /)
         endif
      endif
 
      dimi = (ionc_lima - ionc_limi)/ ionc_pasi + 1
      dimj = (ionc_ljma - ionc_ljmi)/ ionc_pasj + 1

      allocate(bid(dimi,dimj))
      count = (/dimi,dimj,1/)
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      if (scale_factor /= 1.e0 .and. present(fill_value)) then
         indi = 1
         do i=ionc_limi,ionc_lima, ionc_pasi
            indj = 1
            do j=ionc_ljmi,ionc_ljma, ionc_pasj
               if (var(i,j) /= fill_value) then
                  bid(indi,indj)=real(nint((var(i,j)- add_offset)/scale_factor))
               else 
                  ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
                  bid(indi,indj)=-HUGE(1_2)-1.0_4
               endif
               indj = indj + 1
            enddo
            indi = indi + 1
         enddo
      else
         indj = 1
         do j=ionc_ljmi,ionc_ljma, ionc_pasj
            indi = 1
            do i=ionc_limi,ionc_lima, ionc_pasi
               bid(indi,indj)=var(i,j)
               indi = indi + 1
            enddo
            indj = indj + 1
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid) 
   endif
   
   return
   end subroutine ionc4_write_xyt_real
   
      
! ***************************************************************
! * subroutine ionc4_write_xyt_int                  *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable entiere 3D = lon, lat, temps    *
! *        a partir d'une variable reelle 2D = lon, lat     *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! *       - l_out_nc4par     :  logique pour I/O parallel  *
! ***************************************************************
   subroutine ionc4_write_xyt_int(nom_fichier,var_name,var,imin,imax,jmin,jmax,nrec,l_out_nc4par)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,nrec
   integer,dimension(imin:imax,jmin:jmax) :: var
   logical,optional,intent(in) :: l_out_nc4par
   
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id, nc_err,dim_time_id,dim_time
   integer :: i,j,indi,indj,dimi,dimj
   integer,DIMENSION(3) :: start,count
   integer,dimension(:,:),allocatable :: bid
!
!******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_xyt_int"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      call ionc4_readbounds(nc_id)
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      if (nrec .eq. 0) then
         nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
         call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         dim_time = nrec
      endif
      
      start  = (/1,1,dim_time/)
      if (present(l_out_nc4par)) then
         if (l_out_nc4par) then
            call ionc4_locbounds(nc_id,imin,imax,jmin,jmax)
            start = (/ int((ionc_limi-ionc_global_imin)/ionc_pasi) + 1 , &
                       int((ionc_ljmi-ionc_global_jmin)/ionc_pasj) + 1 , &
                       dim_time /)
         endif
      endif

      dimi = (ionc_lima - ionc_limi)/ ionc_pasi + 1
      dimj = (ionc_ljma - ionc_ljmi)/ ionc_pasj + 1

      allocate(bid(dimi,dimj))
      count = (/dimi,dimj,1/)
      
      indi = 1
      do i=ionc_limi,ionc_lima, ionc_pasi
         indj = 1
         do j=ionc_ljmi,ionc_ljma, ionc_pasj
            bid(indi,indj)=var(i,j)
            indj = indj + 1
         enddo
         indi = indi + 1
      enddo
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid) 
   endif
   
   return
   end subroutine ionc4_write_xyt_int
   
! ***************************************************************
! * subroutine ionc4_write_xy_double                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 2D = lon, lat    *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - l_out_nc4par     :  logique pour I/O parallel
! ***************************************************************
   subroutine ionc4_write_xy_double(nom_fichier,var_name,var,imin,imax,jmin,jmax,fill_value,l_out_nc4par)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,jmin,jmax
   real(kind=8) ,dimension(imin:imax,jmin:jmax) :: var
   real(kind=8),intent(in),optional :: fill_value
   logical,optional,intent(in) :: l_out_nc4par 
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id, nc_err
   integer :: i,j,indi,indj,dimi,dimj
   real(kind=8) :: add_offset,scale_factor
   integer,DIMENSION(2) :: start,count
   real(kind=8),dimension(:,:),allocatable :: bid
   
!******** FONCTIONS **************
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_xy_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      call ionc4_readbounds(nc_id)
      
      start  = (/1,1/)
      if (present(l_out_nc4par)) then
         if (l_out_nc4par) then
            call ionc4_locbounds(nc_id,imin,imax,jmin,jmax)
            start = (/ int((ionc_limi-ionc_global_imin)/ionc_pasi) + 1 , &
                       int((ionc_ljmi-ionc_global_jmin)/ionc_pasj) + 1 /)
         endif
      endif

      dimi = (ionc_lima - ionc_limi)/ ionc_pasi + 1
      dimj = (ionc_ljma - ionc_ljmi)/ ionc_pasj + 1

      allocate(bid(dimi,dimj))
      count = (/dimi,dimj/)
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.d0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.d0
      endif
      
      if (scale_factor /= 1.d0 .and. present(fill_value)) then      
         indi = 1
         do i=ionc_limi,ionc_lima, ionc_pasi
            indj = 1
            do j=ionc_ljmi,ionc_ljma, ionc_pasj
               if (var(i,j)/= fill_value) then
                  bid(indi,indj)=REAL(nint((var(i,j)- add_offset)/scale_factor),8)
               else
                  ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
                  bid(indi,indj)=-HUGE(1_2)-1.0_8
               endif
               indj = indj + 1
            enddo
            indi = indi + 1
         enddo
      else
         indi = 1
         do i=ionc_limi,ionc_lima, ionc_pasi
            indj = 1
            do j=ionc_ljmi,ionc_ljma, ionc_pasj
               bid(indi,indj)=var(i,j)
               indj = indj + 1
            enddo
            indi = indi + 1
         enddo
      endif

      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid)
   endif
   
   return
   end subroutine ionc4_write_xy_double

! ***************************************************************
! * subroutine ionc4_write_xy_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 2D = lon, lat    *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - l_out_nc4par     :  logique pour I/O parallel
! ***************************************************************
   subroutine ionc4_write_xy_real(nom_fichier,var_name,var,imin,imax,jmin,jmax,fill_value,l_out_nc4par)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,jmin,jmax
   real(kind=4),dimension(imin:imax,jmin:jmax) :: var
   real(kind=4),intent(in),optional :: fill_value
   logical,optional,intent(in) :: l_out_nc4par 
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id, nc_err
   integer :: i,j,indi,indj,dimi,dimj
   real(kind=4) :: add_offset,scale_factor
   integer,DIMENSION(2) :: start,count
   real(kind=4),dimension(:,:),allocatable :: bid
   
!******** FONCTIONS **************
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_xy_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      call ionc4_readbounds(nc_id)

      start  = (/1,1/)
      if (present(l_out_nc4par)) then
         if (l_out_nc4par) then
            call ionc4_locbounds(nc_id,imin,imax,jmin,jmax)
            start = (/ int((ionc_limi-ionc_global_imin)/ionc_pasi) + 1 , &
                       int((ionc_ljmi-ionc_global_jmin)/ionc_pasj) + 1 /)
         endif
      endif

      dimi = (ionc_lima - ionc_limi)/ ionc_pasi + 1
      dimj = (ionc_ljma - ionc_ljmi)/ ionc_pasj + 1

      allocate(bid(dimi,dimj))
      count = (/dimi,dimj/)
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      if (scale_factor /= 1.e0 .and. present(fill_value)) then     
         indi = 1
         do i=ionc_limi,ionc_lima, ionc_pasi
            indj = 1
            do j=ionc_ljmi,ionc_ljma, ionc_pasj
               if (var(i,j) /= fill_value) then
                  bid(indi,indj)=real(nint((var(i,j)- add_offset)/scale_factor))
               else
                  ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
                  bid(indi,indj)=-HUGE(1_2)-1.0_4
               endif
               indj = indj + 1
            enddo
            indi = indi + 1
         enddo
      else
         indi = 1
         do i=ionc_limi,ionc_lima, ionc_pasi
            indj = 1
            do j=ionc_ljmi,ionc_ljma, ionc_pasj
               bid(indi,indj)=var(i,j)
               indj = indj + 1
            enddo
            indi = indi + 1
         enddo
      endif

      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid) 
   endif
   
   return

   end subroutine ionc4_write_xy_real   
   
! ***************************************************************
! * subroutine ionc4_write_sig_double                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 4D = lon, lat, z, temps    *
! *        a partir d'une variable reelle 3D = z , lon, lat     *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - kmax        :  borne max dim k*
! ***************************************************************
   subroutine ionc4_write_sig_double(nom_fichier,var_name,var,kmax)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: kmax
   real(kind=8) ,dimension(0:kmax+1) :: var
   
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id, nc_err
   integer :: k,indk,dimk
   real(kind=4) :: add_offset,scale_factor
   integer,DIMENSION(1) :: start,count
   real(kind=8),dimension(:),allocatable :: bid
   
!******** FONCTIONS **************
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_sig_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      call ionc4_readbounds(nc_id)
      dimk = (ionc_lkma - ionc_lkmi)/ ionc_pask + 1
      allocate(bid(dimk))
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      start  = (/1/)
      count = (/dimk/)
      
      if (scale_factor /= 1.e0 .and. add_offset /= 0.e0) then
         indk = 1
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            bid(indk)=REAL(nint((var(k)- add_offset)/scale_factor),8)
            indk= indk + 1
         enddo
      else
         indk = 1
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            bid(indk)=var(k)
            indk= indk + 1
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid)
   endif
   
   return
   end subroutine ionc4_write_sig_double
   
! ***************************************************************
! * subroutine ionc4_write_sig_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 4D = lon, lat, z, temps    *
! *        a partir d'une variable reelle 3D = z , lon, lat     *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - kmax        :  borne max dim k*
! ***************************************************************
   subroutine ionc4_write_sig_real(nom_fichier,var_name,var,kmax,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: kmax
   real(kind=4),dimension(0:kmax+1) :: var
   real(kind=4),optional,intent(in) :: fill_value
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id, nc_err
   integer :: k,indk,dimk
   real(kind=4) :: add_offset,scale_factor
   integer,DIMENSION(1) :: start,count
   real(kind=4),dimension(:),allocatable :: bid
   
!******** FONCTIONS **************
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_sig_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      call ionc4_readbounds(nc_id)
      dimk = (ionc_lkma - ionc_lkmi)/ ionc_pask + 1
      allocate(bid(dimk))
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      start  = (/1/)
      count = (/dimk/)
      
      if (scale_factor /= 1.e0 .and. present(fill_value)) then
         indk = 1
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            if (var(k) /= fill_value) then
               bid(indk)=real(nint((var(k)- add_offset)/scale_factor))
            else
               ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
               bid(indk)=-HUGE(1_2)-1.0_4
            endif
            indk= indk + 1
         enddo
      else
         indk = 1
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            bid(indk)=var(k)
            indk= indk + 1
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid)
   endif
   
   return
   end subroutine ionc4_write_sig_real
   
! ***************************************************************
! * subroutine ionc4_write_sigw_double                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 4D = lon, lat, z, temps    *
! *        a partir d'une variable reelle 3D = z , lon, lat     *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - kmax        :  borne max dim k*
! ***************************************************************
   subroutine ionc4_write_sigw_double(nom_fichier,var_name,var,kmax)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: kmax
   real(kind=8) ,dimension(0:kmax+1) :: var
   
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id, nc_err
   integer :: k,indk,dimk
   real(kind=4) :: add_offset,scale_factor
   integer,DIMENSION(1) :: start,count
   real(kind=8),dimension(:),allocatable :: bid
   
!******** FONCTIONS **************
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_sig_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      call ionc4_readbounds(nc_id)
      IF (ionc_lkmi==1 .AND. ionc_pask==1) THEN
          dimk = (ionc_lkma - (ionc_lkmi-1))/ ionc_pask + 1
          allocate(bid(0:dimk))
          bid(0)=var(0)
      ELSE
          dimk = (ionc_lkma - ionc_lkmi)/ ionc_pask + 1
          allocate(bid(dimk))
      END IF
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      start  = (/1/)
      count = (/dimk/)
      
      if (scale_factor /= 1.e0 .and. add_offset /= 0.e0) then
         indk = 1
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            bid(indk)=REAL(nint((var(k)- add_offset)/scale_factor),8)
            indk= indk + 1
         enddo
      else
         indk = 1
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            bid(indk)=var(k)
            indk= indk + 1
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid)
   endif
   
   return
   end subroutine ionc4_write_sigw_double
   
! ***************************************************************
! * subroutine ionc4_write_sigw_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 4D = lon, lat, z, temps    *
! *        a partir d'une variable reelle 3D = z , lon, lat     *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - kmax        :  borne max dim k*
! ***************************************************************
   subroutine ionc4_write_sigw_real(nom_fichier,var_name,var,kmax,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: kmax
   real(kind=4),dimension(0:kmax+1) :: var
   real(kind=4),optional,intent(in) :: fill_value
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id, nc_err
   integer :: k,indk,dimk
   real(kind=4) :: add_offset,scale_factor
   integer,DIMENSION(1) :: start,count
   real(kind=4),dimension(:),allocatable :: bid
   
!******** FONCTIONS **************
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_sig_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else

      call ionc4_readbounds(nc_id)
      IF (ionc_lkmi==1 .AND. ionc_pask==1) THEN
          dimk = (ionc_lkma - (ionc_lkmi-1))/ ionc_pask + 1
          allocate(bid(0:dimk))
          bid(0)=var(0)
      ELSE
          dimk = (ionc_lkma - ionc_lkmi)/ ionc_pask + 1
          allocate(bid(dimk))
      END IF

      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)

      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif

      start  = (/1/)
      count = (/dimk/)

      if (scale_factor /= 1.e0 .and. present(fill_value)) then
         indk = 1
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            if (var(k) /= fill_value) then
               bid(indk)=real(nint((var(k)- add_offset)/scale_factor))
            else
               ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
               bid(indk)=-HUGE(1_2)-1.0_4
            endif
            indk= indk + 1
         enddo
      else
         indk = 1
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            bid(indk)=var(k)
            indk= indk + 1
         enddo
      endif

      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid)
   endif

   return
   end subroutine ionc4_write_sigw_real
   
! ***************************************************************
! * SUBROUTINE ionc4_createfile_traj                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    26/09/06                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)          *
!                       30/07/13 par A. Thevenin (CERFACS)      *
! *      *
! * Role : creer le fichier Netcdf qui contiendra les donnees   *
! *        non maillees (trajectoires)                          *
! *         initialise les unites de temps et l'origine          *
! *        aux valeurs par defaut                               *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier : chemin complet du fichier a *
! *      creer*
! *      - dim_traj : nombre de trajectoires             *
! *                         si dim_traj=0, pas de dimension traj*
! *      - dim_time : nombre de pas de temps que le *
! *           fichier contiendra si dim_time =0  *
! *           la dimension time est unlimited    *
! *           si dim_time = -1 pas de dim temps  *
! *      - dim_z   :                                     *
! *                         si dim_z=0, pas de var z -> traj 2D *
! ***************************************************************
   subroutine ionc4_createfile_traj(nom_fichier,dim_traj,dim_time,dim_z, &
                                    l_out_nc4par,comm_active)

#ifdef MPI   
   include 'mpif.h'
#endif

! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier
   integer          :: dim_time,dim_traj,dim_z
   logical, intent(in), optional  :: l_out_nc4par
   integer, intent(in), optional  :: comm_active

! ******** VARIABLES DE TRAVAIL **************
   
   logical :: l_useless
   integer :: communicator
   integer :: nc_id,traj_id,time_id,z_id
   
! ** Declaration des IDS identifiant les dimensions des variables.
   integer :: dim_traj_id,dim_time_id
   
! ** Declarations des tableaux qui contiendront les IDS des dimensions des variables
   integer,dimension(2) :: tab_dim_traj
   integer,dimension(1) :: tab_dim_time
   integer, allocatable, dimension(:) :: tab_traj
   
   integer :: nc_err,indice,cmode,cache_size
   
! ******** FIN DES DECLARATIONS **************
   
   ionc_rout = "ionc4_createfile_traj"
   
! *** Recuperation des dimensions
   cache_size = 10*dim_traj   
   if (dim_z > 0) cache_size = cache_size*dim_z

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      
      ionc_nfich = ionc_nfich + 1
      if (ionc_nfich .gt. ionc_longtabfich) then
         call ionc4_err(ionc_errcreation,ionc_rout, ' ' ,nom_fichier)
         stop
      else
         ionc_nomfich(ionc_nfich) = nom_fichier
#ifdef MPI
         if (present(comm_active)) then
            communicator=comm_active
         else
#ifdef key_oasis
            communicator=MPI_COMM_MARS
#else
            communicator=MPI_COMM_WORLD
#endif
         endif
         if (present(l_out_nc4par) .and. l_out_nc4par) then
            !print*,'Using parallel I/O features !!!'
            cmode = ior(NF90_NETCDF4, NF90_CLASSIC_MODEL)
            cmode = ior(cmode, NF90_MPIIO)
            nc_err = nf90_create(nom_fichier, cmode, nc_id, cache_size=cache_size, &
                                 comm=communicator, info=MPI_INFO_NULL)
         else
            nc_err = nf90_create(nom_fichier, nf90_share, nc_id)
         endif
#else
         if (present(comm_active)) communicator=comm_active  ! useless for portability only
         if (present(l_out_nc4par)) l_useless = l_out_nc4par  ! useless for portability only
         nc_err = nf90_create(nom_fichier, nf90_share, nc_id)
#endif
         call ionc4_err(nc_err,ionc_rout, 'nf90_create',nom_fichier)
         
         ionc_idfich(ionc_nfich) = nc_id
         
         if (dim_traj > 0) then
            nc_err = nf90_def_dim(nc_id,trim(ionc_nomtraj),dim_traj,dim_traj_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',trim(ionc_nomtraj))
         end if
         
         if (dim_time /= -1) then
            nc_err = nf90_def_dim(nc_id,trim(ionc_nomtime),dim_time,dim_time_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_dim',trim(ionc_nomtime))
         end if
         
! ******** On definit les variables  **********
         
         if (dim_traj > 0) then
         
            tab_dim_traj(1) = dim_traj_id
            nc_err = nf90_def_var(nc_id, trim(ionc_nomtraj), nf90_int, (/dim_traj_id/), traj_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomtraj))

            if (dim_time /= -1) then
               tab_dim_time(1) = dim_time_id
               nc_err = nf90_def_var(nc_id, trim(ionc_nomtime), nf90_double, tab_dim_time, time_id)
               call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomtime))
               
               tab_dim_traj(1) = dim_traj_id
               tab_dim_traj(2) = dim_time_id
               if (dim_z > 0) then
                  nc_err = nf90_def_var(nc_id, trim(ionc_nomz), nf90_real, tab_dim_traj, z_id)
                  call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomz))
               end if
#ifdef MPI
               if (present(l_out_nc4par) .and. l_out_nc4par) then
                  nc_err = nf90_var_par_access(nc_id,time_id,nf90_collective)
                  call ionc4_err(nc_err,ionc_rout,'nf90_var_par_access',trim(ionc_nomtime))
               end if
#endif
            else
               tab_dim_traj(1) = dim_traj_id
               if (dim_z > 0) then
                  nc_err = nf90_def_var(nc_id, trim(ionc_nomz), nf90_real, tab_dim_traj, z_id)
                  call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomz))
               end if
            end if
         else
            if (dim_time /= -1) then
               tab_dim_time(1) = dim_time_id
               nc_err = nf90_def_var(nc_id, trim(ionc_nomtime), nf90_double, tab_dim_time, time_id)
               call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomtime))
#ifdef MPI
               if (present(l_out_nc4par)) then
                 if (l_out_nc4par) then
                   nc_err = nf90_var_par_access(nc_id,time_id,nf90_collective)
                   call ionc4_err(nc_err,ionc_rout,'nf90_var_par_access',trim(ionc_nomtime))
                 end if
               end if
#endif
               tab_dim_time(1) = dim_time_id
               if (dim_z > 0) then
                  nc_err = nf90_def_var(nc_id, trim(ionc_nomz), nf90_real, tab_dim_time, z_id)
                  call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomz))
               end if
            else
               tab_dim_time(1) = dim_time_id
               if (dim_z > 0) then
                  nc_err = nf90_def_var(nc_id, trim(ionc_nomz), nf90_real,tab_dim_traj, z_id)
                  call ionc4_err(nc_err,ionc_rout,'nf90_def_var',trim(ionc_nomz))
               end if
            end if
            
         end if
      end if
      
! ******** On definit les attributs des variables  **********
      if (dim_traj > 0) then
         nc_err = nf90_put_att(nc_id, traj_id,'long_name', ionc_nomtraj)
         call ionc4_err(nc_err,ionc_rout, 'nf90_put_att','TRAJ:long_name')
         nc_err = nf90_put_att(nc_id, traj_id,'axis','X')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TRAJ:axis')
      end if

      if (dim_z > 0) then
         !nc_err = nf90_put_att(nc_id, z_id,'units', 'level')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:units')
         
         nc_err = nf90_put_att(nc_id, z_id, 'long_name ', 'sigma level')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:long_name')
         
         !nc_err = nf90_put_att(nc_id, z_id, 'standard_name ', 'z_grid_index')
         !call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:standard_name')
         
         nc_err = nf90_put_att(nc_id, z_id,'axis','Z')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')

         nc_err = nf90_put_att(nc_id,z_id,'c_grid_axis_shift',0.)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','Z:axis')
         
      end if

      if (dim_time /= -1) then
         nc_err = nf90_put_att(nc_id, time_id, 'long_name', ionc_longnamet)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:long_name')

         nc_err = nf90_put_att(nc_id, time_id, 'standard_name','time')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:standard_name')
         
         nc_err = nf90_put_att(nc_id, time_id, 'units',  ionc_unitst)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:units')

         nc_err = nf90_put_att(nc_id, time_id,'axis','T')
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att_text','TIME:axis')

         nc_err = nf90_put_att(nc_id, time_id, 'time_origin', ionc_originet)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:time_origin')
      end if
      
! ******** On peut maintenant quitter le mode definition *****************
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',' ')

      if (dim_traj > 0) then
         allocate(tab_traj(dim_traj))
         do indice=1, dim_traj
            tab_traj(indice) = indice
         enddo
         nc_err = nf90_put_var(nc_id,traj_id,tab_traj)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_var','TRAJ')
         deallocate(tab_traj)
      endif
      
   else
      call ionc4_err(ionc_errcreation,ionc_rout, ' ' ,nom_fichier)
   end if
   
   return
   end subroutine ionc4_createfile_traj

! ***************************************************************
! * subroutine ionc4_createvar_traj_int                  *
! *      *
! * auteur         :    vgarnier                                *
! * org            :    IFREMER                                 *
! * date creation  :    27/09/11                                *
! * derniere modif :                                            *
! *      *
! * Role : creer une variable entiere 1D = nb traj  ou 2D = nb_traj,t *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - units       :  unite        *
! *       - long_name   :  nom long de la variable        *
! *       - fill_value  :  valeur manquante        *
! *       - standard_name   :  nom standard       *
! *       - valid_min   :  valeur min        *
! *       - valid_max   :  valeur max        *
! *       - dims : entier definissant les dimensions
! *                    1 : variable 1D
! *                    2 : variable 2D
! ***************************************************************
   subroutine ionc4_createvar_traj_int(nom_fichier,var_name,units,long_name, &
                                   fill_value,standard_name,valid_min,valid_max,ndims,&
                                   l_out_nc4par)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),intent(in):: nom_fichier,var_name,units,long_name
   integer,intent(in):: fill_value
   character(len=*),intent(in),optional :: standard_name
   integer,intent(in),optional:: valid_min,valid_max
   integer,intent(in),optional :: ndims
   logical,intent(in),optional :: l_out_nc4par

! ******** VARIABLES DE TRAVAIL **************
   logical :: l_useless
   integer :: nc_id,dim_traj_id,dim_time_id,var_id,nc_err,ndim_var,nf_type,i,dimlen
   integer,dimension(:),allocatable :: tab_dim_var,chunksizes
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_createvar_traj_int"

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      if (present(ndims)) then
         ndim_var = ndims
      else
         ndim_var = 2
      endif
      allocate(tab_dim_var(ndim_var))
      allocate(chunksizes(ndim_var))

      nf_type = nf90_int

      nc_err =  nf90_inq_dimid(nc_id,trim(ionc_nomtraj),dim_traj_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtraj))
      tab_dim_var(1) = dim_traj_id

      if (ndim_var == 2) then
         nc_err =  nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
         tab_dim_var(2) = dim_time_id
      endif

#ifdef MPI
      if (present(l_out_nc4par)) then
         if (l_out_nc4par) then
            do i = 1,ndim_var
               nc_err = nf90_inquire_dimension(nc_id, tab_dim_var(i),len=dimlen)
               call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(var_name))
               chunksizes(i)=max(dimlen,1)
            enddo
            nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id,chunksizes=chunksizes)
            nc_err = nf90_var_par_access(nc_id,var_id,nf90_collective)
         else
            nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
         endif
      else
         nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
      endif
#else
      if (present(l_out_nc4par)) l_useless = l_out_nc4par  ! useless for portability only
      nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
#endif

      call ionc4_err(nc_err,ionc_rout,'nf90_def_var', var_name)

! ******** On definit maintenant les attributs de cette variable *************  
      nc_err = nf90_put_att(nc_id, var_id, 'long_name', long_name)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:long_name',var_name)

      nc_err = nf90_put_att(nc_id, var_id, 'units', units)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:units',var_name)

      if (present(standard_name)) then
         nc_err = nf90_put_att(nc_id, var_id, 'standard_name', standard_name)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:standard_name',var_name)
      endif

      if (present(valid_min)) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_min', valid_min)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
      endif

      if (present(valid_max)) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_max', valid_max)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
      endif

      nc_err = nf90_put_att(nc_id,var_id,'_FillValue', fill_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)

      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif

   if (allocated(tab_dim_var)) deallocate(tab_dim_var)
   if (allocated(chunksizes)) deallocate(chunksizes)
   return
   end subroutine ionc4_createvar_traj_int

! ***************************************************************
! * subroutine ionc4_createvar_traj_real                  *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    16/10/06                                *
! * derniere modif :    01/03/2011 RRamel Alyotech              *
! *      *
! * Role : creer une variable reelle 1D = nb traj  ou 2D = nb_traj,t *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - units       :  unite        *
! *       - long_name   :  nom long de la variable        *
! *       - fill_value  :  valeur manquante        *
! *       - standard_name   :  nom standard       *
! *       - valid_min   :  valeur min        *
! *       - valid_max   :  valeur max        *
! *       - dims : entier definissant les dimensions
! *                    1 : variable 1D
! *                    2 : variable 2D
! ***************************************************************
   subroutine ionc4_createvar_traj_real(nom_fichier,var_name,units,long_name, &
                                   fill_value,standard_name,valid_min,valid_max,ndims,&
                                   l_pack,l_out_nc4par)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),intent(in):: nom_fichier,var_name,units,long_name
   real(kind=4),intent(in):: fill_value
   character(len=*),intent(in),optional :: standard_name
   real(kind=4),intent(in),optional:: valid_min,valid_max
   integer,intent(in),optional :: ndims
   logical,intent(in),optional :: l_pack,l_out_nc4par
   
! ******** VARIABLES DE TRAVAIL **************
   logical :: l_useless
   integer :: nc_id,dim_traj_id,dim_time_id,var_id,nc_err,ndim_var,nf_type,i,dimlen
   integer,dimension(:),allocatable :: tab_dim_var,chunksizes
   real(kind=4):: scale_factor,add_offset
   logical :: packing
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_createvar_traj_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      if (present(ndims)) then
         ndim_var = ndims
      else
         ndim_var = 2
      endif   
      allocate(tab_dim_var(ndim_var))
      allocate(chunksizes(ndim_var))
! activation de la compression en short
      packing = .false.
      nf_type = nf90_real
      if (present(l_pack)) then
         if (l_pack .and. &
             present(valid_min) .and. &
             present(valid_max) ) then 
               packing = .true.
               nf_type = nf90_short
         endif
      endif

      nc_err =  nf90_inq_dimid(nc_id,trim(ionc_nomtraj),dim_traj_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtraj))
      tab_dim_var(1) = dim_traj_id
      
      if (ndim_var == 2) then
         nc_err =  nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
         tab_dim_var(2) = dim_time_id
      endif

#ifdef MPI
      if (present(l_out_nc4par)) then
         if (l_out_nc4par) then
            do i = 1,ndim_var
               nc_err = nf90_inquire_dimension(nc_id, tab_dim_var(i),len=dimlen)
               call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(var_name))
               chunksizes(i)=max(dimlen,1)
            enddo
            nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id,chunksizes=chunksizes)
            nc_err = nf90_var_par_access(nc_id,var_id,nf90_collective)
         else
            nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
         endif
      else
         nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
      endif
#else
      if (present(l_out_nc4par)) l_useless = l_out_nc4par  ! useless for portability only
      nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
#endif

      call ionc4_err(nc_err,ionc_rout,'nf90_def_var', var_name)
      
! ******** On definit maintenant les attributs de cette variable *************    
      nc_err = nf90_put_att(nc_id, var_id, 'long_name', long_name)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:long_name',var_name)
      
      nc_err = nf90_put_att(nc_id, var_id, 'units', units)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:units',var_name)
      
      if (present(standard_name)) then
         nc_err = nf90_put_att(nc_id, var_id, 'standard_name', standard_name)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:standard_name',var_name)
      endif
      
! ****** Calcul du scale factor et de l'add_offset *********************
      if ( packing ) then
         scale_factor = (valid_max - valid_min) / (2**nb_bits - 2)
         add_offset   = (valid_max + valid_min) / 2
         nc_err = nf90_put_att(nc_id, var_id, 'scale_factor',scale_factor)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:scale_factor',var_name)
         nc_err = nf90_put_att(nc_id, var_id, 'add_offset',add_offset)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:add_offset',var_name)
      endif

      if (present(valid_min)) then
         if ( packing ) then
            nc_err = nf90_put_att(nc_id,var_id,'valid_min', nint((valid_min-add_offset)/scale_factor,2))
         else
            nc_err = nf90_put_att(nc_id,var_id,'valid_min', valid_min)
         endif
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
      endif

      if (present(valid_max)) then
         if ( packing ) then
            nc_err = nf90_put_att(nc_id,var_id,'valid_max', nint((valid_max-add_offset)/scale_factor,2))
         else
            nc_err = nf90_put_att(nc_id,var_id,'valid_max', valid_max)
         endif
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
      endif

      if ( packing ) then
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue', -HUGE(1_2)-1_2)
      else
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue', fill_value)
      endif
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif

   if (allocated(tab_dim_var)) deallocate(tab_dim_var)
   if (allocated(chunksizes)) deallocate(chunksizes)
   return
   end subroutine ionc4_createvar_traj_real

! ***************************************************************
! * subroutine ionc4_createvar_traj_double                  *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    16/10/06                                *
! * derniere modif :    01/03/2011 RRamel Alyotech              *
! *      *
! * Role : creer une variable reelle 1D = nb traj  ou 2D = nb_traj,t *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - units       :  unite        *
! *       - long_name   :  nom long de la variable        *
! *       - fill_value  :  valeur manquante        *
! *       - standard_name   :  nom standard       *
! *       - valid_min   :  valeur min        *
! *       - valid_max   :  valeur max        *
! *       - dims : entier definissant les dimensions
! *                    1 : variable 1D
! *                    2 : variable 2D
! ***************************************************************
   subroutine ionc4_createvar_traj_double(nom_fichier,var_name,units,long_name, &
                                   fill_value,standard_name,valid_min,valid_max,ndims,&
                                   l_pack,l_out_nc4par)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),intent(in):: nom_fichier,var_name,units,long_name
   real(kind=8),intent(in):: fill_value
   character(len=*),intent(in),optional :: standard_name
   real(kind=8),intent(in),optional:: valid_min,valid_max
   integer,intent(in),optional :: ndims
   logical,intent(in),optional :: l_pack,l_out_nc4par
   
! ******** VARIABLES DE TRAVAIL **************
   logical :: l_useless
   integer :: nc_id,dim_traj_id,dim_time_id,var_id,nc_err,ndim_var,nf_type,i,dimlen
   integer,dimension(:),allocatable :: tab_dim_var,chunksizes
   real(kind=8):: scale_factor,add_offset
   logical :: packing
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_createvar_traj_double"
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',trim(var_name))
      if (present(ndims)) then
         ndim_var = ndims
      else
         ndim_var = 2
      endif   
      allocate(tab_dim_var(ndim_var))
      allocate(chunksizes(ndim_var))
! activation de la compression en short
      packing = .false.
      nf_type = nf90_double
      if (present(l_pack)) then
         if (l_pack .and. &
             present(valid_min) .and. &
             present(valid_max) ) then 
               packing = .true.
               nf_type = nf90_short
         endif
      endif

      nc_err =  nf90_inq_dimid(nc_id,trim(ionc_nomtraj),dim_traj_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtraj))
      tab_dim_var(1) = dim_traj_id

      if (ndim_var == 2) then
         nc_err =  nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
         tab_dim_var(2) = dim_time_id
      endif

#ifdef MPI
      if (present(l_out_nc4par)) then
         if (l_out_nc4par) then
            do i = 1,ndim_var
               nc_err = nf90_inquire_dimension(nc_id, tab_dim_var(i),len=dimlen)
               call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(var_name))
               chunksizes(i)=max(dimlen,1)
            enddo
            nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id,chunksizes=chunksizes)
            nc_err = nf90_var_par_access(nc_id,var_id,nf90_collective)
         else
            nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
         endif
      else
         nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
      endif
#else
      if (present(l_out_nc4par)) l_useless = l_out_nc4par  ! useless for portability only
      nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
#endif

      call ionc4_err(nc_err,ionc_rout,'nf90_def_var', var_name)
      
! ******** On definit maintenant les attributs de cette variable *************    
      nc_err = nf90_put_att(nc_id, var_id, 'long_name', long_name)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:long_name',var_name)
      
      nc_err = nf90_put_att(nc_id, var_id, 'units', units)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:units',var_name)
      
      if (present(standard_name)) then
         nc_err = nf90_put_att(nc_id, var_id, 'standard_name', standard_name)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:standard_name',var_name)
      endif
      
! ****** Calcul du scale factor et de l'add_offset *********************
      if ( packing ) then
         scale_factor = (valid_max - valid_min) / (2**nb_bits - 2)
         add_offset   = (valid_max + valid_min) / 2
         nc_err = nf90_put_att(nc_id, var_id, 'scale_factor',scale_factor)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:scale_factor',var_name)
         nc_err = nf90_put_att(nc_id, var_id, 'add_offset',add_offset)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:add_offset',var_name)
      endif

      if (present(valid_min)) then
         if ( packing ) then
            nc_err = nf90_put_att(nc_id,var_id,'valid_min', nint((valid_min-add_offset)/scale_factor,2))
         else
            nc_err = nf90_put_att(nc_id,var_id,'valid_min', valid_min)
         endif
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
      endif

      if (present(valid_max)) then
         if ( packing ) then
            nc_err = nf90_put_att(nc_id,var_id,'valid_max', nint((valid_max-add_offset)/scale_factor,2))
         else
            nc_err = nf90_put_att(nc_id,var_id,'valid_max', valid_max)
         endif
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
      endif

      if ( packing ) then
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue', -HUGE(1_2)-1_2)
      else
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue', fill_value)
      endif
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif

   if (allocated(tab_dim_var)) deallocate(tab_dim_var)
   if (allocated(chunksizes)) deallocate(chunksizes)

   return
   end subroutine ionc4_createvar_traj_double
   
! ***************************************************************
! * subroutine ionc4_write_traj_int                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    16/10/06                                *
! * derniere modif :   22/12/10 par rramel (ALYOTECH)                   *
! * derniere modif :   24/03/14 par mhuret                   *
! *      *
! * Role : ecrit une variable integer 1D = traj                  *
! *        a partir d'une variable integer 1D = traj             *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! ***************************************************************
   subroutine ionc4_write_traj_int(nom_fichier,var_name,var,imin,imax,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),                   intent(in) :: nom_fichier,var_name
   integer,                            intent(in) :: imin,imax
   integer, dimension(imin:imax),      intent(in) :: var
   integer,                            intent(in), optional :: fill_value

! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,i
   integer(kind=4), dimension(:), allocatable :: bid
   integer,         dimension(1) :: start,count

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_traj_int"

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   endif

   nc_err = nf90_inq_varid(nc_id,var_name,var_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
   
   allocate(bid(imin:imax))

   do i=imin,imax
      bid(i)=var(i)
   enddo
   
   if ( imax >= imin ) then
      start = (/imin/)
      count = (/imax-imin+1/)
   else
      start = (/1/)
      count = (/0/)
   endif
   nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
   call ionc4_err(nc_err,ionc_rout,' nf90_put_var',var_name)
   deallocate(bid)

   end subroutine ionc4_write_traj_int

! ***************************************************************
! * subroutine ionc4_write_traj_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    16/10/06                                *
! * derniere modif :   22/12/10 par rramel (ALYOTECH)                   *
! *      *
! * Role : ecrit une variable reelle 1D = traj                  *
! *        a partir d'une variable reelle 1D = traj             *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! ***************************************************************
   subroutine ionc4_write_traj_real(nom_fichier,var_name,var,imin,imax,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),                   intent(in) :: nom_fichier,var_name
   integer,                            intent(in) :: imin,imax
   real(kind=4), dimension(imin:imax), intent(in) :: var
   real(kind=4),                       intent(in), optional :: fill_value

! ******** VARIABLES DE TRAVAIL **************
   integer      :: nc_id,var_id,nc_err,i
   real(kind=4) :: add_offset,scale_factor
   real(kind=4), dimension(:), allocatable :: bid
   integer,      dimension(1) :: start,count

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_traj_real"

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   endif

   nc_err = nf90_inq_varid(nc_id,var_name,var_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)

   nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
   if (nc_err.ne.nf90_noerr) then
      scale_factor = 1.e0
   endif
   nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
   if (nc_err.ne.nf90_noerr) then
      add_offset = 0.e0
   endif

   allocate(bid(imin:imax))

   if (scale_factor /= 1.e0 .and. present(fill_value)) then
      do i=imin,imax
         if (bid(i) /= fill_value) then
            bid(i)=REAL(nint((var(i)- add_offset)/scale_factor),4)
         else
            ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
            bid(i) = -HUGE(1_2)-1.0_4
         endif
      enddo
   else
      do i=imin,imax
         bid(i)=var(i)
      enddo
   endif

   if ( imax >= imin ) then
      start = (/imin/)
      count = (/imax-imin+1/)
   else
      start = (/1/)
      count = (/0/)
   endif
   nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
   call ionc4_err(nc_err,ionc_rout,' nf90_put_var',var_name)
   deallocate(bid)

   end subroutine ionc4_write_traj_real
   
! ***************************************************************
! * subroutine ionc4_write_traj_double                  *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    16/10/06                                *
! * derniere modif :   22/12/10 par rramel (ALYOTECH)                   *
! *      *
! * Role : ecrit une variable reelle 1D = traj                  *
! *        a partir d'une variable reelle 1D = traj             *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - fill_value  :  valeur_manquante
! ***************************************************************
   subroutine ionc4_write_traj_double(nom_fichier,var_name,var,imin,imax,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),                   intent(in) :: nom_fichier,var_name
   integer,                            intent(in) :: imin,imax
   real(kind=8), dimension(imin:imax), intent(in) :: var
   real(kind=8),                       intent(in), optional :: fill_value

! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,i
   real(kind=8) :: add_offset,scale_factor
   real(kind=8), dimension(:), allocatable :: bid
   integer,      dimension(1) :: start,count

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_traj_double"

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   endif

   nc_err = nf90_inq_varid(nc_id,var_name,var_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)

   nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
   if (nc_err.ne.nf90_noerr) then
      scale_factor = 1.d0
   endif
   nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
   if (nc_err.ne.nf90_noerr) then
      add_offset = 0.d0
   endif

   allocate(bid(imin:imax))

   if (scale_factor /= 1.d0 .and. present(fill_value)) then
      do i=imin,imax
         if (bid(i) /= fill_value) then
            bid(i)=REAL(nint((var(i)- add_offset)/scale_factor),8)
         else
            ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
            bid(i) = -HUGE(1_2)-1.0_8
         endif
      enddo
   else
      do i=imin,imax
         bid(i)=var(i)
      enddo
   endif

   if ( imax >= imin ) then
      start = (/imin/)
      count = (/imax-imin+1/)
   else
      start = (/1/)
      count = (/0/)
   endif
   nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
   call ionc4_err(nc_err,ionc_rout,' nf90_put_var',var_name)
   deallocate(bid)

   end subroutine ionc4_write_traj_double

! ***************************************************************
! * subroutine ionc4_write_trajt_int                    *
! *      *
! * auteur         :    vgarnier                                *
! * org            :    IFREMER                                 *
! * date creation  :    27/09/11                                *
! * derniere modif :                                            *
! *      *
! * Role : ecrit une variable entiere 1D = traj                 *
! *        a partir d'une variable entiere 1D = traj            *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - nrec       :   numero d'enregistrement *
! *                        0 si t est illimite*
! ***************************************************************
   subroutine ionc4_write_trajt_int(nom_fichier,var_name,var,imin,imax,nrec,fill_value)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),intent(in) :: nom_fichier,var_name
   integer,intent(in) :: imin,imax,nrec
   integer,dimension(imin:imax),intent(in) :: var
   integer,intent(in),optional :: fill_value

! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,dim_traj_id,dim_time_id,dim_time,i
   integer,dimension(:),allocatable :: bid
   integer,dimension(2) :: start,count

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_trajt_int"

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   endif

   nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtraj),dim_traj_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','traj')

   nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))

   if (nrec .eq. 0) then
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
   else
      dim_time = nrec
   endif

   nc_err = nf90_inq_varid(nc_id,var_name,var_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)

   allocate(bid(imin:imax))

   do i=imin,imax
      bid(i)=var(i)
   enddo

   if ( imax >= imin ) then
      start = (/imin,dim_time/)
      count = (/imax-imin+1,1/)
   else
      start = (/1,1/)
      count = (/0,0/)
   endif
   nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
   call ionc4_err(nc_err,ionc_rout,' nf90_put_var',var_name)
   deallocate(bid)

   end subroutine ionc4_write_trajt_int


! ***************************************************************
! * subroutine ionc4_write_trajt_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    16/10/06                                *
! * derniere modif :   22/12/10 par rramel (ALYOTECH)                   *
! *      *
! * Role : ecrit une variable reelle 1D = traj                  *
! *        a partir d'une variable reelle 1D = traj             *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - nrec       :   numero d'enregistrement *
! *                        0 si t est illimite*
! ***************************************************************
   subroutine ionc4_write_trajt_real(nom_fichier,var_name,var,imin,imax,nrec,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),intent(in) :: nom_fichier,var_name
   integer,intent(in) :: imin,imax,nrec
   real(kind=4),dimension(imin:imax),intent(in) :: var
   real(kind=4),intent(in),optional :: fill_value

! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,dim_traj_id,dim_time_id,dim_time,i
   real(kind=4) :: add_offset,scale_factor
   real(kind=4),dimension(:),allocatable :: bid
   integer,dimension(2) :: start,count

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_trajt_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   endif

   nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtraj),dim_traj_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','traj')

   nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))

   if (nrec .eq. 0) then
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
   else
      dim_time = nrec
   endif

   nc_err = nf90_inq_varid(nc_id,var_name,var_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)

   nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
   if (nc_err.ne.nf90_noerr) then
      scale_factor = 1.e0
   endif
   nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
   if (nc_err.ne.nf90_noerr) then
      add_offset = 0.e0
   endif

   allocate(bid(imin:imax))

   if (scale_factor /= 1.e0 .and. present(fill_value)) then
      do i=imin,imax
         if (bid(i) /= fill_value) then
            bid(i)=real(nint((var(i)- add_offset)/scale_factor))
         else
            ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
            bid(i) = -HUGE(1_2)-1.0_4
         endif
      enddo
   else
      do i=imin,imax
         bid(i)=var(i)
      enddo
   endif

   if ( imax >= imin ) then
      start = (/imin,dim_time/)
      count = (/imax-imin+1,1/)
   else
      start = (/1,1/)
      count = (/0,0/)
   endif
   nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
   call ionc4_err(nc_err,ionc_rout,' nf90_put_var',var_name)
   deallocate(bid)

   end subroutine ionc4_write_trajt_real
   
! ***************************************************************
! * subroutine ionc4_write_trajt_double                *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    16/10/06                                *
! * derniere modif :   22/12/10 par rramel (ALYOTECH)                   *
! *      *
! * Role : ecrit une variable reelle 1D = traj                  *
! *        a partir d'une variable reelle 1D = traj             *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - nrec       :   numero d'enregistrement *
! *                        0 si t est illimite*
! ***************************************************************
   subroutine ionc4_write_trajt_double(nom_fichier,var_name,var,imin,imax,nrec,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),intent(in) :: nom_fichier,var_name
   integer,intent(in) :: imin,imax,nrec
   real(kind=8),dimension(imin:imax),intent(in) :: var
   real(kind=8),intent(in),optional :: fill_value

! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,dim_traj_id,dim_time_id,dim_time,i
   real(kind=8) :: add_offset,scale_factor
   real(kind=8),dimension(:),allocatable :: bid
   integer,dimension(2) :: start,count

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_trajt_double"

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   endif

   nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtraj),dim_traj_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','traj')
   
   nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
   
   if (nrec .eq. 0) then
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
   else
      dim_time = nrec
   endif
   
   nc_err = nf90_inq_varid(nc_id,var_name,var_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
   
   nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
   if (nc_err.ne.nf90_noerr) then
      scale_factor = 1.d0
   endif
   nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
   if (nc_err.ne.nf90_noerr) then
      add_offset = 0.d0
   endif

   allocate(bid(imin:imax))

   if (scale_factor /= 1.d0 .and. present(fill_value)) then
      do i=imin,imax
         if (bid(i) /= fill_value) then
            bid(i)=real(nint((var(i)- add_offset)/scale_factor))
         else
            ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
            bid(i) = -HUGE(1_2)-1.0_8
         endif
      enddo
   else
      do i=imin,imax
         bid(i)=var(i)
      enddo
   endif

   if ( imax >= imin ) then
      start = (/imin,dim_time/)
      count = (/imax-imin+1,1/)
   else
      start = (/1,1/)
      count = (/0,0/)
   endif
   nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
   call ionc4_err(nc_err,ionc_rout,' nf90_put_var',var_name)
   deallocate(bid)

   end subroutine ionc4_write_trajt_double

! ***************************************************************
! * subroutine ionc4_read_traj_int                             *
! *                                                             *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/01/12                                *
! * derniere modif :    25/01/12 par jfleroux                   *
! * derniere modif :    24/03/14 par mhuret                     *
! *                                                             *
! * Role : lit une variable integer 1D = traj                    *
! *                                                             *
! * Parametres :                                                *
! *       entree:- nom_fichier :  nom du fichier                *
! *              - var_name    :  nom de la variable            *
! *              - imin        :  borne min dim i               *
! *              - imax        :  borne max dim i               *
! *       sortie:- var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_traj_int(nom_fichier,var_name,var,imin,imax)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax
   integer :: var(imin:imax)
   
! ******** DECLARATIONS EXTERNES *************
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,nc_err_fill
   integer, dimension(1) :: start,count
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_traj_int"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   endif

   nc_err = nf90_inq_varid(nc_id,var_name,var_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
 
   start = (/imin/)
   count = (/imax-imin+1/)
   nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
   call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)

   end subroutine ionc4_read_traj_int

! ***************************************************************
! * subroutine ionc4_read_traj_real                             *
! *                                                             *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/01/12                                *
! * derniere modif :    25/01/12 par jfleroux                   *
! *                                                             *
! * Role : lit une variable reelle 1D = traj                    *
! *                                                             *
! * Parametres :                                                *
! *       entree:- nom_fichier :  nom du fichier                *
! *              - var_name    :  nom de la variable            *
! *              - imin        :  borne min dim i               *
! *              - imax        :  borne max dim i               *
! *       sortie:- var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_traj_real(nom_fichier,var_name,var,imin,imax,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax
   real(kind=4) :: var(imin:imax)
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** DECLARATIONS EXTERNES *************
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,dim_traj_id,dim_traj,nc_err_fill
   integer ,dimension(1) :: start,count
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
   real(kind=4),dimension(:),allocatable :: tab_add_offset,tab_scale_factor
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_traj_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   endif

   nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtraj),dim_traj_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','traj')
   
   nc_err = nf90_inquire_dimension(nc_id,dim_traj_id,len=dim_traj)
   call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtraj))
   
   nc_err = nf90_inq_varid(nc_id,var_name,var_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
   
   nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
   if (nc_err.ne.nf90_noerr) then
      scale_factor = 1.e0
   endif
   nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
   if (nc_err.ne.nf90_noerr) then
      add_offset = 0.e0
   endif
   nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
   
   allocate(tab_add_offset(dim_traj))
   allocate(tab_scale_factor(dim_traj))
   tab_add_offset(:) = add_offset
   tab_scale_factor(:)=scale_factor
   
   start = (/imin/)
   count = (/imax-imin+1/)
   nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
   call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
   
   var(:)=var(:)*tab_scale_factor+tab_add_offset
  
   if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
      where (var(:)==fillvalue_read*scale_factor+add_offset ) var(:)=fillvar
   endif

   deallocate(tab_add_offset)
   deallocate(tab_scale_factor)

   end subroutine ionc4_read_traj_real
   
! ***************************************************************
! * subroutine ionc4_read_traj_double                           *
! *                                                             *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/01/12                                *
! * derniere modif :    25/01/12 par jfleroux                   *
! *                                                             *
! * Role : lit une variable double 1D = traj                    *
! *                                                             *
! * Parametres :                                                *
! *       entree:- nom_fichier :  nom du fichier                *
! *              - var_name    :  nom de la variable            *
! *              - imin        :  borne min dim i               *
! *              - imax        :  borne max dim i               *
! *       sortie:- var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_traj_double(nom_fichier,var_name,var,imin,imax,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax
   real(kind=8) :: var(imin:imax)
   real(kind=8),optional,intent(in) :: fillvar
   
! ******** DECLARATIONS EXTERNES *************
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,dim_traj_id,dim_traj,nc_err_fill
   integer, dimension(1) :: start,count
   real(kind=8) :: add_offset,scale_factor,fillvalue_read
   real(kind=8),dimension(:),allocatable :: tab_add_offset,tab_scale_factor

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_traj_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   endif
      
   nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtraj),dim_traj_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','traj')
   
   nc_err = nf90_inquire_dimension(nc_id,dim_traj_id,len=dim_traj)
   call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtraj))

   nc_err = nf90_inq_varid(nc_id,var_name,var_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
   
   nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
   if (nc_err.ne.nf90_noerr) then
      scale_factor = 1.d0
   endif
   nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
   if (nc_err.ne.nf90_noerr) then
      add_offset = 0.d0
   endif
   nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
   
   allocate(tab_add_offset(dim_traj))
   allocate(tab_scale_factor(dim_traj))
   tab_add_offset(:) = add_offset
   tab_scale_factor(:)=scale_factor
   
   start = (/ imin /)
   count = (/ imax-imin+1 /)
   nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
   call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
   
   var(:)=var(:)*tab_scale_factor+tab_add_offset
   
   if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
      where (var(:)==fillvalue_read*scale_factor+add_offset ) var(:)=fillvar
   endif

   deallocate(tab_add_offset)
   deallocate(tab_scale_factor)

   end subroutine ionc4_read_traj_double

! ***************************************************************
! * subroutine ionc4_read_trajt_int                             *
! *                                                             *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/01/12                                *
! * derniere modif :    25/01/12 par jfleroux                   *
! *                                                             *
! * Role : lit une variable double 2D = traj , temps            *
! *                                                             *
! * Parametres :                                                *
! *       entree:- nom_fichier :  nom du fichier                *
! *              - var_name    :  nom de la variable            *
! *              - imin        :  borne min dim i               *
! *              - imax        :  borne max dim i               *
! *              - nrec        :  numero d'enrtegistrement      *
! *       sortie:- var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_trajt_int(nom_fichier,var_name,var,imin,imax,nrec)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,nrec
   integer :: var(imin:imax)
   
! ******** DECLARATIONS EXTERNES *************
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   integer, dimension(2) :: start,count
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_trajt_int"
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   endif
      
   nc_err = nf90_inq_varid(nc_id,var_name,var_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
        
   start = (/ imin,nrec /)
   count = (/ imax-imin+1,1 /)
   nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
   call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)

   end subroutine ionc4_read_trajt_int
   
! ***************************************************************
! * subroutine ionc4_read_trajt_real                            *
! *                                                             *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/01/12                                *
! * derniere modif :    25/01/12 par jfleroux                   *
! *                                                             *
! * Role : lit une variable reelle 2D = traj, temps             *
! *                                                             *
! * Parametres :                                                *
! *       entree:- nom_fichier :  nom du fichier                *
! *              - var_name    :  nom de la variable            *
! *              - imin        :  borne min dim i               *
! *              - imax        :  borne max dim i               *
! *              - nrec        :  numero d'enrtegistrement      *
! *       sortie:- var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_trajt_real(nom_fichier,var_name,var,imin,imax,nrec,fillvar)
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,nrec
   real(kind=4) :: var(imin:imax)
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** DECLARATIONS EXTERNES *************
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,dim_traj_id,dim_traj,nc_err_fill
   integer, dimension(2) :: start,count
   real(kind=4) :: add_offset,scale_factor
   real(kind=4),dimension(:),allocatable :: tab_add_offset,tab_scale_factor,fillvalue_read
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_trajt_real"

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   endif
      
   nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtraj),dim_traj_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','traj')
   
   nc_err = nf90_inquire_dimension(nc_id,dim_traj_id,len=dim_traj)
   call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtraj))

   nc_err = nf90_inq_varid(nc_id,var_name,var_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
   
   nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
   if (nc_err.ne.nf90_noerr) then
      scale_factor = 1.e0
   endif
   nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
   if (nc_err.ne.nf90_noerr) then
      add_offset = 0.e0
   endif
   nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
   
   allocate(tab_add_offset(dim_traj))
   allocate(tab_scale_factor(dim_traj))
   tab_add_offset(:) = add_offset
   tab_scale_factor(:)=scale_factor
   
   start = (/ imin,nrec /)
   count = (/ imax-imin+1,1 /)
   nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
   call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
   
   var(:)=var(:)*tab_scale_factor+tab_add_offset
   
   if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
      where (var(:)==fillvalue_read*scale_factor+add_offset ) var(:)=fillvar
   endif

   deallocate(tab_add_offset)
   deallocate(tab_scale_factor)

   end subroutine ionc4_read_trajt_real

! ***************************************************************
! * subroutine ionc4_read_trajt_double                          *
! *                                                             *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/01/12                                *
! * derniere modif :    25/01/12 par jfleroux                   *
! *                                                             *
! * Role : lit une variable double 2D = traj , temps            *
! *                                                             *
! * Parametres :                                                *
! *       entree:- nom_fichier :  nom du fichier                *
! *              - var_name    :  nom de la variable            *
! *              - imin        :  borne min dim i               *
! *              - imax        :  borne max dim i               *
! *              - nrec        :  numero d'enrtegistrement      *
! *       sortie:- var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_trajt_double(nom_fichier,var_name,var,imin,imax,nrec,fillvar)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,nrec
   real(kind=8) :: var(imin:imax)
   real(kind=8),optional,intent(in) :: fillvar
   
! ******** DECLARATIONS EXTERNES *************
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,dim_traj_id,dim_traj,nc_err_fill
   integer, dimension(2) :: start,count
   real(kind=8) :: add_offset,scale_factor,fillvalue_read
   real(kind=8),dimension(:),allocatable :: tab_add_offset,tab_scale_factor
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_trajt_double"

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   endif
      
   nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtraj),dim_traj_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','traj')
   
   nc_err = nf90_inquire_dimension(nc_id,dim_traj_id,len=dim_traj)
   call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtraj))

   nc_err = nf90_inq_varid(nc_id,var_name,var_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
   
   nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
   if (nc_err.ne.nf90_noerr) then
      scale_factor = 1.d0
   endif
   nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
   if (nc_err.ne.nf90_noerr) then
      add_offset = 0.d0
   endif
   nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
   
   allocate(tab_add_offset(dim_traj))
   allocate(tab_scale_factor(dim_traj))
   tab_add_offset(:) = add_offset
   tab_scale_factor(:)=scale_factor
      
   start = (/ imin,nrec /)
   count = (/ imax-imin+1,1 /)
   nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
   call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
   
   var(:)=var(:)*tab_scale_factor+tab_add_offset
   
   if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
      where (var(:)==fillvalue_read*scale_factor+add_offset ) var(:)=fillvar
   endif

   deallocate(tab_add_offset)
   deallocate(tab_scale_factor)

   end subroutine ionc4_read_trajt_double
   
! ***************************************************************
! * subroutine ionc4_read_dimtraj                  *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/06/07                                *
! * derniere modif :    25/06/07 par jfleroux                   *
! *      *
! * Role : lit le nombre de trajectoire du fichier              *
! *      *
! * Parametres :      *
! *       entree:- nomfich   :  nom du fichier            *
! *       sortie:- dim  :  dimension traj                        *
! ***************************************************************
   subroutine ionc4_read_dimtraj(nom_fichier,dim)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   integer :: dim
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: dim_traj_id,nc_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_dimtraj"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   endif

   nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtraj),dim_traj_id)
   call ionc4_err(nc_err,ionc_rout, 'nf90_inq_dimid',trim(ionc_nomtraj))
   
   nc_err = nf90_inquire_dimension(nc_id,dim_traj_id,len=dim)
   call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtraj))
      
   end subroutine ionc4_read_dimtraj
   
! ***************************************************************
! * subroutine ionc4_read_lon_double                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    09/01/02                                *
! * derniere modif :    09/01/02 par jfleroux                   *
! *      *
! * Role : lit la variable lon                              *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_lon_double(nom_fichier,var,imin,imax)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier
   integer :: imin,imax
   real(kind=8),dimension(imin:imax) ::  var
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer:: nc_id,lon_id,nc_err
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_lon_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,trim(ionc_nomlon),lon_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',trim(ionc_nomlon))
      nc_err = nf90_get_var(nc_id,lon_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_var',trim(ionc_nomlon))
   endif
   
   return
   end subroutine ionc4_read_lon_double
   
! ***************************************************************
! * subroutine ionc4_read_lon_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    09/01/02                                *
! * derniere modif :    09/01/02 par jfleroux                   *
! *      *
! * Role : lit la variable lon                              *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_lon_real(nom_fichier,var,imin,imax)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier
   integer :: imin,imax
   real(kind=4),dimension(imin:imax) ::  var
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer:: nc_id,lon_id,nc_err
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_lon_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,trim(ionc_nomlon),lon_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',trim(ionc_nomlon))
      nc_err = nf90_get_var(nc_id,lon_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_var',trim(ionc_nomlon))
   endif
   
   return
   end subroutine ionc4_read_lon_real
   
! ***************************************************************
! * subroutine ionc4_read_lat_double                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    09/01/02                                *
! * derniere modif :    09/01/02 par jfleroux                   *
! *      *
! * Role : lit la variable lat                              *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_lat_double(nom_fichier,var,imin,imax)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier
   integer :: imin,imax
   real(kind=8),dimension(imin:imax) ::  var
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer:: nc_id,lat_id,nc_err
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_lat_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,trim(ionc_nomlat),lat_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',trim(ionc_nomlat))
      nc_err = nf90_get_var(nc_id,lat_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_var',trim(ionc_nomlat))
   endif
   
   return
   end subroutine ionc4_read_lat_double
   
! ***************************************************************
! * subroutine ionc4_read_lat_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    09/01/02                                *
! * derniere modif :    09/01/02 par jfleroux                   *
! *      *
! * Role : lit la variable lat                              *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_lat_real(nom_fichier,var,imin,imax)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier
   integer :: imin,imax
   real(kind=4),dimension(imin:imax) ::  var
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer:: nc_id,lat_id,nc_err
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_lat_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,trim(ionc_nomlat),lat_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',trim(ionc_nomlat))
      nc_err = nf90_get_var(nc_id,lat_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_var',trim(ionc_nomlat))
   endif
   
   return
   end subroutine ionc4_read_lat_real
   
! ***************************************************************
! * subroutine ionc4_write_int                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    21/02/01 par jfleroux                   *
! *      *
! * Role : ecrit une variable integer               *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! ***************************************************************
   subroutine ionc4_write_int(nom_fichier,var_name,var)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: var
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_int"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_put_var(nc_id,var_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
   endif
   
   return
   end subroutine ionc4_write_int
   
! ***************************************************************
! * subroutine ionc4_read_int                           *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    21/02/01 par jfleroux                   *
! *      *
! * Role : lit une variable entiere dans le fichier               *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       sortie:- var         :  variable                *
! ***************************************************************
   subroutine ionc4_read_int(nom_fichier,var_name,var)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: var
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_int"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      nc_err = nf90_get_var(nc_id,var_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
   endif
   
   return
   end subroutine ionc4_read_int
   
! ***************************************************************
! * subroutine ionc4_write_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    21/02/01 par jfleroux                   *
! *      *
! * Role : ecrit une variable entiere               *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! ***************************************************************
   subroutine ionc4_write_real(nom_fichier,var_name,var)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   real(kind=4) :: var
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_put_var(nc_id,var_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
   endif
   
   return
   end subroutine ionc4_write_real
   
! ***************************************************************
! * subroutine ionc4_read_real                           *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    21/02/01 par jfleroux                   *
! *      * 
! * Role : lit une variable reelle dans le fichier               *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       sortie:- var         :  variable                *
! ***************************************************************
   subroutine ionc4_read_real(nom_fichier,var_name,var)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   real(kind=4) :: var
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      nc_err = nf90_get_var(nc_id,var_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
   endif
   
   return
   end subroutine ionc4_read_real
   
! ***************************************************************
! * subroutine ionc4_write_double                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    21/02/01 par jfleroux                   *
! *      *
! * Role : ecrit une variable entiere               *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! ***************************************************************
   subroutine ionc4_write_double(nom_fichier,var_name,var)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   real(kind=8) :: var
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_put_var(nc_id,var_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
   endif
   
   return
   end subroutine ionc4_write_double
   
! ***************************************************************
! * subroutine ionc4_read_double                           *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    21/02/01 par jfleroux                   *
! *      *
! * Role : lit une variable reelle dans le fichier               *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       sortie:- var         :  variable                *
! ***************************************************************
   subroutine ionc4_read_double(nom_fichier,var_name,var)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   real(kind=8) :: var
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      nc_err = nf90_get_var(nc_id,var_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
   endif
   
   return
   end subroutine ionc4_read_double
   
! ***************************************************************
! * subroutine ionc4_read_time                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/07/01                                *
! * derniere modif :    27/07/01 par jfleroux, ajout            *
! *      *
! * Role : lit le temps              *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - nrec        :  numero d'enregistrement*
! *       sortie:- temps       :  variable temps             *
! ***************************************************************
   subroutine ionc4_read_time(nom_fichier,nrec,temps)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*)  nom_fichier
   real(kind=8) :: temps
   integer :: nrec
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,dim_time_id,dim_time,nc_err
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_time"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      
      if (nrec .gt. dim_time) then
         call ionc4_err(ionc_err_t,ionc_rout,' ',trim(ionc_nomtime))
      endif
      
      nc_err = nf90_inq_varid(nc_id,trim(ionc_nomtime),var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',trim(ionc_nomtime))
      
      nc_err = nf90_get_var(nc_id,var_id,temps,start=(/nrec/))
      call ionc4_err(nc_err,ionc_rout,'nf90_get_var',trim(ionc_nomtime))
      
   endif
   
   return
   end subroutine ionc4_read_time
   
! ***************************************************************
! * fonction ionc4_read_dimt                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    27/07/01                                *
! * derniere modif :                *
! *      *
! * Role : lit la dimension du temps              *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! * Resultat   :      *
! *      sortie:- ionc4_read_dimt                *
! ***************************************************************
   function ionc4_read_dimt(nom_fichier)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   integer :: ionc4_read_dimt
   character(len=*) :: nom_fichier
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,dim_time_id,dim_time,nc_err
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_dimt"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
   endif
   
   ionc4_read_dimt = dim_time
   
   return
   end function ionc4_read_dimt
   
! ***************************************************************
! * subroutine ionc4_write_t_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    21/02/01 par jfleroux                   *
! *      *
! * Role : ecrit une variable reelle 1D = temps        *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! ***************************************************************
   subroutine ionc4_write_t_real(nom_fichier,var_name,var,nrec,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier,var_name
   real(kind=4) :: var
   real(kind=4),intent(in),optional :: fill_value
   integer :: nrec
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,dim_time_id,dim_time,nc_err
   real(kind=4) :: scale_factor,add_offset,bid
   
! ******** FONCTIONS **************
   integer trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_t_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      if (nrec .eq. 0) then
         nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
         call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         dim_time = nrec
      endif
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor', scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset', add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      if (scale_factor /= 1.e0 .and. present(fill_value)) then
         if (var /= fill_value) then
            bid = real(nint((var - add_offset)/scale_factor))
         else
            ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
            bid = -HUGE(1_2)-1.0_4
         endif
      else
         bid = var
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=(/dim_time/))
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',var_name)
   endif
   
   return
   end subroutine ionc4_write_t_real
   
! ***************************************************************
! * subroutine ionc4_read_t_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    21/02/01 par jfleroux                   *
! *      *
! * Role : lit une variable reelle 1D = temps              *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - nrec        :  numero d'enregistrement        *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_t_real(nom_fichier,var_name,var,nrec)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier,var_name
   real(kind=4) :: var
   integer :: nrec
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,dim_time_id,dim_time,nc_err
   real :: scale_factor,add_offset
   real(kind=4) :: bid
   
! ******** FONCTIONS **************
   integer trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_t_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      
      if (nrec .gt. dim_time) then
         call ionc4_err(ionc_err_t,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         nc_err = nf90_inq_varid(nc_id,var_name,var_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
         nc_err = nf90_get_att(nc_id,var_id,'scale_factor', scale_factor)
         if (nc_err.ne.nf90_noerr) then
            scale_factor = 1.e0
         endif
         nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
         if (nc_err.ne.nf90_noerr) then
            add_offset = 0.e0
         endif
         
         nc_err = nf90_get_var(nc_id,var_id,bid,start=(/nrec/))
         call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
         var=bid*scale_factor+add_offset
      endif
   endif
   
   return
   end subroutine ionc4_read_t_real
   
! ***************************************************************
! * subroutine ionc4_write_t_double                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    21/02/01 par jfleroux                   *
! *      *
! * Role : ecrit une variable reelle 1D = temps        *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! ***************************************************************
   subroutine ionc4_write_t_double(nom_fichier,var_name,var,nrec,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier,var_name
   real(kind=8) :: var
   real(kind=8),intent(in),optional :: fill_value
   integer :: nrec
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,dim_time_id,dim_time,nc_err
   real(kind=8) :: bid,scale_factor,add_offset
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_t_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      if (nrec .eq. 0) then
         nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
         call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         dim_time = nrec
      endif
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor', scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.d0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset', add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.d0
      endif
      
      if (scale_factor /= 1.d0 .and. present(fill_value)) then
         if (var /= fill_value) then
            bid = REAL(nint((var - add_offset)/scale_factor),8)
         else
            ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
            bid = -HUGE(1_2)-1.0_8
         endif
      else
         bid = var
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=(/dim_time/))
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',var_name)
   endif
   
   return
   end subroutine ionc4_write_t_double
   
! ***************************************************************
! * subroutine ionc4_read_t_double                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    04/01/11 par rramel (Alyotech)                   *
! *      *
! * Role : lit une variable reelle 1D = temps              *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - nrec        :  numero d'enregistrement        *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_t_double(nom_fichier,var_name,var,nrec)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier,var_name
   real(kind=8) :: var
   integer :: nrec
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,dim_time_id,dim_time,nc_err
   real :: scale_factor,add_offset
   real(kind=8) :: bid
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_t_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      
      if (nrec .gt. dim_time) then
         call ionc4_err(ionc_err_t,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         nc_err = nf90_inq_varid(nc_id,var_name,var_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
         nc_err = nf90_get_att(nc_id,var_id,'scale_factor', scale_factor)
         if (nc_err.ne.nf90_noerr) then
            scale_factor = 1.e0
         endif
         nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
         if (nc_err.ne.nf90_noerr) then
            add_offset = 0.e0
         endif
         
         nc_err = nf90_get_var(nc_id,var_id,bid,start=(/nrec/))
         call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
         var=bid*scale_factor+add_offset
      endif
   endif
   
   return
   end subroutine ionc4_read_t_double
! ***************************************************************
! * subroutine ionc4_init_nomlat            *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    01/03/06                                *
! * derniere modif :    01/03/06 par jfleroux                   *
! *      *
! * Role : initialisation du nom de la dimension j              *
! *      *
! ***************************************************************
   subroutine ionc4_init_nomlat(nomlat)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) ::   nomlat
   
! ******** VARIABLES DE TRAVAIL **************
   
! ******** FIN DES DECLARATIONS **************
   
   ionc_nomlat = nomlat
   
   return
   end subroutine ionc4_init_nomlat
   
   
! ***************************************************************
! * subroutine ionc4_init_nomlon            *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    01/03/06                                *
! * derniere modif :    01/03/06 par jfleroux                   *
! *      *
! * Role : initialisation du nom de la dimension i              *
! *      *
! ***************************************************************
   subroutine ionc4_init_nomlon(nomlon)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) ::   nomlon
   
! ******** VARIABLES DE TRAVAIL **************
   
! ******** FIN DES DECLARATIONS **************
   
   ionc_nomlon = nomlon
   
   return
   end subroutine ionc4_init_nomlon
   
   
! ***************************************************************
! * subroutine ionc4_init_nomz            *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    01/03/06                                *
! * derniere modif :    01/03/06 par jfleroux                   *
! *      *
! * Role : initialisation du nom de la dimension z              *
! *      *
! ***************************************************************
   subroutine ionc4_init_nomz(nomz)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) ::   nomz
   
! ******** VARIABLES DE TRAVAIL **************
   
! ******** FIN DES DECLARATIONS **************
   
   ionc_nomz = nomz
   
   return
   end subroutine ionc4_init_nomz
   
   
! ***************************************************************
! * subroutine ionc4_init_nomtime            *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    01/03/06                                *
! * derniere modif :    01/03/06 par jfleroux                   *
! *      *
! * Role : initialisation du nom de la dimension k              *
! *      *
! ***************************************************************
   subroutine ionc4_init_nomtime(nomtime)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) ::   nomtime
   
! ******** VARIABLES DE TRAVAIL **************
   
! ******** FIN DES DECLARATIONS **************
   
   ionc_nomtime = nomtime
   
   return
   end subroutine ionc4_init_nomtime
   
! ***************************************************************
! * subroutine ionc4_init_nomtraj            *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    27/02/07                                *
! * derniere modif :    27/02/07 par jfleroux                   *
! *      *
! * Role : initialisation du nom de la dimension traj              *
! *      *
! ***************************************************************
   subroutine ionc4_init_nomtraj(nomtraj)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) ::   nomtraj
   
! ******** VARIABLES DE TRAVAIL **************
   
! ******** FIN DES DECLARATIONS **************
   
   ionc_nomtraj = nomtraj
   
   return
   end subroutine ionc4_init_nomtraj
   
! ***************************************************************
! * subroutine ionc4_read_var_dim                  *
! *      *
! * auteur         :    vgarnier                                *
! * org            :    IFREMER                                 *
! * date creation  :    17/06/15                                *
! *                                                             *
! * Role : ajoute l'attribut valid_range a une variable         *
! *      *
! * Parametres :                              *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - dimvar      :  dimensions de la variable            *
! ***************************************************************
   subroutine ionc4_read_var_dim(nom_fichier,var_name,dimvar)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),intent(in) ::   nom_fichier
   character(len=*),intent(in) :: var_name
   integer,intent(out) :: dimvar
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_var_dim"
  
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)

      nc_err = nf90_inquire_variable(nc_id,var_id,ndims=dimvar)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_variable',var_name)      
   endif

   return
   end subroutine ionc4_read_var_dim

! ***************************************************************
! * subroutine ionc4_vatt_range                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    26/02/01                                *
! * derniere modif :    04/01/11 par rramel (Alyotech)                   *
! *      *
! * Role : ajoute l'attribut valid_range a une variable         *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - valid_min   :  valeur valide min              *
! *       - valid_max   :  valeur valide max              *
! ***************************************************************
   subroutine ionc4_vatt_range(nom_fichier,var_name,valid_min,valid_max)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier
   character(len=*) :: var_name
   real(kind=4) :: valid_min, valid_max
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,nc_type
   real(kind=4) ::  range(2),scale_factor,add_offset
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_range"
   
   range(1) = valid_min
   range(2) = valid_max
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      nc_err = nf90_inquire_variable(nc_id,var_id,xtype=nc_type)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_variable',var_name)      

      if (nc_type .eq. nf90_short) then
         if (scale_factor /= 1.e0 .and. add_offset /= 0.e0) then
            nc_err = nf90_put_att(nc_id,var_id,'valid_range', &
                                 (/ nint((valid_min-add_offset)/scale_factor,2),&
                                    nint((valid_max-add_offset)/scale_factor,2) /))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_range',var_name)
         else
            nc_err = nf90_put_att(nc_id,var_id,'valid_range', &
                                           (/nint(valid_min,2),nint(valid_max,2)/))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_range',var_name)
         endif   
      elseif (nc_type .eq. nf90_int) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_range',&
                                           (/nint(valid_min),nint(valid_max)/))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_range',var_name)
      elseif (nc_type .eq. nf90_real) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_range',&
                                           (/valid_min,valid_max/))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_range',var_name)
      elseif (nc_type .eq. nf90_double) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_range',&
                                           (/REAL(valid_min,8),REAL(valid_max,8)/))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_range',var_name)
      endif

      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_range
   
! ***************************************************************
! * subroutine ionc4_vatt_fill_double                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    26/02/01                                *
! * derniere modif :    26/02/01 par jfleroux                   *
! *      *
! * Role : ajoute l'attribut _FillValue a une variable          *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - fillvalue   :  valeur de remplissage          *
! ***************************************************************
   subroutine ionc4_vatt_fill_double(nom_fichier,var_name,fillvalue)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier,var_name
   real(kind=8) :: fillvalue
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,nc_type,var_id,nc_err
   real(kind=4) :: scale_factor,add_offset
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_fill_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_inquire_variable(nc_id,var_id,xtype=nc_type)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_variable',var_name)

      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      if (nc_type .eq. nf90_short) then
         if (scale_factor /= 1.e0 .and. add_offset /= 0.e0) then
            nc_err = nf90_put_att(nc_id,var_id,'_FillValue',nint((fillvalue-add_offset)/scale_factor,2))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)
         else
            nc_err = nf90_put_att(nc_id,var_id,'_FillValue',nint(fillvalue,2))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)
         endif   
      elseif (nc_type .eq. nf90_int) then
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue',nint(fillvalue))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)
      elseif (nc_type .eq. nf90_real) then
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue',real(fillvalue,4))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)
      elseif (nc_type .eq. nf90_double) then
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue',fillvalue)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)
      endif
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_fill_double
   
! ***************************************************************
! * subroutine ionc4_vatt_fill_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    26/02/01                                *
! * derniere modif :    26/02/01 par jfleroux                   *
! *      *
! * Role : ajoute l'attribut _FillValue a une variable          *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - fillvalue   :  valeur de remplissage          *
! ***************************************************************
   subroutine ionc4_vatt_fill_real(nom_fichier,var_name,fillvalue)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier,var_name
   real(kind=4) :: fillvalue
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,nc_type,var_id,nc_err
   real(kind=4) :: scale_factor,add_offset
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_fill_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)

      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      nc_err = nf90_inquire_variable(nc_id,var_id,xtype=nc_type)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      if (nc_type .eq. nf90_short) then
         if (scale_factor /= 1.e0 .and. add_offset /= 0.e0) then
            nc_err = nf90_put_att(nc_id,var_id,'_FillValue',nint((fillvalue-add_offset)/scale_factor,2))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)
         else
            nc_err = nf90_put_att(nc_id,var_id,'_FillValue',nint(fillvalue,2))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)
         endif   
      elseif (nc_type .eq. nf90_int) then
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue',nint(fillvalue))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)
      elseif (nc_type .eq. nf90_real) then
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue',fillvalue)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)
      elseif (nc_type .eq. nf90_double) then
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue',REAL(fillvalue,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)
      endif
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_fill_real
   
! ***************************************************************
! * subroutine ionc4_vatt_fformat                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    26/02/01                                *
! * derniere modif :    26/02/01 par jfleroux                   *
! *      *
! * Role : ajoute l'attribut Fortran_format a une variable      *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - fformat     :  format fortran                 *
! ***************************************************************
   subroutine ionc4_vatt_fformat(nom_fichier,var_name,fformat)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name,fformat
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_fformat"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_put_att(nc_id, var_id, 'FORTRAN_format',  fformat)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:FORTRAN_format',var_name)
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_fformat
   
! ***************************************************************
! * subroutine ionc4_vatt_missing                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    26/02/01                                *
! * derniere modif :    26/02/01 par jfleroux                   *
! *      *
! * Role : ajoute l'attribut missing_value a une variable       *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - mvalue      :  valeur manquante               *
! ***************************************************************
   subroutine ionc4_vatt_missing(nom_fichier,var_name,mvalue)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   real(kind=4) :: mvalue,add_offset,scale_factor
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,nc_type,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_missing"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_inquire_variable(nc_id,var_id,xtype=nc_type)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_variable',var_name)

      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      if (nc_type .eq. nf90_short) then
         if (scale_factor /= 1.e0 .and. add_offset /= 0.e0) then
            nc_err = nf90_put_att(nc_id,var_id,'missing_value',nint((mvalue-add_offset)/scale_factor,2))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:missing_value',var_name)
         else
            nc_err = nf90_put_att(nc_id,var_id,'missing_value',nint(mvalue,2))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:missing_value',var_name)
         endif      
      elseif (nc_type .EQ. nf90_int) then
         nc_err = nf90_put_att(nc_id,var_id,'missing_value',nint(mvalue))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:missing_value',var_name)
      elseif (nc_type .EQ. nf90_real) then
         nc_err = nf90_put_att(nc_id,var_id,'missing_value',mvalue)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:missing_value',var_name)
      elseif (nc_type .EQ. nf90_double) then
         nc_err = nf90_put_att(nc_id,var_id,'missing_value',REAL(mvalue,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:missing_value',var_name)
      endif
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_missing
   
! ***************************************************************
! * subroutine ionc4_write_lon_double                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    09/01/02                                *
! * derniere modif :    09/01/02 par jfleroux                   *
! *      *
! * Role : ecrit la variable lon                      *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! ***************************************************************
   subroutine ionc4_write_lon_double(nom_fichier,var,imin,imax)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   integer :: imin,imax
   real(kind=8) :: var(imin:imax)
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,lon_id,nc_err
   
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_lon_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,trim(ionc_nomlon),lon_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',trim(ionc_nomtime))
      
      nc_err = nf90_put_var(nc_id,lon_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomlon))
   endif
   
   return
   end subroutine ionc4_write_lon_double
   
! ***************************************************************
! * subroutine ionc4_write_lon_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/03/03                                *
! * derniere modif :    25/03/03 par jfleroux                   *
! *      *
! * Role : ecrit la variable lon a partir d'une variable real   =           *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! ***************************************************************
   subroutine ionc4_write_lon_real(nom_fichier,var,imin,imax)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier
   integer :: imin,imax
   real(kind=4) :: var(imin:imax)
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,lon_id,nc_err
   
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_lon_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,trim(ionc_nomlon),lon_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',trim(ionc_nomtime))
      
      nc_err = nf90_put_var(nc_id,lon_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomlon))
   endif
   
   return
   end subroutine ionc4_write_lon_real
   
! ***************************************************************
! * subroutine ionc4_write_lat_double                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    09/01/02                                *
! * derniere modif :    09/01/02 par jfleroux                   *
! *      *
! * Role : ecrit la variable lat                      *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var         :  variable                *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j        *
! ***************************************************************
   subroutine ionc4_write_lat_double(nom_fichier,var,jmin,jmax)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   integer :: jmin,jmax
   real(kind=8) :: var(jmin:jmax)
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,lat_id,nc_err
   
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_lat_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,trim(ionc_nomlat),lat_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',trim(ionc_nomtime))
      nc_err = nf90_put_var(nc_id,lat_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomlat))
   endif
   
   return
   end subroutine ionc4_write_lat_double
   
! ***************************************************************
! * subroutine ionc4_write_lat_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/03/03                                *
! * derniere modif :    25/03/03 par jfleroux                   *
! *      *
! * Role : ecrit la variable lat a partir d'une variable real   *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var         :  variable                *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j        *
! ***************************************************************
   subroutine ionc4_write_lat_real(nom_fichier,var,jmin,jmax)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   integer :: jmin,jmax
   real(kind=4) :: var(jmin:jmax)
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,lat_id,nc_err
   
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_lat_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,trim(ionc_nomlat),lat_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',trim(ionc_nomtime))
      nc_err = nf90_put_var(nc_id,lat_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(ionc_nomlat))
   endif
   
   return
   end subroutine ionc4_write_lat_real

! ***************************************************************
! * subroutine ionc4_read_xy_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    21/02/01 par jfleroux                   *
! *      *
! * Role : lit une variable reelle 2D = lon,lat               *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_xy_real(nom_fichier,var_name,var,imin,imax,jmin,jmax,fillvar)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,jmin,jmax
   real(kind=4),dimension(imin:imax,jmin:jmax) :: var
   real(kind=4),optional,intent(in) :: fillvar

! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,i,j,nc_err_fill
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
      
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_xy_real"

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
        
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
            
      nc_err = nf90_get_var(nc_id,var_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
      do j=jmin,jmax
         do i= imin,imax
            var(i,j)=var(i,j)*scale_factor+add_offset
         end do
      end do

      if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
         where (var(:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:)=fillvar
      endif

   endif

   return
   end subroutine ionc4_read_xy_real

! ***************************************************************
! * subroutine ionc4_read_xy_double                  *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    21/02/01 par jfleroux                   *
! *      *
! * Role : lit une variable double 2D = lon,lat               *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_xy_double(nom_fichier,var_name,var,imin,imax,jmin,jmax,fillvar)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,jmin,jmax
   real(kind=8),dimension(imin:imax,jmin:jmax) :: var
   real(kind=8),optional,intent(in) :: fillvar

! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,i,j,nc_err_fill
   real(kind=8) :: add_offset,scale_factor,fillvalue_read
      
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_xy_double"

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
        
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.d0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.d0
      endif
      nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)

      nc_err = nf90_get_var(nc_id,var_id,var)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
      do j=jmin,jmax
         do i= imin,imax
            var(i,j)=var(i,j)*scale_factor+add_offset
         end do
      end do

      if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
         where (var(:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:)=fillvar
      endif

   endif

   return
   end subroutine ionc4_read_xy_double
   
! ***************************************************************
! * subroutine ionc4_read_subxy_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    15/03/02                                *
! * derniere modif :    15/03/02 par jfleroux                   *
! *      *
! * Role : lit une sous matrice d'une variable reelle 2D        *
! *                                     lon, lat                *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j    *
! *       - imin0        : indice du premier element en i
! *       - jmin0      :  indice du premier element en j
! *     sortie: - var         :  variable                      *
! *
! ***************************************************************
   subroutine ionc4_read_subxy_real(nom_fichier,var_name,var,&
                                    imin,imax,jmin,jmax,imin0,jmin0,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,imin0,jmin0
   real(kind=4) :: var(imin:imax,jmin:jmax)
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,i,j,nc_err_fill
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
   integer,dimension(2) :: start,count
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_subxy_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
      
      start = (/imin+1-imin0,jmin+1-jmin0/)
      count = (/imax - imin + 1,jmax - jmin + 1/)
      
      nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
      do j=jmin,jmax
         do i= imin,imax
            var(i,j)=var(i,j)*scale_factor+add_offset
         end do
      end do

      if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
         where (var(:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:)=fillvar
      endif

   endif
   
   return
   end subroutine ionc4_read_subxy_real
   
! ***************************************************************
! * subroutine ionc4_read_subxy_double                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    15/03/02                                *
! * derniere modif :    15/03/02 par jfleroux                   *
! *      *
! * Role : lit une sous matrice d'une variable double 2D        *
! *                                     lon, lat                *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - imin0       :  indice du premier element en i
! *       - jmin0       :  indice du premier element en j
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_subxy_double(nom_fichier,var_name,var,&
                                      imin,imax,jmin,jmax,imin0,jmin0,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,imin0,jmin0
   real(kind=8) :: var(imin:imax,jmin:jmax)
   real(kind=8),optional,intent(in) :: fillvar
   
! ******** DECLARATIONS EXTERNES *************
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,i,j,nc_err_fill
   real(kind=8) :: add_offset,scale_factor,fillvalue_read
   integer,dimension(2) :: start,count
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_subxy_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.d0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.d0
      endif
      nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
      
      start = (/imin+1-imin0,jmin+1-jmin0/)
      count = (/imax - imin + 1,jmax - jmin + 1/)
      
      nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
      do j=jmin,jmax
         do i= imin,imax
            var(i,j)=var(i,j)*scale_factor+add_offset
         end do
      end do

      if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
         where (var(:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:)=fillvar
      endif
      
   endif
   
   return
   end subroutine ionc4_read_subxy_double
    
! ****************************************************************
! * subroutine ionc4_createvar_sta_real                          *
! *                                                              *
! * auteur         :    pgarreau,jfleroux                        *
! * org            :    IFREMER                                  *
! * date creation  :    xx/xx/xx                                 *
! * derniere modif :    21/12/10 par rramel                      *
! *                     genericite et nf90                       *
! * Role : creer une variable reelle simple precision            *
! *        pour les fichiers stations                            *
! *                                                              *
! * Parametres :                                                 *
! *       entree:- nom_fichier :  nom du fichier                 *
! *       - var_name    :  nom de la variable                    *
! *       - units       :  unite                                 *
! *       - long_name   :  nom long de la variable               *
! *       - dims        :  dimensions de la variable (optionnel) *
! *                        n'importe quelle combinaison de 'xyzt'*
! *                        (eg :'xt', 'xy', 'xyt', 'z', ... )    *
! *                        Si absent creation d'un scalaire      *
! *       - standard_name: nom standard (optionnel)              *
! *       - valid_min    : valeur minimale acceptable (optionnel)*
! *       - valid_max    : valeur maximale acceptable (optionnel)*
! *       - fill_value   : valeur manquante*
! *       - l_pack       : logique pour compression en short   *
! *
! ****************************************************************
   subroutine ionc4_createvar_sta_real(nom_fichier,var_name,units,&
                               long_name,standard_name,  &
                               valid_min,valid_max,fill_value,&
                               dims,l_pack)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),intent(in) :: nom_fichier,var_name,units,long_name,dims
   character(len=*),intent(in),optional :: standard_name
   real(kind=4),intent(in),optional :: valid_min,valid_max
   real(kind=4),intent(in) :: fill_value
   logical,intent(in),optional :: l_pack
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,dim_id,indice
   integer :: var_id,nc_err,ndim_var,i,j,nf_type
   integer,DIMENSION(:),ALLOCATABLE :: tab_dim_var
   logical :: l_dim,packing
   character(len=1),dimension(3) :: tabdim
   character(len=20),dimension(3) :: ionc_nom
   real(kind=4) :: scale_factor,add_offset

! ******** FONCTIONS **************
   integer trim

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_createvar_sta_real"

#ifdef key_gfortran
   tabdim = (/'s','z','t'/)
!cval because size of ionc_nomz different from ionc_nomlat and ionc_nomlon ? see comionc4
   ionc_nom(1) = ionc_nomsta
   ionc_nom(2) = ionc_nomz
   ionc_nom(3) = ionc_nomtime
#else
   tabdim = (/'s','z','t'/)
   ionc_nom = (/ionc_nomsta,ionc_nomz,ionc_nomtime/)
#endif

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nf_type = nf90_real
      packing = .false.
! on teste si on utilise la compression avec scale_fator et add_offset
      if (present(l_pack)) then
         if( l_pack             .and. &
             present(valid_min) .and. &
             present(valid_max)        ) packing = .true.
      endif

      if ( packing ) nf_type = nf90_short

      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')

      ndim_var = len(dims)
      allocate(tab_dim_var(ndim_var))

      indice = 0
      do j =1,size(tabdim)
         l_dim = .false.
         do i=1,ndim_var
            if (dims(i:i) == tabdim(j)) then
               l_dim = .true.
                exit
            endif
         enddo
         if (l_dim) then
            indice = indice + 1
            nc_err =  nf90_inq_dimid(nc_id,trim(ionc_nom(j)),dim_id)
            call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nom(j)))
            tab_dim_var(indice) = dim_id
         endif
      enddo

      nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_def_var',var_name)

! ******** On definit maintenant les attributs de cette variable *************
      nc_err = nf90_put_att(nc_id, var_id, 'long_name', long_name)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:long_name',var_name)
      
      if (present(standard_name)) then
         nc_err = nf90_put_att(nc_id, var_id, 'standard_name',standard_name)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:standard_name',var_name)
      endif
      
      nc_err = nf90_put_att(nc_id, var_id, 'units',units)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:units',var_name)
      
! ****** Calcul du scale factor et de l'add_offset *********************
      if ( packing ) then
         scale_factor = (valid_max - valid_min) / (2**nb_bits - 2)
         add_offset   = (valid_max + valid_min) / 2
         nc_err = nf90_put_att(nc_id, var_id, 'scale_factor',scale_factor)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:scale_factor',var_name)
         nc_err = nf90_put_att(nc_id, var_id, 'add_offset',add_offset)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:add_offset',var_name)
      endif

      if (present(valid_min)) then
         if ( packing ) then
            nc_err = nf90_put_att(nc_id,var_id,'valid_min', nint((valid_min-add_offset)/scale_factor,2))
         else
            nc_err = nf90_put_att(nc_id,var_id,'valid_min', valid_min)
         endif
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
      endif

      if (present(valid_max)) then
         if ( packing ) then
            nc_err = nf90_put_att(nc_id,var_id,'valid_max', nint((valid_max-add_offset)/scale_factor,2))
         else
            nc_err = nf90_put_att(nc_id,var_id,'valid_max', valid_max)
         endif
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
      endif

      if ( packing ) then
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue', -HUGE(1_2)-1_2)
      else
         nc_err = nf90_put_att(nc_id,var_id,'_FillValue', fill_value)
      endif
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:_FillValue',var_name)

      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif

   if (allocated(tab_dim_var)) deallocate(tab_dim_var)
   return
   end subroutine ionc4_createvar_sta_real
   
! ****************************************************************
! * subroutine ionc4_createvar_sta_char                          *
! *                                                              *
! * auteur         :    vgarnier                                 *
! * org            :    IFREMER                                  *
! * date creation  :    2012/06/29                               *
! * derniere modif :                                             *
! * Role : creer une variable character pour les fichiers stations *
! *                                                              *
! * Parametres :                                                 *
! *       entree:- nom_fichier :  nom du fichier                 *
! *       - var_name    :  nom de la variable                    *
! *       - long_name   :  nom long de la variable               *
! *       - dims        :  dimensions de la variable (optionnel) *
! *                        n'importe quelle combinaison de 'xyzt'*
! *                        (eg :'xt', 'xy', 'xyt', 'z', ... )    *
! *                        Si absent creation d'un scalaire      *
! *       - standard_name: nom standard (optionnel)              *
! *
! ****************************************************************
   subroutine ionc4_createvar_sta_char(nom_fichier,var_name,&
                               long_name,dims)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),intent(in) :: nom_fichier,var_name,long_name,dims
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,dim_id,indice
   integer :: var_id,nc_err,ndim_var,i,j,nf_type,dim_len_id
   integer,DIMENSION(:),ALLOCATABLE :: tab_dim_var
   logical :: l_dim
   character(len=1),dimension(1) :: tabdim
   character(len=20),dimension(1) :: ionc_nom

! ******** FONCTIONS **************
   integer trim

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_createvar_sta_char"

   tabdim = (/'s'/)
   ionc_nom = (/ionc_nomsta/)

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nf_type = nf90_char
! on teste si on utilise la compression avec scale_fator et add_offset

      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')

      nc_err = nf90_def_dim(nc_id,'namelen',30,dim_len_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_def_dim','namelen')

      ndim_var = len(dims)
      allocate(tab_dim_var(ndim_var+1))

      nc_err =  nf90_inq_dimid(nc_id,trim(ionc_nom(1)),dim_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nom(1)))

      tab_dim_var(1) = dim_len_id
      tab_dim_var(2) = dim_id

      nc_err = nf90_def_var(nc_id, var_name, nf_type,tab_dim_var, var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_def_var',var_name)

! ******** On definit maintenant les attributs de cette variable *************
      nc_err = nf90_put_att(nc_id, var_id, 'long_name', long_name)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:long_name',var_name)
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif

   if (allocated(tab_dim_var)) deallocate(tab_dim_var)
   return
   end subroutine ionc4_createvar_sta_char
   
! ***************************************************************
! * subroutine ionc4_write_xt_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    24/01/02                                *
! * derniere modif :    24/01/02 par jfleroux                   *
! *      *
! * Role : ecrit une variable reelle 2D = lon,  temps           *
! *        a partir d'une variable reelle 2D = lon              *
! *        et de l'indice du temps                              *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! ***************************************************************
   subroutine ionc4_write_xt_real(nom_fichier,var_name,var,imin,imax,nrec)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,nrec
   real(kind=4) :: var(imin:imax)
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,dim_time_id,dim_time,dimi
   real(kind=4) :: add_offset,scale_factor
   real(kind=4),dimension(:),allocatable :: bid
   integer,dimension(2) :: start,count
   integer :: i,indi
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_xt_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      call ionc4_readbounds(nc_id)
      dimi = (ionc_lima - ionc_limi)/ ionc_pasi + 1
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      if (nrec .eq. 0) then
         nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
         call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         dim_time = nrec
      endif
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      allocate(bid(dimi))
      start = (/1,dim_time/)
      count = (/dimi,1/)
      
      if (scale_factor /= 1.e0 .and. add_offset /= 0.e0) then 
         indi = 1
         do i=ionc_limi,ionc_lima, ionc_pasi
            bid(indi)=real(nint((var(i)- add_offset)/scale_factor))
            indi = indi + 1
         enddo
      else
         indi = 1
         do i=ionc_limi,ionc_lima, ionc_pasi
            bid(indi)=var(i)
            indi = indi + 1
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid)
   endif
   
   return
   end subroutine ionc4_write_xt_real
   
   
! ***************************************************************
! * subroutine ionc4_read_xt_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    28/02/02                                *
! * derniere modif :    28/02/02 par jfleroux                   *
! *      *
! * Role : lit une variable reelle 2D = lon, temps              *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - nrec        :  numero d'enrtegistrement       *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_xt_real(nom_fichier,var_name,var,imin,imax,nrec,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,nrec
   real(kind=4) :: var(imin:imax)
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** DECLARATIONS EXTERNES *************
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id,nc_err,i,nc_err_fill
   integer,dimension(2) :: start,count
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_xt_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
      
      start(1) = 1
      count(1) = imax - imin + 1
      
      start(2)=nrec
      count(2) = 1
      
      nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_var', var_name)
      
      do i=imin,imax
         var(i)=var(i)*scale_factor+add_offset
      end do

      if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
         where (var(:)==fillvalue_read*scale_factor+add_offset ) var(:)=fillvar
      endif

   endif
   
   return
   end subroutine ionc4_read_xt_real
   
! ***************************************************************
! * subroutine ionc4_write_stat_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    24/01/02                                *
! * derniere modif :    24/01/02 par jfleroux                   *
! *      *
! * Role : ecrit une variable reelle 2D = lon,  temps           *
! *        a partir d'une variable reelle 2D = lon              *
! *        et de l'indice du temps                              *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! ***************************************************************
   subroutine ionc4_write_stat_real(nom_fichier,var_name,var,imin,imax,nrec,fill_value)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,nrec
   real(kind=4) :: var(imin:imax)
   real(kind=4),intent(in),optional :: fill_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   integer :: dim_time_id,dim_sta_id
   integer :: dim_time, dim_sta,i
   integer,dimension(2) :: start,count
   real(kind=4) :: add_offset,scale_factor
   real(kind=4),dimension(:),allocatable :: bid
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_stat_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_dimid(nc_id,'sta',dim_sta_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','sta')
      
      nc_err = nf90_inquire_dimension(nc_id,dim_sta_id,len=dim_sta)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension','sta')
      allocate(bid(dim_sta))
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      if (nrec .eq. 0) then
         nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
         call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         dim_time = nrec
      endif
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      start(1) = 1
      start(2) = dim_time
      count(1) = dim_sta
      count(2) = 1
      
      if (scale_factor /= 1.e0 .and. present(fill_value)) then
         do i=1, dim_sta
            if (var(i) /= fill_value) then
               bid(i)=real(nint((var(i)- add_offset)/scale_factor))
            else
               ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
               bid(i)=-HUGE(1_2)-1.0_4
            endif
         enddo
      else
         do i=1, dim_sta
            bid(i)=var(i)
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_vara_real',trim(var_name))
      deallocate(bid) 
   endif
   
   return
   end subroutine ionc4_write_stat_real

! ***************************************************************
! * subroutine ionc4_write_stat_double                          *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    24/01/02                                *
! * derniere modif :    24/01/02 par jfleroux                   *
! *      *
! * Role : ecrit une variable reelle 2D = lon,  temps           *
! *        a partir d'une variable reelle 2D = lon              *
! *        et de l'indice du temps                              *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! ***************************************************************
   subroutine ionc4_write_stat_double(nom_fichier,var_name,var,imin,imax,nrec,fill_value)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,nrec
   real(kind=8) :: var(imin:imax)
   real(kind=8),intent(in),optional :: fill_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   integer :: dim_time_id,dim_sta_id
   integer :: dim_time, dim_sta,i
   integer,dimension(2) :: start,count
   real(kind=8) :: add_offset,scale_factor
   real(kind=8),dimension(:),allocatable :: bid
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_stat_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_dimid(nc_id,'sta',dim_sta_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','sta')
      
      nc_err = nf90_inquire_dimension(nc_id,dim_sta_id,len=dim_sta)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension','sta')
      allocate(bid(dim_sta))
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      if (nrec .eq. 0) then
         nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
         call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         dim_time = nrec
      endif
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.d0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.d0
      endif
      
      start(1) = 1
      start(2) = dim_time
      count(1) = dim_sta
      count(2) = 1
      
      if (scale_factor /= 1.d0 .and. present(fill_value)) then
         do i=1, dim_sta
            if (var(i) /= fill_value) then
               bid(i)=REAL(nint((var(i)- add_offset)/scale_factor),8)
            else
               ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
               bid(i)=-HUGE(1_2)-1.0_8
            endif
         enddo
      else
         do i=1, dim_sta
            bid(i)=var(i)
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid) 
   endif
   
   return
   end subroutine ionc4_write_stat_double
      
! ***************************************************************
! * subroutine ionc4_read_stat_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    28/02/02                                *
! * derniere modif :    28/02/02 par jfleroux                   *
! *      *
! * Role : lit une variable reelle 2D = sta, temps              *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - nrec        :  numero d'enrtegistrement       *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_stat_real(nom_fichier,var_name,var,imin,imax,nrec,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,nrec
   real(kind=4) :: var(imin:imax)
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** DECLARATIONS EXTERNES *************
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,i,nc_err_fill
   integer ,dimension(2) :: start,count
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_stat_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
      
      start(1) = 1
      count(1) = imax - imin + 1
      start(2)=nrec
      count(2) = 1
      
      nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_vara_real',var_name)
      
      do i=imin,imax
         var(i)=var(i)*scale_factor+add_offset
      end do

      if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
         where (var(:)==fillvalue_read*scale_factor+add_offset ) var(:)=fillvar
      endif

   endif
   
   return
   end subroutine ionc4_read_stat_real
   
! ***************************************************************
! * subroutine ionc4_read_stat_double                  *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    28/02/02                                *
! * derniere modif :    28/02/02 par jfleroux                   *
! *      *
! * Role : lit une variable reelle 2D = sta, temps              *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - nrec        :  numero d'enrtegistrement       *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_stat_double(nom_fichier,var_name,var,imin,imax,nrec,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,nrec
   real(kind=8) :: var(imin:imax)
   real(kind=8),optional,intent(in) :: fillvar
   
! ******** DECLARATIONS EXTERNES *************
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,i,nc_err_fill
   integer ,dimension(2) :: start,count
   real(kind=8) :: add_offset,scale_factor,fillvalue_read
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_stat_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.d0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.d0
      endif
      nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
      
      start(1) = 1
      count(1) = imax - imin + 1
      start(2)=nrec
      count(2) = 1
      
      nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_vara_real',var_name)
      
      do i=imin,imax
         var(i)=var(i)*scale_factor+add_offset
      end do

      if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
         where (var(:)==fillvalue_read*scale_factor+add_offset ) var(:)=fillvar
      endif

   endif
 
   return
   end subroutine ionc4_read_stat_double
   
! ***************************************************************
! * subroutine ionc4_write_sta_char                   *
! *      *
! * auteur         :    vgarnier                               *
! * org            :    IFREMER                                *
! * date creation  :    2012/06/29                             *
! * derniere modif :                                           *
! *      *
! * Role : ecrit une variable chain 1D = nb stations           *
! *        a partir d'une variable chain 1D = nom des stations *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! ***************************************************************
   subroutine ionc4_write_sta_char(nom_fichier,var_name,var,imin,imax)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax
   character(len=*) :: var(imin:imax)
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   integer :: dim_sta_id,dim_len_id,dim_sta
   real(kind=4) :: add_offset,scale_factor
   integer :: i
   character(len=30),dimension(:),allocatable :: bid
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_sta_char"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_dimid(nc_id,'sta',dim_sta_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','sta')
   
      nc_err = nf90_inquire_dimension(nc_id,dim_sta_id,len=dim_sta)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension','sta')
      
      allocate(bid(dim_sta))
      do i=1, dim_sta
         bid(i)=TRIM(var(i))
      enddo
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
     
      nc_err = nf90_put_var(nc_id,var_id,bid)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',TRIM(var_name))
      deallocate(bid) 
   endif
   
   return
   end subroutine ionc4_write_sta_char

! ***************************************************************
! * subroutine ionc4_write_sta_real                   *
! *      *
! * auteur         :    vgarnier                                *
! * org            :    IFREMER                                 *
! * date creation  :    24/01/02                                *
! * derniere modif :    24/01/02 par jfleroux                   *
! *      *
! * Role : ecrit une variable reelle 1D = nb stations           *
! *        a partir d'une variable reelle 1D = nb stations      *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! ***************************************************************
   subroutine ionc4_write_sta_real(nom_fichier,var_name,var,imin,imax,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax
   real(kind=4) :: var(imin:imax)
   real(kind=4),intent(in),optional :: fill_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   integer :: dim_sta_id,dim_sta
   real(kind=4) :: add_offset,scale_factor
   integer :: i
   real(kind=4),dimension(:),allocatable :: bid
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_sta_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_dimid(nc_id,'sta',dim_sta_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','sta')
      
      nc_err = nf90_inquire_dimension(nc_id,dim_sta_id,len=dim_sta)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension','sta')
      allocate(bid(dim_sta))
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      if (scale_factor /= 1.e0 .and. present(fill_value)) then
         do i=1, dim_sta
            if (var(i) /= fill_value) then
               bid(i)=real(nint((var(i)- add_offset)/scale_factor))
            else
               ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
               bid(i)=-HUGE(1_2)-1.0_4
            endif
         enddo
      else
         do i=1, dim_sta
            bid(i)=var(i)
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid) 
   endif
   
   return
   end subroutine ionc4_write_sta_real

! ***************************************************************
! * subroutine ionc4_write_sta_double                   *
! *      *
! * auteur         :    vgarnier                                *
! * org            :    IFREMER                                 *
! * date creation  :    24/01/02                                *
! * derniere modif :    24/01/02 par jfleroux                   *
! *      *
! * Role : ecrit une variable reelle 1D = nb stations           *
! *        a partir d'une variable reelle 1D = nb stations      *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! ***************************************************************
   subroutine ionc4_write_sta_double(nom_fichier,var_name,var,imin,imax,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax
   real(kind=8) :: var(imin:imax)
   real(kind=8),intent(in),optional :: fill_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   integer :: dim_sta_id,dim_sta
   real(kind=8) :: add_offset,scale_factor
   integer :: i
   real(kind=8),dimension(:),allocatable :: bid
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_sta_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_dimid(nc_id,'sta',dim_sta_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','sta')
      
      nc_err = nf90_inquire_dimension(nc_id,dim_sta_id,len=dim_sta)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension','sta')
      allocate(bid(dim_sta))
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.d0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.d0
      endif
      
      if (scale_factor /= 1.d0 .and. present(fill_value)) then
         do i=1, dim_sta
            if (var(i) /= fill_value) then
               bid(i)=REAL(nint((var(i)- add_offset)/scale_factor),8)
            else
               ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
               bid(i)=-HUGE(1_2)-1.0_8
            endif
         enddo
      else
         do i=1, dim_sta
            bid(i)=var(i)
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid) 
   endif
   
   return
   end subroutine ionc4_write_sta_double
   
! ***************************************************************
! * subroutine ionc4_write_zstat_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 3D = nb stations, z, temps *
! *        a partir d'une variable reelle 2D = z , nb stations  *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - kmax        :  borne max dim k*
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! ***************************************************************
   subroutine ionc4_write_zstat_real(nom_fichier,var_name,var,imin,imax,kmax,nrec,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,kmax,nrec
   real(kind=4) :: var(kmax,imin:imax)
   real(kind=4),intent(in),optional :: fill_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   integer :: dim_time_id,dim_sta_id,dim_k_id
   integer :: dim_time,dim_sta,dim_k
   real(kind=4) :: add_offset,scale_factor
   integer,dimension(3) :: start,count
   integer :: i,k
   real(kind=4),dimension(:,:),allocatable :: bid
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_zstat_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,'sta',dim_sta_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','sta')
      
      nc_err = nf90_inquire_dimension(nc_id,dim_sta_id,len=dim_sta)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension','sta')
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomz),dim_k_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomz))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_k_id,len=dim_k)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomz))
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      if (nrec .eq. 0) then
         nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
         call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         dim_time = nrec
      endif
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      allocate(bid(dim_sta,dim_k))
      start  = (/1,1,dim_time/)
      count = (/dim_sta,dim_k,1/)
      
      if (scale_factor /= 1.e0 .and. present(fill_value)) then
         do k=1,kmax
            do i=1,dim_sta
               if (var(k,i) /= fill_value) then
                  bid(i,k)=real(nint((var(k,i)- add_offset)/scale_factor))
               else
                  ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
                  bid(i,k)=-HUGE(1_2)-1.0_4
               endif
            enddo
         enddo
      else
         do k=1,kmax
            do i=1,dim_sta
               bid(i,k)=var(k,i)
            enddo
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid) 
   endif
   
   return
   end subroutine ionc4_write_zstat_real

! ***************************************************************
! * subroutine ionc4_write_zstat_double                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 3D = nb stations, z, temps *
! *        a partir d'une variable reelle 2D = z , nb stations  *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - kmax        :  borne max dim k*
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! ***************************************************************
   subroutine ionc4_write_zstat_double(nom_fichier,var_name,var,imin,imax,kmax,nrec,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,kmax,nrec
   real(kind=8) :: var(kmax,imin:imax)
   real(kind=8),intent(in),optional :: fill_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   integer :: dim_time_id,dim_sta_id,dim_k_id
   integer :: dim_time,dim_sta,dim_k
   real(kind=8) :: add_offset,scale_factor
   integer,dimension(3) :: start,count
   integer :: i,k
   real(kind=8),dimension(:,:),allocatable :: bid
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_zstat_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,'sta',dim_sta_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','sta')
      
      nc_err = nf90_inquire_dimension(nc_id,dim_sta_id,len=dim_sta)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension','sta')
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomz),dim_k_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomz))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_k_id,len=dim_k)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomz))
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      if (nrec .eq. 0) then
         nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
         call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         dim_time = nrec
      endif
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.d0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.d0
      endif
      
      allocate(bid(dim_sta,dim_k))
      start  = (/1,1,dim_time/)
      count = (/dim_sta,dim_k,1/)
      
      if (scale_factor /= 1.d0 .and. present(fill_value)) then
         do k=1,kmax
            do i=1,dim_sta
               if (var(k,i) /= fill_value) then
                  bid(i,k)=REAL(nint((var(k,i)- add_offset)/scale_factor),8)
               else
                  ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
                  bid(i,k)=-HUGE(1_2)-1.0_8
               endif
            enddo
         enddo
      else
         do k=1,kmax
            do i=1,dim_sta
               bid(i,k)=var(k,i)
            enddo
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid) 
   endif
   
   return
   end subroutine ionc4_write_zstat_double
   
! ***************************************************************
! * subroutine ionc4_read_zstat_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : lit une variable reelle 3D = nb stations, z, temps   *
! *        dans une variable reelle 2D = z , nb stations        *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - kmax        :  borne max dim k*
! *       - nrec        :  numero d'enregistrement*
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_zstat_real(nom_fichier,var_name,var,imin,imax,kmax,nrec,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,kmax,nrec
   real(kind=4) :: var(kmax,imin:imax)
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,dim_time_id,dim_time
   integer :: var_id,nc_err,nc_err_fill
   integer :: i,k,indi
   integer,dimension(3) :: start,count
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
   real(kind=4),dimension(imax - imin + 1, kmax) :: bid
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_zstat_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      
      if (nrec .gt. dim_time) then
         call ionc4_err(ionc_err_t,ionc_rout,' ',trim(ionc_nomtime))
      else
         nc_err = nf90_inq_varid(nc_id,var_name,var_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
         nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
         if (nc_err.ne.nf90_noerr) then
            scale_factor = 1.e0
         endif
         nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
         if (nc_err.ne.nf90_noerr) then
            add_offset = 0.e0
         endif
         nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
         
         start(1)  = 1
         count(1) = imax - imin + 1
         start(2)  = 1
         count(2) = kmax
         start(3)  = nrec
         count(3) = 1
         
         nc_err = nf90_get_var(nc_id,var_id,bid,start=start,count=count)
         
         do k = 1,kmax
            indi = 1
            do i=imin,imax
               var(k,i)=bid(indi,k)*scale_factor+add_offset
               indi = indi + 1
            end do
         end do

         if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
            where (var(:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:)=fillvar
         endif
         
      endif
   endif
   
   return
   end subroutine ionc4_read_zstat_real
   
! ***************************************************************
! * subroutine ionc4_write_zsta_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    04/10/06                                *
! * derniere modif :    04/10/06 par jfleroux                   *
! *      *
! * Role : ecrit une variable reelle 2D = nb stations, z        *
! *        a partir d'une variable reelle 2D = z , nbstation    *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - kmax        :  borne max dim k*
! ***************************************************************
   subroutine ionc4_write_zsta_real(nom_fichier,var_name,var,imin,imax,kmax,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,kmax
   real(kind=4) :: var(kmax,imin:imax)
   real(kind=4),intent(in),optional :: fill_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   integer :: dim_sta_id,dim_k_id,dim_sta,dim_k
   real(kind=4) :: add_offset,scale_factor
   integer ,dimension(2):: start,count
   integer :: i,k
   real(kind=4) ,dimension(:,:),allocatable :: bid
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_zsta_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,'sta',dim_sta_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','sta')
      
      nc_err = nf90_inquire_dimension(nc_id,dim_sta_id,len=dim_sta)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension','sta')
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomz),dim_k_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomz))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_k_id,len=dim_k)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomz))
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      allocate(bid(dim_sta,dim_k))
      start = (/1,1/)
      count = (/dim_sta,dim_k/)
      
      if (scale_factor /= 1.e0 .and. present(fill_value)) then
         do k=1,kmax
            do i=1,dim_sta
               if (var(k,i) /= fill_value) then
                  bid(i,k)=real(nint((var(k,i)- add_offset)/scale_factor))
               else
                  ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
                  bid(i,k)=-HUGE(1_2)-1.0_4
               endif
            enddo
         enddo
      else
         do k=1,kmax
            do i=1,dim_sta
               bid(i,k)=var(k,i)
            enddo
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid) 
   endif
   
   return
   end subroutine ionc4_write_zsta_real

! ***************************************************************
! * subroutine ionc4_write_zsta_double                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    04/10/06                                *
! * derniere modif :    04/10/06 par jfleroux                   *
! *      *
! * Role : ecrit une variable reelle 2D = nb stations, z        *
! *        a partir d'une variable reelle 2D = z , nbstation    *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - kmax        :  borne max dim k*
! ***************************************************************
   subroutine ionc4_write_zsta_double(nom_fichier,var_name,var,imin,imax,kmax,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,kmax
   real(kind=8) :: var(kmax,imin:imax)
   real(kind=8),intent(in),optional :: fill_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   integer :: dim_sta_id,dim_k_id,dim_sta,dim_k
   real(kind=8) :: add_offset,scale_factor
   integer ,dimension(2):: start,count
   integer :: i,k
   real(kind=8) ,dimension(:,:),allocatable :: bid
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_zsta_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,'sta',dim_sta_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid','sta')
      
      nc_err = nf90_inquire_dimension(nc_id,dim_sta_id,len=dim_sta)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension','sta')
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomz),dim_k_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomz))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_k_id,len=dim_k)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomz))
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.d0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.d0
      endif
      
      allocate(bid(dim_sta,dim_k))
      start = (/1,1/)
      count = (/dim_sta,dim_k/)
      
      if (scale_factor /= 1.d0 .and. present(fill_value)) then
         do k=1,kmax
            do i=1,dim_sta
               if (var(k,i) /= fill_value) then
                  bid(i,k)=REAL(nint((var(k,i)-add_offset)/scale_factor),8)
               else
                  ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
                  bid(i,k)=-HUGE(1_2)-1.0_8
               endif
            enddo
         enddo
      else
         do k=1,kmax
            do i=1,dim_sta
               bid(i,k)=var(k,i)
            enddo
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid) 
   endif
   
   return
   end subroutine ionc4_write_zsta_double
   
! ***************************************************************
! * subroutine ionc4_read_zsta_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    04/10/06                                *
! * derniere modif :    04/10/06 par jfleroux                   *
! *      *
! * Role : lit une variable reelle 2D = nb stations, z          *
! *        dans une variable reelle 2D = z , nb stations        *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - kmax        :  borne max dim k*
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_zsta_real(nom_fichier,var_name,var,imin,imax,kmax,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,kmax
   real(kind=4) :: var(kmax,imin:imax)
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,nc_err_fill
   integer,dimension(2) :: start,count
   integer :: i,k,indi
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
   real(kind=4),dimension(imax - imin + 1, kmax) :: bid
   
! ******** FONCTIONS **************
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_zsta_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
      
      start(1) = 1
      count(1) = imax - imin + 1
      start(2)= 1
      count(2) = kmax
      
      nc_err = nf90_get_var(nc_id,var_id,bid,start=start,count=count)
      
      do k = 1,kmax
         indi = 1
         do i=imin,imax
            var(k,i)=bid(indi,k)*scale_factor+add_offset
            indi = indi + 1
         end do
      end do
      if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
         where (var(:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:)=fillvar
      endif
      
   endif
   
   return
   end subroutine ionc4_read_zsta_real
   
! ***************************************************************
! * function ionc4_fnlsta                                        *
! *                                                             *
! * auteur         :    vgarnier                                *
! * org            :    IFREMER                                 *
! * date creation  :    06/05/02                                *
! * derniere modif :    16/05/02 par jfleroux                   *
! *                                                             *
! * Role : lit l'indice max des stations                        *
! *                                                             *
! * Parametres :                                                *
! *     entree: - nom_fichier :  nom du fichier                 *
! * Resultat :                                                  *
! *     indice nombre de stations max                           *
! ***************************************************************
   function ionc4_fnlsta(nom_fichier)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   integer :: ionc4_fnlsta 
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: dim_sta_id,lima,nc_err,nc_id
   
! ******** FIN DES DECLARATIONS **************
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,'ionc4_fnlsta', ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_dimid(nc_id,'sta',dim_sta_id)
      call ionc4_err(nc_err,'ionc4_fnlsta', 'nf90_inq_dimid','sta')
      
      nc_err = nf90_inquire_dimension(nc_id,dim_sta_id,len=lima)
      call ionc4_err(nc_err,'ionc4_fnlsta','nf90_inquire_dimension','sta')
      
      ionc4_fnlsta = lima
   end if
   end function ionc4_fnlsta
   
! ***************************************************************
! * subroutine ionc4_read_xyt_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : lit une variable reelle 3D = lon, lat, temps         *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - nrec        :  numero d'enrtegistrement       *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_xyt_real(nom_fichier,var_name,var,imin,imax,jmin,jmax,nrec,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,nrec
   real(kind=4),dimension(imin:imax,jmin:jmax) :: var
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,dim_time_id,dim_time,nc_err_fill
   integer :: i,j
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
   integer,dimension(3) :: start,count
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_xyt_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      
      if (nrec .gt. dim_time) then
         call ionc4_err(ionc_err_t,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         
         nc_err = nf90_inq_varid(nc_id,var_name,var_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
         nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
         if (nc_err.ne.nf90_noerr) then
            scale_factor = 1.e0
         endif
         nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
         if (nc_err.ne.nf90_noerr) then
            add_offset = 0.e0
         endif
         nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
         
         start  = (/1, 1, nrec/)
         count = (/imax - imin + 1, jmax - jmin + 1, 1/)
                 
         nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
         call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
         
         do j=jmin,jmax
            do i=imin,imax
               var(i,j)=var(i,j)*scale_factor+add_offset
            end do
         end do

         if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
            where (var(:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:)=fillvar
         endif

      endif
   endif
   
   return
   end subroutine ionc4_read_xyt_real
   
! ***************************************************************
! * subroutine ionc4_read_xyt_double                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : lit une variable reelle 3D = lon, lat, temps         *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - nrec        :  numero d'enrtegistrement       *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_xyt_double(nom_fichier,var_name,var,imin,imax,jmin,jmax,nrec,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,nrec
   real(kind=8),dimension(imin:imax,jmin:jmax) :: var
   real(kind=8),optional,intent(in) :: fillvar
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,dim_time_id,dim_time,nc_err_fill
   integer :: i,j
   real(kind=8) :: add_offset,scale_factor,fillvalue_read
   integer,dimension(3) :: start,count
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_xyt_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      
      if (nrec .gt. dim_time) then
         call ionc4_err(ionc_err_t,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         
         nc_err = nf90_inq_varid(nc_id,var_name,var_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
         nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
         if (nc_err.ne.nf90_noerr) then
            scale_factor = 1.d0
         endif
         nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
         if (nc_err.ne.nf90_noerr) then
            add_offset = 0.d0
         endif
         nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
         
         start  = (/1, 1, nrec/)
         count = (/imax - imin + 1, jmax - jmin + 1, 1/)
         
         nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
         call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
         
         do j=jmin,jmax
            do i=imin,imax
               var(i,j)=var(i,j)*scale_factor+add_offset
            end do
         end do

         if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
            where (var(:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:)=fillvar
         endif

      endif
   endif
   
   return
   end subroutine ionc4_read_xyt_double

! ***************************************************************
! * subroutine ionc4_read_xyt_int                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : lit une variable entier 3D = lon, lat, temps         *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - nrec        :  numero d'enrtegistrement       *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_xyt_int(nom_fichier,var_name,var,imin,imax,jmin,jmax,nrec)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,nrec
   integer,dimension(imin:imax,jmin:jmax) :: var
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,dim_time_id,dim_time
   integer,dimension(3) :: start,count
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_xyt_int"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      
      if (nrec .gt. dim_time) then
         call ionc4_err(ionc_err_t,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         
         nc_err = nf90_inq_varid(nc_id,var_name,var_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
         start  = (/1, 1, nrec/)
         count = (/imax - imin + 1, jmax - jmin + 1, 1/)
         
         nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
         call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
         
      endif
   endif
   
   return
   end subroutine ionc4_read_xyt_int
   
! ***************************************************************
! * subroutine ionc4_read_subxyt_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    15/03/02                                *
! * derniere modif :    15/03/02 par jfleroux                   *
! *      *
! * Role : lit une sous matrice d'une variable reelle 3D        *
! *                                     lon, lat, temps         *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - nrec        :  numero d'enrtegistrement       *
! *       - imin0        : indice du premier element en i
! *       - jmin0      :  indice du premier element en j
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_subxyt_real(nom_fichier,var_name,var,&
                                     imin,imax,jmin,jmax,nrec,imin0,jmin0,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,imin0,jmin0
   real(kind=4),dimension(imin:imax,jmin:jmax) :: var
   integer :: nrec
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,nc_err_fill
   integer :: dim_time_id,dim_time,i,j
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
   integer,dimension(3) :: start,count
   
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_subxyt_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      
      if (nrec .gt. dim_time) then
         call ionc4_err(ionc_err_t,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         
         nc_err = nf90_inq_varid(nc_id,var_name,var_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
         nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
         if (nc_err.ne.nf90_noerr) then
            scale_factor = 1.e0
         endif
         nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
         if (nc_err.ne.nf90_noerr) then
            add_offset = 0.e0
         endif
         nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
         
         start(1) = imin + 1 - imin0
         count(1) = imax - imin + 1
         start(2) = jmin + 1 - jmin0
         count(2) = jmax - jmin + 1
         
         start(3)=nrec
         count(3) = 1
         
         nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
         call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)
         
         do j=jmin,jmax
            do i=imin,imax
               var(i,j)=var(i,j)*scale_factor+add_offset
            end do
         end do

         if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
            where (var(:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:)=fillvar
         endif

      endif
   endif
   
   return
   end subroutine ionc4_read_subxyt_real

! ***************************************************************
! * subroutine ionc4_read_subxyt_double               *
! *      *
! * auteur         :    rramel                                   *
! * org            :    ALYOTECH                                 *
! * date creation  :    15/01/11                                *
! * derniere modif :    15/01/11 par rramel                   *
! *      *
! * Role : lit une sous matrice d'une variable double 3D        *
! *                                     lon, lat, temps         *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - nrec        :  numero d'enrtegistrement       *
! *       - imin0        : indice du premier element en i
! *       - jmin0      :  indice du premier element en j
! *       sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_subxyt_double(nom_fichier,var_name,var,&
                                     imin,imax,jmin,jmax,nrec,imin0,jmin0,fillvar)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,imin0,jmin0
   real(kind=8),dimension(imin:imax,jmin:jmax) :: var
   integer :: nrec
   real(kind=8),optional,intent(in) :: fillvar

! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,nc_err_fill
   integer :: dim_time_id,dim_time,i,j
   real(kind=8) :: add_offset,scale_factor,fillvalue_read
   integer,dimension(3) :: start,count


! ******** FONCTIONS **************
   integer :: trim

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_subxyt_double"

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else

      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))

      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))

      if (nrec .gt. dim_time) then
         call ionc4_err(ionc_err_t,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else

         nc_err = nf90_inq_varid(nc_id,var_name,var_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)

         nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
         if (nc_err.ne.nf90_noerr) then
            scale_factor = 1.d0
         endif
         nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
         if (nc_err.ne.nf90_noerr) then
            add_offset = 0.d0
         endif
         nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)

         start(1) = imin + 1 - imin0
         count(1) = imax - imin + 1
         start(2) = jmin + 1 - jmin0
         count(2) = jmax - jmin + 1

         start(3)=nrec
         count(3) = 1

         nc_err = nf90_get_var(nc_id,var_id,var,start=start,count=count)
         call ionc4_err(nc_err,ionc_rout,'nf90_get_var',var_name)

         do j=jmin,jmax
            do i=imin,imax
               var(i,j)=var(i,j)*scale_factor+add_offset
            end do
         end do

         if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
            where (var(:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:)=fillvar
         endif

      endif
   endif

   return
   end subroutine ionc4_read_subxyt_double
   
! ***************************************************************
! * subroutine ionc4_read_sublxyt_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    12/12/03                                *
! * derniere modif :    12/12/03 par jfleroux                   *
! *      *
! * Role : lit une sous matrice d'une variable reelle 3D        *
! * lon, lat, temps dans une variable(ilmin:ilmin+(imax-imin),  *
! *            jlmim:jlmin+(jmax-jmin))  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - nrec        :  numero d'enrtegistrement       *
! *       - ilmin       :  borne imin variable de lecture *
! *       - jlmin       :  borne jmin variable de lecture *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_sublxyt_real(nom_fichier,var_name,var, &
                                      imin,imax,jmin,jmax,nrec,ilmin,jlmin,fillvar)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,ilmin,jlmin
   real(kind=4) :: var(ilmin:ilmin+(imax-imin),jlmin:jlmin+(jmax-jmin))
   integer :: nrec
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,dim_time_id,dim_time,i,j,nc_err_fill
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
   real(kind=4) :: bid(imin:imax,jmin:jmax)
   integer,dimension(3):: start,count
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_sublxyt_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      
      if (nrec .gt. dim_time) then
         call ionc4_err(ionc_err_t,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         
         nc_err = nf90_inq_varid(nc_id,var_name,var_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
         nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
         if (nc_err.ne.nf90_noerr) then
            scale_factor = 1.e0
         endif
         nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
         if (nc_err.ne.nf90_noerr) then
            add_offset = 0.e0
         endif
         nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
         
         start(1) = imin+1
         count(1) = imax - imin + 1
         start(2) = jmin+1
         count(2) = jmax - jmin + 1
         
         start(3)=nrec
         count(3) = 1
         
         nc_err = nf90_get_var(nc_id,var_id,bid,start=start,count=count)
         call ionc4_err(nc_err,ionc_rout,'nf90_get_var', trim(var_name))
         
         do j=jmin,jmax
            do i=imin,imax
               var(ilmin+(i-imin),jlmin+(j-jmin)) = &
               bid(i,j)*scale_factor+add_offset
            end do
         end do

         if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
            where (var(:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:)=fillvar
         endif

      endif
   endif
   
   return
   end subroutine ionc4_read_sublxyt_real
   
! ***************************************************************
! * subroutine ionc4_write_zxy_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 3D = lon, lat, z           *
! *        a partir d'une variable reelle 3D = z , lon, lat     *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - kmax        :  borne max dim k*
! ***************************************************************
   subroutine ionc4_write_zxy_real(nom_fichier,var_name,var,imin,imax,jmin,jmax,kmax)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,kmax
   real(kind=4),dimension(kmax,imin:imax,jmin:jmax) :: var
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   real(kind=4) :: add_offset,scale_factor
   integer :: dimi,dimj,dimk
   integer :: i,j,k,indi,indj,indk
   integer,dimension(3) :: start,count
   real(kind=4),dimension(:,:,:),allocatable :: bid
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_zxy_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      call ionc4_readbounds(nc_id)
      
      dimi = (ionc_lima - ionc_limi)/ ionc_pasi + 1
      dimj = (ionc_ljma - ionc_ljmi)/ ionc_pasj + 1
      dimk = (ionc_lkma - ionc_lkmi)/ ionc_pask + 1
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      allocate(bid(dimi,dimj,dimk))
      start = (/1,1,1/)
      count = (/dimi,dimj,dimk/)
      
      if (scale_factor /= 1.e0 .and. add_offset /=0.e0) then
         indk = 1
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            indi = 1
            do i=ionc_limi,ionc_lima, ionc_pasi
               indj = 1
               do j=ionc_ljmi,ionc_ljma, ionc_pasj
                  bid(indi,indj,indk)=real(nint((var(i,j,k)- add_offset)/scale_factor))
                  indj = indj + 1
               enddo
               indi = indi + 1
            enddo
            indk = indk + 1
         enddo
      else
         indk = 1
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            indi = 1
            do i=ionc_limi,ionc_lima, ionc_pasi
               indj = 1
               do j=ionc_ljmi,ionc_ljma, ionc_pasj
                  bid(indi,indj,indk)=var(i,j,k)
                  indj = indj + 1
               enddo
               indi = indi + 1
            enddo
            indk = indk + 1
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid) 
   endif
   
   return
   end subroutine ionc4_write_zxy_real
   
   
! ***************************************************************
! * subroutine ionc4_read_zxy_real                  *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : lit une variable reelle 3D = lon, lat, z             *
! *        dans une variable reelle 3D = z , lon, lat           *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - kmax        :  borne max dim k*
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_zxy_real(nom_fichier,var_name,var,imin,imax,jmin,jmax,kmax,fillvar)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,kmax
   real(kind=4),dimension(kmax,imin:imax,jmin:jmax) :: var
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id,nc_err,nc_err_fill
   integer :: i,j,k,indi, indj
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
   integer,dimension(3) :: start,count
   real(kind=4),dimension(imax-imin+1,jmax-jmin+1,kmax) :: bid
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_zxy_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
      
      start(1) = 1
      count(1) = imax - imin + 1
      start(2) = 1
      count(2) = jmax - jmin + 1
      start(3)= 1
      count(3) = kmax
      
      nc_err = nf90_get_var(nc_id,var_id,bid,start=start,count=count)
      
      do k = 1,kmax
         indj = 1
         do j=jmin,jmax
            indi = 1
            do i=imin,imax
               var(i,j,k)=bid(indi,indj,k)*scale_factor+add_offset
               indi = indi + 1
            end do
            indj = indj + 1
         end do
      end do

      if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
         where (var(:,:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:,:)=fillvar
      endif

   endif
   
   return
   end subroutine ionc4_read_zxy_real
   
! ***************************************************************
! * subroutine ionc4_write_zxy_double                  *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 3D = lon, lat, z           *
! *        a partir d'une variable reelle 3D = z , lon, lat     *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - kmax        :  borne max dim k*
! ***************************************************************
   subroutine ionc4_write_zxy_double(nom_fichier,var_name,var,imin,imax,jmin,jmax,kmax)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,kmax
   real(kind=8),dimension(kmax,imin:imax,jmin:jmax) :: var
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   real(kind=4) :: add_offset,scale_factor
   integer :: dimi,dimj,dimk
   integer :: i,j,k,indi,indj,indk
   integer,dimension(3) :: start,count
   real(kind=8),dimension(:,:,:),allocatable :: bid
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_zxy_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      call ionc4_readbounds(nc_id)
      
      dimi = (ionc_lima - ionc_limi)/ ionc_pasi + 1
      dimj = (ionc_ljma - ionc_ljmi)/ ionc_pasj + 1
      dimk = (ionc_lkma - ionc_lkmi)/ ionc_pask + 1
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      allocate(bid(dimi,dimj,dimk))
      start = (/1,1,1/)
      count = (/dimi,dimj,dimk/)
      
      if (scale_factor /= 1.e0 .and. add_offset /=0.e0) then
         indk = 1
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            indi = 1
            do i=ionc_limi,ionc_lima, ionc_pasi
               indj = 1
               do j=ionc_ljmi,ionc_ljma, ionc_pasj
                  bid(indi,indj,indk)=REAL(nint((var(i,j,k)- add_offset)/scale_factor),8)
                  indj = indj + 1
               enddo
               indi = indi + 1
            enddo
            indk = indk + 1
         enddo
      else
         indk = 1
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            indi = 1
            do i=ionc_limi,ionc_lima, ionc_pasi
               indj = 1
               do j=ionc_ljmi,ionc_ljma, ionc_pasj
                  bid(indi,indj,indk)=var(i,j,k)
                  indj = indj + 1
               enddo
               indi = indi + 1
            enddo
            indk = indk + 1
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid) 
   endif
   
   return
   end subroutine ionc4_write_zxy_double
   
! ***************************************************************
! * subroutine ionc4_read_zxy_double                 *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : lit une variable reelle 3D = lon, lat, z             *
! *        dans une variable reelle 3D = z , lon, lat           *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - kmax        :  borne max dim k*
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_zxy_double(nom_fichier,var_name,var,imin,imax,jmin,jmax,kmax,fillvar)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,kmax
   real(kind=8),dimension(kmax,imin:imax,jmin:jmax) :: var
   real(kind=8),optional,intent(in) :: fillvar
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,var_id,nc_err,nc_err_fill
   integer :: i,j,k,indi, indj
   real(kind=8) :: add_offset,scale_factor,fillvalue_read
   integer,dimension(3) :: start,count
   real(kind=8),dimension(imax-imin+1,jmax-jmin+1,kmax) :: bid
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_zxy_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.d0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.d0
      endif
      nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
      
      start(1) = 1
      count(1) = imax - imin + 1
      start(2) = 1
      count(2) = jmax - jmin + 1
      start(3)= 1
      count(3) = kmax
      
      nc_err = nf90_get_var(nc_id,var_id,bid,start=start,count=count)
      
      do k = 1,kmax
         indj = 1
         do j=jmin,jmax
            indi = 1
            do i=imin,imax
               var(i,j,k)=bid(indi,indj,k)*scale_factor+add_offset
               indi = indi + 1
            end do
            indj = indj + 1
         end do
      end do

      if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
         where (var(:,:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:,:)=fillvar
      endif

   endif
   
   return
   end subroutine ionc4_read_zxy_double

! ***************************************************************
! * subroutine ionc4_write_xyz_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 3D = lon, lat, z           *
! *        a partir d'une variable reelle 3D = lon, lat,z       * 
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - kmax        :  borne max dim j*
! ***************************************************************
      subroutine ionc4_write_xyz_real(nom_fichier,var_name,var,imin,imax,jmin,jmax,kmax)

! ******** PARAMETRES DE LA SUBROUTINE *******
      character(len=*) :: nom_fichier,var_name
      integer :: imin,imax,jmin,jmax,kmax
      real(kind=4),dimension(imin:imax,jmin:jmax,1:kmax) :: var

! ******** VARIABLES DE TRAVAIL **************
      integer :: nc_id,var_id,nc_err
      real(kind=4) :: add_offset,scale_factor
      integer  :: dimi,dimj,dimk
      
! ******** FIN DES DECLARATIONS **************
      ionc_rout = "ionc4_write_xyz_real"

      call ionc4_corres(nom_fichier,nc_id)
      if (nc_id .eq. 0) then
          call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
          stop
      else

         call ionc4_readbounds(nc_id)         
        
         dimi = (ionc_lima - ionc_limi)/ ionc_pasi + 1
         dimj = (ionc_ljma - ionc_ljmi)/ ionc_pasj + 1
         dimk = (ionc_lkma - ionc_lkmi)/ ionc_pask + 1

         nc_err = nf90_inq_varid(nc_id,var_name,var_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
         nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
         if (nc_err.ne.nf90_noerr) then
             scale_factor = 1.e0
         endif
         nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
         if (nc_err.ne.nf90_noerr) then
             add_offset = 0.e0
         endif

         call ionc4_putr4_xyz_real(nc_id,var_id,var,imin,imax,jmin,jmax, &
                                      kmax,dimi-1,dimj-1,dimk-1, &
                                      scale_factor,add_offset)
     
      endif

      return
      end subroutine ionc4_write_xyz_real

! ***************************************************************
! * subroutine ionc4_put_xyz_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    05/03/01                                *
! * derniere modif :    05/03/01 par jfleroux                   *
! *      *
! * Role :                                                      * 
! *      *
! * Parametres :      *
! ***************************************************************
      subroutine ionc4_putr4_xyz_real(nc_id,var_id,var, &
                                             imin,imax,jmin,jmax,kmax, &
                                             dimi,dimj,dimk,scale_factor,add_offset)

! ******** PARAMETRES DE LA SUBROUTINE *******

      integer :: nc_id,var_id
      integer :: imin,imax,jmin,jmax,kmax
      integer :: dimi,dimj,dimk
      real(kind=4) :: var(imin:imax,jmin:jmax,1:kmax)
      real(kind=4) :: add_offset,scale_factor

! ******** VARIABLES DE TRAVAIL **************
      integer :: nc_err
      integer,dimension(3)::  start,count
      integer :: i,j,k,indi,indj,indk
      real(kind=4) :: bid(0:dimi,0:dimj,0:dimk)
      
! ******** FIN DES DECLARATIONS **************
      ionc_rout = "ionc4_putr4_xyz_real"

      start(1) = 1
      start(2) = 1
      start(3) = 1

      count(1) = dimi+1
      count(2) = dimj+1
      count(3) = dimk+1
      
      if (scale_factor /= 1.e0 .and. add_offset /= 0.e0) then
         indk = 0
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            indi = 0
            do i=ionc_limi,ionc_lima, ionc_pasi
               indj = 0
               do j=ionc_ljmi,ionc_ljma, ionc_pasj
                  bid(indi,indj,indk)=real(nint((var(i,j,k)- add_offset)/scale_factor))
                  indj = indj + 1
               enddo
               indi = indi + 1
            enddo
            indk = indk + 1
         enddo
      else
         indk = 0
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            indi = 0
            do i=ionc_limi,ionc_lima, ionc_pasi
               indj = 0
               do j=ionc_ljmi,ionc_ljma, ionc_pasj
                  bid(indi,indj,indk)=var(i,j,k)
                  indj = indj + 1
               enddo
               indi = indi + 1
            enddo
            indk = indk + 1
         enddo
      endif
 
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',' ')

      return
      end subroutine ionc4_putr4_xyz_real

! ***************************************************************
! * subroutine ionc4_write_sublxyt_real                   *
! *      *
! * auteur         :    jfleroux                              *
! * org            :    IFREMER                                 *
! * date creation  :    16/12/03                                *
! * derniere modif :    01/02/2011 par rramel (ALYOTECH)                 *
! *      *
! * Role : ecrit une variable reelle 3D = lon, lat, temps    *
! *        a partir d'une variable reelle 2D =  lon, lat     *
! *        (ilmin:ilmin+(imax-imin),                            *
! *         jlmim:jlmin+(jmax-jmin))                            *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! *       - ilmin       :  borne min var d'entree        *
! *       - jlmin       :  borne min var d'entree        *
! ***************************************************************
   subroutine ionc4_write_sublxyt_real(nom_fichier,var_name,var, &
                                       imin,imax,jmin,jmax,nrec, &
                                       ilmin,jlmin)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,jmin,jmax
   integer :: ilmin,jlmin
   real(kind=4) :: var( ilmin:ilmin+(imax-imin),  &
                        jlmin:jlmin+(jmax-jmin))
   integer :: nrec

! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,dim_time_id,dim_time
   integer :: i,j,dimi,dimj,indi,indj
   real(kind=4) :: add_offset,scale_factor
   real(kind=4) :: varl(imin:imax,jmin:jmax)
   integer,dimension(3) :: start,count
   real(kind=4),dimension(:,:),allocatable :: bid

! ******** FONCTIONS **************
   integer :: trim

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_sublxyt_real"

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else

      call ionc4_readbounds(nc_id)

      dimi = (ionc_lima - ionc_limi)/ ionc_pasi + 1
      dimj = (ionc_ljma - ionc_ljmi)/ ionc_pasj + 1

      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))

      if (nrec .eq. 0) then
         nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
         call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         dim_time = nrec
      endif

      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)

      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif

      do j=jmin,jmax
         do i=imin,imax
            varl(i,j)=var(ilmin+(i-imin),jlmin+(j-jmin))
         end do
      end do

      allocate(bid(dimi,dimj))
      start = (/1,1,dim_time/)
      count = (/dimi,dimj,1/)

      if (scale_factor /=1.e0 .and. add_offset /= 0.e0) then
         indi = 1
         do i=ionc_limi,ionc_lima, ionc_pasi
            indj = 1
            do j=ionc_ljmi,ionc_ljma, ionc_pasj
               bid(indi,indj)=real(nint((varl(i,j)- add_offset)/scale_factor))
               indj = indj + 1
            enddo
            indi = indi + 1
         enddo
      else
         indi = 1
         do i=ionc_limi,ionc_lima, ionc_pasi
            indj = 1
            do j=ionc_ljmi,ionc_ljma, ionc_pasj
               bid(indi,indj)=varl(i,j)
               indj = indj + 1
            enddo
            indi = indi + 1
         enddo
      endif

      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid)

   endif

   return

   end subroutine ionc4_write_sublxyt_real

! ***************************************************************
! * subroutine ionc4_write_sublzxyt_real                   *
! *      *
! * auteur         :    jfleroux                              *
! * org            :    IFREMER                                 *
! * date creation  :    16/12/03                                *
! * derniere modif :    16/12/03 par jfleroux                   *
! *      *
! * Role : ecrit une variable reelle 4D = lon, lat, z, temps    *
! *        a partir d'une variable reelle 3D = z , lon, lat     *
! *        (klmin:klmin+(kmax-kmin),                            *
! *         ilmin:ilmin+(imax-imin),                            *
! *         jlmim:jlmin+(jmax-jmin))                            *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - kmin        :  borne min dim k        *
! *       - kmax        :  borne max dim k*
! *       - nrec        :  numero d'enregistrement*
! *                        0 si t est illimite*
! *       - ilmin       :  borne min var d'entree        *
! *       - jlmin       :  borne min var d'entree        *
! *       - klmin       :  borne min var d'entree        *
! ***************************************************************
   subroutine ionc4_write_sublzxyt_real(nom_fichier,var_name,var, &
                                        imin,imax,jmin,jmax,kmin,kmax,nrec, &
                                        ilmin,jlmin,klmin)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,kmin,kmax
   integer :: ilmin,jlmin,klmin
   real(kind=4) :: var(klmin:klmin+(kmax-kmin),    &
                       ilmin:ilmin+(imax-imin),  &
                       jlmin:jlmin+(jmax-jmin))
   integer :: nrec
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,dim_time_id,dim_time
   integer :: i,j,k,dimi,dimj,dimk,indi,indj,indk
   real(kind=4) :: add_offset,scale_factor
   real(kind=4) :: varl(kmax,imin:imax,jmin:jmax)
   integer,dimension(4) :: start,count
   real(kind=4),dimension(:,:,:),allocatable :: bid
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_sublzxyt_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      call ionc4_readbounds(nc_id)
      
      dimi = (ionc_lima - ionc_limi)/ ionc_pasi + 1
      dimj = (ionc_ljma - ionc_ljmi)/ ionc_pasj + 1
      dimk = (ionc_lkma - ionc_lkmi)/ ionc_pask + 1
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      if (nrec .eq. 0) then
         nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
         call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         dim_time = nrec
      endif
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      
      do k=0,kmax-1
         do j=jmin,jmax
            do i=imin,imax
               varl(i,j,k)=var(klmin+(k-kmin),ilmin+(i-imin),jlmin+(j-jmin))
            end do
         end do
      end do
      
      allocate(bid(dimi,dimj,dimk))
      start = (/1,1,1,dim_time/)
      count = (/dimi,dimj,dimk,1/)
      
      if (scale_factor /=1.e0 .and. add_offset /= 0.e0) then
         indk = 1
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            indi = 1
            do i=ionc_limi,ionc_lima, ionc_pasi
               indj = 1
               do j=ionc_ljmi,ionc_ljma, ionc_pasj
                  bid(indi,indj,indk)=real(nint((varl(i,j,k)- add_offset)/scale_factor))
                  indj = indj + 1
               enddo
               indi = indi + 1
            enddo
            indk = indk + 1
         enddo
      else
         indk = 1
         do k=ionc_lkmi,ionc_lkma, ionc_pask
            indi = 1
            do i=ionc_limi,ionc_lima, ionc_pasi
               indj = 1
               do j=ionc_ljmi,ionc_ljma, ionc_pasj
                  bid(indi,indj,indk)=varl(i,j,k)
                  indj = indj + 1
               enddo
               indi = indi + 1
            enddo
            indk = indk + 1
         enddo
      endif
      
      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid)

   endif

   return

   end subroutine ionc4_write_sublzxyt_real
   
! ***************************************************************
! * subroutine ionc4_read_zxyt_real                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : lit une variable reelle 4D = lon, lat, z, temps      *
! *        dans une variable reelle 3D = z , lon, lat           *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - kmax        :  borne max dim k*
! *       - nrec        :  numero d'enregistrement*
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_zxyt_real(nom_fichier,var_name,var,imin,imax,jmin,jmax,kmax,nrec,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,kmax,nrec
   real(kind=4),dimension(kmax,imin:imax,jmin:jmax) :: var
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,dim_time_id,dim_time
   integer :: var_id,nc_err,nc_err_fill
   integer,dimension(4) :: start,count
   integer :: i,j,k,indi, indj
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
   real(kind=4),dimension(imax-imin+1,jmax-jmin+1,kmax) :: bid
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_zxyt_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      
      if (nrec .gt. dim_time) then
         call ionc4_err(ionc_err_t,ionc_rout,' ',trim(ionc_nomtime))
      else
         
         nc_err = nf90_inq_varid(nc_id,var_name,var_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
         nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
         if (nc_err.ne.nf90_noerr) then
            scale_factor = 1.e0
         endif
         nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
         if (nc_err.ne.nf90_noerr) then
            add_offset = 0.e0
         endif
         nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
         
         start(1) = 1
         count(1) = imax - imin + 1
         start(2) = 1
         count(2) = jmax - jmin + 1
         start(3)= 1
         count(3) = kmax
         start(4)= nrec
         count(4) = 1
         
         nc_err = nf90_get_var(nc_id,var_id,bid,start=start,count=count)
         
         do k = 1,kmax
            indj = 1
            do j=jmin,jmax
               indi = 1
               do i=imin,imax
                  var(i,j,k)=bid(indi,indj,k)*scale_factor+add_offset
                  indi = indi + 1
               end do
               indj = indj + 1
            end do
         end do

         if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
            where (var(:,:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:,:)=fillvar
         endif

      endif
   endif
   
   return
   end subroutine ionc4_read_zxyt_real
   
! ***************************************************************
! * subroutine ionc4_read_subzxyt_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    09/12/03                                *
! * derniere modif :    09/12/03 par jfleroux                   *
! *      *
! * Role : lit une sous matrice (i,j,k) d'une variable reelle 4D*
! *                                     lon, lat, z, temps      *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - kmin        :  borne min dim k        *
! *       - kmax        :  borne max dim k*
! *       - nrec        :  numero d'enrtegistrement       *
! *       - imin0        : indice du premier element en i
! *       - jmin0      :  indice du premier element en j
! *       - kmin0      :  indice du premier element en k
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_subzxyt_real(nom_fichier,var_name,var, &
                                      imin,imax,jmin,jmax,kmin,kmax,nrec,&
                                      imin0,jmin0,kmin0,fillvar)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier
   character(len=*) :: var_name
   integer :: imin,imax,jmin,jmax, kmin, kmax,imin0,jmin0,kmin0
   real(kind=4),dimension(imin:imax,jmin:jmax,kmin:kmax),intent(inout) :: var
   integer :: nrec
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id
   integer :: var_id
   integer :: nc_err,nc_err_fill
   integer :: dim_time_id
   integer :: dim_time
   integer :: i,j,k
   integer :: indi, indj, indk
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
   real(kind=4),dimension(imax-imin+1,jmax-jmin+1,kmax-kmin+1) :: bid
   integer,dimension(4) :: start,count
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_subzxyt_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      
      if (nrec .gt. dim_time) then
         call ionc4_err(ionc_err_t,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         
         nc_err = nf90_inq_varid(nc_id,var_name,var_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
         nc_err = nf90_get_att(nc_id,var_id,'scale_factor', scale_factor)
         if (nc_err.ne.nf90_noerr) then
            scale_factor = 1.e0
         endif
         nc_err = nf90_get_att(nc_id,var_id,'add_offset', add_offset)
         if (nc_err.ne.nf90_noerr) then
            add_offset = 0.e0
         endif
         nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
         
         start = (/imin+1-imin0,jmin+1-jmin0,kmin+1-kmin0,nrec/)
         count = (/imax-imin+1,jmax-jmin+1,kmax-kmin+ 1,1/)
         
         nc_err = nf90_get_var(nc_id,var_id,bid,start=start,count=count)
         
         call ionc4_err(nc_err,ionc_rout,'nf90_get_var', var_name)
         
         indk = 1
         do k = kmin,kmax
            indj = 1
            do j=jmin,jmax
               indi = 1
               do i=imin,imax
                  var(i,j,k)=bid(indi,indj,indk)*scale_factor+add_offset
                  indi = indi + 1
               end do
               indj = indj + 1
            end do
            indk = indk + 1
         end do

         if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
            where (var(:,:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:,:)=fillvar
         endif
         
      endif
   endif
   
   return
   end subroutine ionc4_read_subzxyt_real

! ***************************************************************
! * subroutine ionc4_read_subzxyt_double                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    09/12/03                                *
! * derniere modif :    09/12/03 par jfleroux                   *
! *      *
! * Role : lit une sous matrice (i,j,k) d'une variable reelle 4D*
! *                                     lon, lat, z, temps      *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - kmin        :  borne min dim k        *
! *       - kmax        :  borne max dim k*
! *       - nrec        :  numero d'enrtegistrement       *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_subzxyt_double(nom_fichier,var_name,var, &
                                        imin,imax,jmin,jmax,kmin,kmax,nrec,&
                                        imin0,jmin0,kmin0,fillvar)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier
   character(len=*) :: var_name
   integer :: imin,imax,jmin,jmax, kmin, kmax,imin0,jmin0,kmin0
   real(kind=8),dimension(imin:imax,jmin:jmax,kmin:kmax),intent(inout) :: var
   integer :: nrec
   real(kind=8),optional,intent(in) :: fillvar
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id
   integer :: var_id
   integer :: nc_err,nc_err_fill
   integer :: dim_time_id
   integer :: dim_time
   integer :: i,j,k
   integer :: indi, indj, indk
   real(kind=8) :: add_offset,scale_factor,fillvalue_read
   real(kind=8),dimension(imax-imin+1,jmax-jmin+1,kmax-kmin+1) :: bid
   integer,dimension(4) :: start,count
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_subzxyt_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      
      if (nrec .gt. dim_time) then
         call ionc4_err(ionc_err_t,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         
         nc_err = nf90_inq_varid(nc_id,var_name,var_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
         nc_err = nf90_get_att(nc_id,var_id,'scale_factor', scale_factor)
         if (nc_err.ne.nf90_noerr) then
            scale_factor = 1.d0
         endif
         nc_err = nf90_get_att(nc_id,var_id,'add_offset', add_offset)
         if (nc_err.ne.nf90_noerr) then
            add_offset = 0.d0
         endif
         nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
         
         start = (/imin+1-imin0,jmin+1-jmin0,kmin+1-kmin0,nrec/)
         count = (/imax-imin+1,jmax-jmin+1,kmax-kmin+1,1/)
         
         nc_err = nf90_get_var(nc_id,var_id,bid,start=start,count=count)
         !print*, 'limites    : ', imin,imax,jmin,jmax
         !PRINT*, 'Dim var 3D : ', SIZE(bid,dim=1), SIZE(bid,dim=2), SIZE(bid,dim=3)
         !print*, 'start      : ', start
         !print*, 'count      : ', count
         !print*, 'somme      : ', start + count
         !print*, ' ' 
         call ionc4_err(nc_err,ionc_rout,'nf90_get_var', var_name)
         
         !!! Lis correctement les valeurs d'entree !!! 
         indk = 1
         do k = 1,kmax
            indj = 1
            do j=jmin,jmax
               indi = 1
               do i=imin,imax
                  var(i,j,k)=bid(indi,indj,indk)*scale_factor+add_offset
                  indi = indi + 1
               end do
               indj = indj + 1
            end do
            indk = indk + 1
         end do


         if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
            where (var(:,:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:,:)=fillvar
         endif

      endif
   endif
   
   return
   end subroutine ionc4_read_subzxyt_double
   
! ***************************************************************
! * subroutine ionc4_read_sublzxyt_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    16/12/03                                *
! * derniere modif :    01/06/06 par jfleroux                   *
! *      *
! * Role : lit une sous matrice d'une variable reelle 3D        *
! * lon, lat, temps dans une variable(klmin:klmin+(kmax-kmin),  *
! *                                   ilmin:ilmin+(imax-imin),  *
! *            jlmim:jlmin+(jmax-jmin))  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - kmin        :  borne min dim k        *
! *       - kmax        :  borne max dim k*
! *       - nrec        :  numero d'enrtegistrement       *
! *       - ilmin       :  borne imin variable de lecture *
! *       - jlmin       :  borne jmin variable de lecture *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_sublzxyt_real(nom_fichier,var_name,var, &
                                       imin,imax,jmin,jmax,kmin,kmax,nrec, &
                                       ilmin,jlmin,klmin,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) ::   nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,kmin,kmax
   integer :: ilmin,jlmin,klmin
   real(kind=4) :: var(klmin:klmin+(kmax-kmin), &
                       ilmin:ilmin+(imax-imin),    &
                       jlmin:jlmin+(jmax-jmin))
   integer :: nrec
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,dim_time_id,dim_time,nc_err_fill
   integer :: i,j,k
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
   real(kind=4),dimension(imin:imax,jmin:jmax,kmin:kmax) :: bid
   integer,dimension(4) :: start,count
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_sublzxyt_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      
      if (nrec .gt. dim_time) then
         call ionc4_err(ionc_err_t,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      else
         
         nc_err = nf90_inq_varid(nc_id,var_name,var_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
         nc_err = nf90_get_att(nc_id,var_id,'scale_factor', scale_factor)
         if (nc_err.ne.nf90_noerr) then
            scale_factor = 1.e0
         endif
         nc_err = nf90_get_att(nc_id,var_id,'add_offset', add_offset)
         if (nc_err.ne.nf90_noerr) then
            add_offset = 0.e0
         endif
         nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
         
         start(1) = imin+1
         count(1) = imax - imin + 1
         start(2) = jmin+1
         count(2) = jmax - jmin + 1
         start(3)=  kmin
         count(3) = kmax - kmin + 1
         
         start(4)=nrec
         count(4) = 1
         
         nc_err = nf90_get_var(nc_id,var_id,bid,start=start,count=count)
         call ionc4_err(nc_err,ionc_rout,'nf90_get_var', trim(var_name))
         
         do k=kmin,kmax
            do j=jmin,jmax
               do i=imin,imax
                  var(klmin+(k-kmin),ilmin+(i-imin),jlmin+(j-jmin))=bid(i,j,k)*scale_factor+add_offset
               end do
            end do
         end do

         if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
            where (var(:,:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:,:)=fillvar
         endif

      endif
   endif
   
   return
   end subroutine ionc4_read_sublzxyt_real
   
! ***************************************************************
! * subroutine ionc4_read_zxyt_double                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : lit une variable double 4D = lon, lat, z, temps      *
! *        dans une variable double 3D = z , lon, lat           *
! *         l'indice du temps est fixe par nrec                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - jmin        :  borne min dim j        *
! *       - jmax        :  borne max dim j*
! *       - kmax        :  borne max dim k*
! *       - nrec        :  numero d'enregistrement*
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_zxyt_double(nom_fichier,var_name,var, &
                                     imin,imax,jmin,jmax,kmax,nrec,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,jmin,jmax,kmax,nrec
   real(kind=8),dimension(kmax,imin:imax,jmin:jmax) :: var
   real(kind=8),optional,intent(in) :: fillvar
   
! ******** VARIABLES DE TRAVAIL **************
   
   integer :: nc_id,dim_time_id,dim_time,var_id,nc_err,nc_err_fill
   integer,dimension(4) :: start,count
   integer :: i,j,k,indi, indj
   real(kind=8) :: add_offset,scale_factor,fillvalue_read
   real(kind=8),dimension(imax-imin+1,jmax-jmin+1,kmax) ::  bid
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_zxyt_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=dim_time)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      
      if (nrec .gt. dim_time) then
         call ionc4_err(ionc_err_t,ionc_rout,' ',trim(ionc_nomtime))
      else
         
         nc_err = nf90_inq_varid(nc_id,var_name,var_id)
         call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
         
         nc_err = nf90_get_att(nc_id,var_id,'scale_factor', scale_factor)
         if (nc_err.ne.nf90_noerr) then
            scale_factor = 1.d0
         endif
         nc_err = nf90_get_att(nc_id,var_id,'add_offset', add_offset)
         if (nc_err.ne.nf90_noerr) then
            add_offset = 0.d0
         endif
         nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
         
         start(1) = 1
         count(1) = imax - imin + 1
         start(2) = 1
         count(2) = jmax - jmin + 1
         start(3)= 1
         count(3) = kmax
         start(4)= nrec
         count(4) = 1
         
         nc_err = nf90_get_var(nc_id,var_id,bid,start=start,count=count)
         
         do k = 1,kmax
            indj = 1
            do j=jmin,jmax
               indi = 1
               do i=imin,imax
                  var(i,j,k)=bid(indi,indj,k)*scale_factor+add_offset
                  indi = indi + 1
               end do
               indj = indj + 1
            end do
         end do

         if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
            where (var(:,:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:,:)=fillvar
         endif

      endif
   endif
   
   return
   end subroutine ionc4_read_zxyt_double
   
! ***************************************************************
! * subroutine ionc4_read_zxt_real                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    05/03/02                                *
! * derniere modif :    05/03/02 par rramel                   *
! *      *
! * Role : lit une variable reelle 3D = lon, z, temps           *
! *        dans une variable reelle 2D = z , lon                *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - imin        :  borne min dim i        *
! *       - imax        :  borne max dim i        *
! *       - kmax        :  borne max dim k*
! *       - nrec        :  numero d'enrtegistrement       *
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_zxt_real(nom_fichier,var_name,var, &
                                  imin,imax,kmax,nrec,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier,var_name
   integer :: imin,imax,kmax,nrec
   real(kind=4),dimension(kmax,imin:imax) :: var
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** DECLARATIONS EXTERNES *************
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,nc_err_fill
   integer,dimension(3) :: start,count
   integer :: i,k,indi
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
   real(kind=4),dimension(imax-imin+1,kmax) :: bid
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_zxt_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor', scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset', add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
      
      start(1) = 1
      count(1) = imax - imin + 1
      
      start(2)= 1
      count(2) = kmax
      
      start(3)=nrec
      count(3) = 1
      
      nc_err = nf90_get_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_var', var_name)
      
      do k = 1,kmax
         indi = 1
         do i=imin,imax
            var(k,i)=bid(indi,k)*scale_factor+add_offset
            indi = indi + 1
         end do
      end do

      if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
         where (var(:,:)==fillvalue_read*scale_factor+add_offset ) var(:,:)=fillvar
      endif

   endif
   
   return
   end subroutine ionc4_read_zxt_real
   
! ***************************************************************
! * subroutine ionc4_title                           *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    23/02/01                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                   *
! *      *
! * Role : attribut titre d'un fichier netcdf                   *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - titre       :  titre du fichier        *
! ***************************************************************
   subroutine ionc4_title(nom_fichier,titre)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,titre
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_title"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_put_att(nc_id,nf90_global,'title',titre )
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:title',' ')
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',' ')
   endif
   
   return
   end subroutine ionc4_title
   
! ***************************************************************
! * subroutine ionc4_history                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    27/02/01                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                    *
! *      *
! * Role : attribut history d'un fichier netcdf                 *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - history    : historique de creation du fichier*
! ***************************************************************
   subroutine ionc4_history(nom_fichier,history)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,history
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_history"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_put_att(nc_id,nf90_global,'history',history )
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:history',' ')
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',' ')
   endif
   
   return
   end subroutine ionc4_history
   
! ***************************************************************
! * subroutine ionc4_tunite                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    27/02/01                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                    *
! *      *
! * Role : modifie l'attribut units de la variable time         *
! *        d'un fichier netcdf                                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - tunite      :  unite de temps                 *
! ***************************************************************
   subroutine ionc4_tunite(nom_fichier,tunite)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,tunite
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_time_id,nc_err
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_tunite"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,trim(ionc_nomtime),var_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_put_att(nc_id, var_time_id, 'units',  tunite)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:units')
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',' ')
   endif
   
   return
   end subroutine ionc4_tunite
   
! ***************************************************************
! * subroutine ionc4_torigine                  *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    27/02/01                                *
! * derniere modif :   22/12/10 par rramel (ALYOTECH)                   *
! *      *
! * Role : modifie l'attribut units de la variable time         *
! *        d'un fichier netcdf                                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - torigine    :  origine du temps               *
! ***************************************************************
   subroutine ionc4_torigine(nom_fichier,torigine)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,torigine
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_time_id,nc_err
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_torigine"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,trim(ionc_nomtime),var_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_put_att(nc_id, var_time_id, 'time_origin',  torigine)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att','TIME:time_origin')
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',' ')
   endif
   
   return
   end subroutine ionc4_torigine
   
! ***************************************************************
! * subroutine ionc4_read_torigine                  *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/10/01                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                  *
! *      *
! * Role : lecture de l'attribut time_origin de la varible time *
! *        d'un fichier netcdf                                  *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       sortie:- torigine    :  origine du temps               *
! ***************************************************************
   subroutine ionc4_read_torigine(nom_fichier,torigine)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,torigine
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_time_id,nc_err
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_torigine"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_varid(nc_id,trim(ionc_nomtime),var_time_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_get_att(nc_id, var_time_id, 'time_origin', torigine)
      call ionc4_err(nc_err,ionc_rout,'nf90_get_att','TIME:time_origin')
      
   endif
   
   return
   end subroutine ionc4_read_torigine
   
! ***************************************************************
! * subroutine ionc4_openr                   *
! *      *
! * auteur         :    pgarreau,jfleroux                       *
! * org            :    IFREMER                                 *
! * date creation  :    xx/xx/xx                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)          *
!                       30/07/13 par A. Thevenin (CERFACS)      *
! *                     13/11/13 par L. Debreu (INRIA) : mode_flag for l_out_nc4par *
! *      *
! * Role : ouverture d'un fichier netcdf en lecture      *
! *      *
! * Parametres :      *
! *       entree:- nomfich   :  nom du fichier            *
! ***************************************************************
   subroutine ionc4_openr(nom_fichier,l_in_nc4par)
#ifdef MPI   
  include 'mpif.h'
#endif
 
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   logical,optional,intent(in) :: l_in_nc4par
   
! ******** VARIABLES DE TRAVAIL **************
   logical :: l_useless
   integer :: nc_id,nc_err,mode_flag
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_openr"
   
   call ionc4_corres(nom_fichier,nc_id)
   mode_flag = nf90_nowrite
   if (nc_id .eq. 0) then
#ifdef MPI
      if (present(l_in_nc4par)) then
         if (l_in_nc4par) then
#ifdef key_oasis
            nc_err = nf90_open(nom_fichier,ior(mode_flag,nf90_mpiio),nc_id, &
                               comm = MPI_COMM_MARS, info = MPI_INFO_NULL)
#else
            nc_err = nf90_open(nom_fichier,ior(mode_flag,nf90_mpiio),nc_id, &
                               comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
#endif
         else
            nc_err = nf90_open(nom_fichier,mode_flag,nc_id)
         endif
      else
         nc_err = nf90_open(nom_fichier,mode_flag,nc_id)
      endif
#else
      if (present(l_in_nc4par)) l_useless = l_in_nc4par
      nc_err = nf90_open(nom_fichier,mode_flag,nc_id)
#endif
      call ionc4_err(nc_err,ionc_rout,'nf90_open',nom_fichier)
      ionc_nfich = ionc_nfich + 1
      ionc_nomfich(ionc_nfich) = nom_fichier
      ionc_idfich(ionc_nfich) = nc_id
   endif
   
   return
   end subroutine ionc4_openr
   
   
! ***************************************************************
! * subroutine ionc4_open                   *
! *      *
! * auteur         :    jfleroux                             *
! * org            :    IFREMER                                 *
! * date creation  :    09/01/02                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)          *
!                       30/07/13 par A. Thevenin (CERFACS)      *
! *      *
! * Role : ouverture d'un fichier netcdf en lecture/ecriture      *
! *      *
! * Parametres :      *
! *       entree:- nomfich   :  nom du fichier            *
! ***************************************************************
   subroutine ionc4_open(nom_fichier,l_in_nc4par)
#ifdef MPI   
  include 'mpif.h'
#endif
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   logical,optional,intent(in) :: l_in_nc4par
   
! ******** VARIABLES DE TRAVAIL **************
   logical :: l_useless
   integer :: nc_id,nc_err,mode_flag
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_open"
   
   call ionc4_corres(nom_fichier,nc_id)
   mode_flag = nf90_write
   if (nc_id .eq. 0) then
#ifdef MPI
      if (present(l_in_nc4par)) then
         if (l_in_nc4par) then
#ifdef key_oasis
            nc_err = nf90_open(nom_fichier,ior(mode_flag,nf90_mpiio),nc_id, &
                               comm = MPI_COMM_MARS, info = MPI_INFO_NULL)
#else
            nc_err = nf90_open(nom_fichier,ior(mode_flag,nf90_mpiio),nc_id, &
                               comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
#endif
         else
            nc_err = nf90_open(nom_fichier,mode_flag,nc_id)
         endif
      else
         nc_err = nf90_open(nom_fichier,mode_flag,nc_id)
      endif
#else
      if (present(l_in_nc4par)) l_useless = l_in_nc4par  ! useless for portability only
      nc_err = nf90_open(nom_fichier,mode_flag,nc_id)
#endif
      call ionc4_err(nc_err,ionc_rout,'nf90_open',nom_fichier)
      ionc_nfich = ionc_nfich + 1
      ionc_nomfich(ionc_nfich) = nom_fichier
      ionc_idfich(ionc_nfich) = nc_id
   endif
   
   return
   end subroutine ionc4_open
   
! ***************************************************************
! * subroutine ionc4_close                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    21/02/01                                *
! * derniere modif :   22/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : fermeture d'un fichier netcdf en lecture      *
! *      *
! * Parametres :      *
! *       entree:- nomfich   :  nom du fichier            *
! ***************************************************************
   subroutine ionc4_close(nom_fichier)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,nc_err,indfich
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_close"
   
   nc_id = 0
   
   if (ionc_nfich .GT. 0) then
      do indfich = 1,ionc_nfich
         if (ionc_nomfich(indfich) .EQ.  nom_fichier) then
            nc_id = ionc_idfich(indfich)
            ionc_nomfich(indfich) = ""
            ionc_idfich(indfich)=0
            exit
         endif
      enddo
   endif
   
   if (nc_id .eq. 0) then
   else
      nc_err = nf90_close(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_close',nom_fichier)
   endif
   
   return
   end subroutine ionc4_close

! ***************************************************************
! * subroutine ionc4_sync                   *
! *      *
! * auteur         :    rramel                                *
! * org            :    ALYOTECH                                 *
! * date creation  :    12/01/11                                *
! *      *
! * Role : Synchronisation d'un fichier (flush)      *
! *      *
! * Parametres :      *
! *       entree:- nomfich   :  nom du fichier            *
! ***************************************************************
   subroutine ionc4_sync(nom_fichier)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier

! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,nc_err

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_sync"
   call ionc4_corres(nom_fichier,nc_id)
   nc_err = nf90_sync(nc_id)
   call ionc4_err(nc_err,ionc_rout,'nf90_sync',nom_fichier)

   return
   end subroutine ionc4_sync
   
! ***************************************************************
! * subroutine ionc4_closeall                   *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    21/02/01                                *
! * derniere modif :   22/12/10 par rramel (ALYOTECH)                *
! *      *
! * Role : fermeture de tout les fichiers netcdf               *
! *      *
! * Parametres :      *
! ***************************************************************
   subroutine ionc4_closeall()
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: ifich,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_closeall"
   
   do ifich = 1,ionc_nfich
      nc_err = nf90_close(ionc_idfich(ifich))
!      call ionc4_err(nc_err,ionc_rout,'nf90_close',ionc_nomfich(ifich))
   enddo
   
   return
   end subroutine ionc4_closeall
   
! ***************************************************************
! * subroutine ionc4_read_dimi                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/10/01                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                   *
! *      *
! * Role : lit la taille de la grille en i                      *
! *      *
! * Parametres :      *
! *       entree:- nomfich   :  nom du fichier            *
! *       sortie:- dimi  :  taille de la grille en i        *
! ***************************************************************
   subroutine ionc4_read_dimi(nom_fichier,dimi)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   integer :: dimi
   
! ******** VARIABLES DE TRAVAIL **************
   
! ******** FONCTIONS **************
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_dimi"

   dimi=ionc4_fnlimax(nom_fichier) 

   return
   end subroutine ionc4_read_dimi
   
! ***************************************************************
! * subroutine ionc4_read_dimj                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/10/01                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                   *
! *      *
! * Role : lit la taille de la grille en j                      *
! *      *
! * Parametres :      *
! *       entree:- nomfich   :  nom du fichier            *
! *       sortie:- dimj  :  taille de la grille en j        *
! ***************************************************************
   subroutine ionc4_read_dimj(nom_fichier,dimj)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   integer :: dimj
   
! ******** VARIABLES DE TRAVAIL **************
   
! ******** FONCTIONS **************
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_dimj"
   
   dimj=ionc4_fnljmax(nom_fichier)

   return
   end subroutine ionc4_read_dimj
   
! ***************************************************************
! * subroutine ionc4_litdimk                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    10/12/03                                *
! * derniere modif :    122/12/10 par rramel (ALYOTECH)                 *
! *      *
! * Role : lit la taille de la grille en k                      *
! *      *
! * Parametres :      *
! *       entree:- nomfich   :  nom du fichier            *
! *       sortie:- dimk  :  taille de la grille en k        *
! ***************************************************************
   subroutine ionc4_read_dimk(nom_fichier,dimk)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   integer :: dimk
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: dim_z_id,nc_id,nc_err
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_dimk"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomz),dim_z_id)
      call ionc4_err(nc_err,ionc_rout, 'nf90_inq_dimid',trim(ionc_nomz))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_z_id,len=dimk)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomz))
      
   endif
   
   return
   end subroutine ionc4_read_dimk
   
! ***************************************************************
! * subroutine ionc4_read_diml                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    10/12/03                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                    *
! *      *
! * Role : lit la taille de la grille en l (temps)              *
! *      *
! * Parametres :      *
! *       entree:- nomfich   :  nom du fichier            *
! *       sortie:- dimk  :  taille de la grille en l (temps)     *
! ***************************************************************
   subroutine ionc4_read_diml(nom_fichier,diml)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   integer :: diml
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: dim_l_id,nc_id,nc_err
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_diml"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_l_id)
      call ionc4_err(nc_err,ionc_rout, 'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_l_id,len=diml)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomtime))
      
   endif
   
   return
   end subroutine ionc4_read_diml
   
! ***************************************************************
! * function ionc4_flimax                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    22/02/01                                *
! * derniere modif :   26/09/11 par jfleroux (IFREMER)          *
! *      *
! * Role : lit l'indice max des longitudes                      *
! *      *
! * Parametres :      *
! *       entree:- nc_id :  id du fichier        *
! *             :- varCoord : nom de la variable longitude        *
! * Resultat :      *
! *       indice longitude max*
! ***************************************************************
   function ionc4_flimax(nc_id)
   
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   integer :: nc_id
   integer :: ionc4_flimax
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: var_lon_id,lima,nc_err,nbDims, i, id,iLon
   integer, dimension(nf90_max_var_dims) :: lonDimIds   
   character (len = 80) :: attCoordinates,attTestYX,attTestXY
   logical :: isPresent = .false.
   CHARACTER(LEN=8) :: ionc_tabnomlat(5)
   CHARACTER(LEN=9) :: ionc_tabnomlon(5)

! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
  
   ionc_tabnomlat(1) = "lat     "
   ionc_tabnomlat(2) = "latitude"
   ionc_tabnomlat(3) = "LAT     "
   ionc_tabnomlat(4) = "LATITUDE"
   ionc_tabnomlat(5) = "nav_lat "
   ionc_tabnomlon(1) = "lon      "
   ionc_tabnomlon(2) = "longitude"
   ionc_tabnomlon(3) = "LON      "
   ionc_tabnomlon(4) = "LONGITUDE"
   ionc_tabnomlon(5) = "nav_lon  "
  
   ! recherche du nom utilise pour la variable longitude
   !**************************************
   do i = 1,size(ionc_tabnomlon)
      nc_err = nf90_inq_varid(nc_id,trim(ionc_tabnomlon(i)),id)
      if (nc_err .eq. NF90_NOERR) then
         isPresent = .true.
         iLon=i
         var_lon_id = id
         exit
      endif 
   end do

   if (.NOT. isPresent ) then
      write(*,*) "None name for the longitude variable is correct !!"
      write(*,*) "They are different from lon,longitude,LON,LONGITUDE"
      stop
   endif
   
   ! recherche du nombre de dimensions de la variable ionc_tabnomlon(iLon)
   !**************************************************************
   nc_err = nf90_inquire_variable(nc_id,var_lon_id,ndims=nbDims)
   call ionc4_err(nc_err,'ionc4_flimax', 'nf90_inquire_variable',trim(ionc_tabnomlon(iLon)))
   
   ! recherche des IDs des dimensions de la variable ionc_tabnomlon(iLon)
   !**************************************************************
   nc_err = nf90_inquire_variable(nc_id,var_lon_id,dimids=lonDimIds(:nbDims))
   call ionc4_err(nc_err,'ionc4_flimax', 'nf90_inquire_variable',trim(ionc_tabnomlon(iLon)))

   if (nbDims .eq. 1) then   
      ! si une seule dimension, on recherche la taille de cette dimension
      !******************************************************************
      nc_err = nf90_inquire_dimension(nc_id,lonDimIds(nbDims),len=lima)
      call ionc4_err(nc_err,'ionc4_flimax','nf90_inquire_dimension',trim(ionc_tabnomlon(iLon)))
   else if (nbDims .eq. 2) then
      ! si 2 dimensions, on recherche si l'attribut coordinates existe pour la variable ionc_tabnomlon(iLon)
      !*********************************************************************************************
      nc_err = nf90_inquire_attribute(nc_id,var_lon_id , "coordinates")
      if (nc_err .eq. NF90_NOERR) then
         ! on recupere la valeur de l'attribut coordinates de la variable ionc_tabnomlon(iLon)
         !****************************************************************************
         nc_err = nf90_get_att(nc_id,var_lon_id , "coordinates", attCoordinates)
         call ionc4_err(nc_err,'ionc4_flimax','nf90_get_att',"coordinates")
         ! on construit 2 variables de test avec ionc_tabnomlon(iLon) et ionc_tabnomlat(iLon)
         !*******************************************************************
         attTestYX = trim(ionc_tabnomlat(iLon))//' '//trim(ionc_tabnomlon(iLon))
         attTestXY = trim(ionc_tabnomlon(iLon))//' '//trim(ionc_tabnomlat(iLon))
         if ( index(attCoordinates,attTestYX) .eq. 1) then
            ! on est dans le cas 'ionc_tabnomlat(iLon) ionc_tabnomlon(iLon)' , on garde la dim 1
            !*******************************************************************
            nc_err = nf90_inquire_dimension(nc_id,lonDimIds(1),len=lima)
            call ionc4_err(nc_err,'ionc4_flimax','nf90_inquire_dimension',trim(ionc_tabnomlon(iLon)))
         else if ( index(attCoordinates,attTestXY) .eq. 1) then
            ! on est dans le cas 'ionc_tabnomlon(iLon) ionc_tabnomlat(iLon)' , on garde la dim 2
            !*******************************************************************
            nc_err = nf90_inquire_dimension(nc_id,lonDimIds(2),len=lima)
            call ionc4_err(nc_err,'ionc4_flimax','nf90_inquire_dimension',trim(ionc_tabnomlon(iLon)))
         else
            write(*,*) 'ERROR ionc4_flimax: coordinates attribute not correct for var ',trim(ionc_tabnomlon(iLon))
            stop
         endif
       else
            ! on garde la dim 1 de la variable longitude par defaut
            !*******************************************************************
            nc_err = nf90_inquire_dimension(nc_id,lonDimIds(1),len=lima)
            call ionc4_err(nc_err,'ionc4_flimax','nf90_inquire_dimension',trim(ionc_tabnomlon(iLon)))
       endif
   else
      write(*,*) 'ERROR ionc4_flimax: number of dimensions not correct for var ',trim(ionc_tabnomlon(iLon))
      stop
   endif
   
   ionc4_flimax = lima
   end function ionc4_flimax
   
! ***************************************************************
! * function ionc4_fnlimax                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    28/02/02                                *
! * derniere modif :   26/09/11 par jfleroux (IFREMER)          *
! *      *
! * Role : lit l'indice max des longitudes                      *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! * Resultat :      *
! *       indice longitude max*
! ***************************************************************
   function ionc4_fnlimax(nom_fichier)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   integer :: ionc4_fnlimax 
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id
   
! ******** FONCTIONS **************
   
! ******** FIN DES DECLARATIONS **************
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,'ionc4_fnlimax', ' ' ,nom_fichier)
      stop
   else
      ionc4_fnlimax = ionc4_flimax(nc_id)
   end if
   end function ionc4_fnlimax
   
! ***************************************************************
! * fonction ionc4_fljmax                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    22/02/01                                *
! * derniere modif :   26/09/11 par jfleroux (IFREMER)          *
! *      *
! * Role : lit l'indice max des latitudes                       *
! *      *
! * Parametres :      *
! *       entree:- nc_id :  id du fichier        *
! * Resultat :      *
! *       indice latitude max*
! ***************************************************************
   function ionc4_fljmax(nc_id)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   integer :: nc_id
   integer :: ionc4_fljmax
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: var_lat_id,ljma,nc_err,nbDims, i, id,iLat
   integer, dimension(nf90_max_var_dims) :: latDimIds   
   character (len = 80) :: attCoordinates,attTestYX,attTestXY
   logical :: isPresent = .false.
   
   CHARACTER(LEN=8) :: ionc_tabnomlat(5)
   CHARACTER(LEN=9) :: ionc_tabnomlon(5)
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************

   ionc_tabnomlat(1) = "lat     "
   ionc_tabnomlat(2) = "latitude"
   ionc_tabnomlat(3) = "LAT     "
   ionc_tabnomlat(4) = "LATITUDE"
   ionc_tabnomlat(5) = "nav_lat "
   ionc_tabnomlon(1) = "lon      "
   ionc_tabnomlon(2) = "longitude"
   ionc_tabnomlon(3) = "LON      "
   ionc_tabnomlon(4) = "LONGITUDE"
   ionc_tabnomlon(5) = "nav_lon  "

   ! recherche du nom utilise pour la variable latitude
   !**************************************
   do i = 1,size(ionc_tabnomlat)
      nc_err = nf90_inq_varid(nc_id,trim(ionc_tabnomlat(i)),id)
      if (nc_err .eq. NF90_NOERR) then
         isPresent = .true.
         iLat=i
         var_lat_id = id
         exit
      endif 
   end do

   if (.NOT. isPresent ) then
      write(*,*) "None name for the latitude variable is correct !!"
      write(*,*) "They are different from lat,latitude,LAT,LATITUDE"
      stop
   endif

   ! recherche du nombre de dimensions de la variable ionc_tabnomlat(iLat)
   !**************************************************************
   nc_err = nf90_inquire_variable(nc_id,var_lat_id,ndims=nbDims)
   call ionc4_err(nc_err,'ionc4_fljmax', 'nf90_inquire_variable',trim(ionc_tabnomlat(iLat)))
   
   ! recherche des IDs des dimensions de la variable ionc_tabnomlat(iLat)
   !**************************************************************
   nc_err = nf90_inquire_variable(nc_id,var_lat_id,dimids=latDimIds(:nbDims))
   call ionc4_err(nc_err,'ionc4_fljmax', 'nf90_inquire_variable',trim(ionc_tabnomlat(iLat)))

   if (nbDims .eq. 1) then   
      ! si une seule dimension, on recherche la taille de cette dimension
      !******************************************************************
      nc_err = nf90_inquire_dimension(nc_id,latDimIds(nbDims),len=ljma)
      call ionc4_err(nc_err,'ionc4_fljmax','nf90_inquire_dimension',trim(ionc_tabnomlat(iLat)))
   else if (nbDims .eq. 2) then
      ! si 2 dimensions, on recherche si l'attribut coordinates existe pour la variable ionc_tabnomlat(iLat)
      !*********************************************************************************************
      nc_err = nf90_inquire_attribute(nc_id,var_lat_id , "coordinates")
      if (nc_err .eq. NF90_NOERR) then
         ! on recupere la valeur de l'attribut coordinates de la variable ionc_tabnomlat(iLat)
         !****************************************************************************
         nc_err = nf90_get_att(nc_id,var_lat_id , "coordinates", attCoordinates)
         call ionc4_err(nc_err,'ionc4_fljmax','nf90_get_att',"coordinates")
         ! on construit 2 variables de test avec ionc_tabnomlat(iLat) et ionc_tabnomlat(iLat)
         !*******************************************************************
         attTestYX = trim(ionc_tabnomlat(iLat))//' '//trim(ionc_tabnomlon(iLat))
         attTestXY = trim(ionc_tabnomlon(iLat))//' '//trim(ionc_tabnomlat(iLat))
         if ( index(attCoordinates,attTestYX) .eq. 1) then
            ! on est dans le cas 'ionc_tabnomlat(iLat) ionc_tabnomlon(iLat)' , on garde la dim 2
            !*******************************************************************
            nc_err = nf90_inquire_dimension(nc_id,latDimIds(2),len=ljma)
            call ionc4_err(nc_err,'ionc4_fljmax','nf90_inquire_dimension',trim(ionc_tabnomlat(iLat)))
         else if ( index(attCoordinates,attTestXY) .eq. 1) then
            ! on est dans le cas 'ionc_tabnomlon(iLat) ionc_tabnomlat(iLat)' , on garde la dim 1
            !*******************************************************************
            nc_err = nf90_inquire_dimension(nc_id,latDimIds(1),len=ljma)
            call ionc4_err(nc_err,'ionc4_fljmax','nf90_inquire_dimension',trim(ionc_tabnomlat(iLat)))
         else
            write(*,*) 'ERROR ionc4_fljmax: coordinates attribute not correct for var ',trim(ionc_tabnomlat(iLat))
            stop
         endif
       else
            ! on garde la dim 2 de la variable latitude par defaut
            !*******************************************************************
            nc_err = nf90_inquire_dimension(nc_id,latDimIds(2),len=ljma)
            call ionc4_err(nc_err,'ionc4_fljmax','nf90_inquire_dimension',trim(ionc_tabnomlat(iLat)))
       endif
   else
      write(*,*) 'ERROR ionc4_fljmax: number of dimensions not correct for var ',trim(ionc_tabnomlat(iLat))
      stop
   endif
   
   ionc4_fljmax = ljma
   end function ionc4_fljmax
   
! ***************************************************************
! * fonction ionc4_fnljmax                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    28/02/02                                *
! * derniere modif :   26/09/11 par jfleroux (IFREMER)          *
! *      *
! * Role : lit l'indice max des latitudes                       *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! * Resultat :      *
! *       indice latitude max*
! ***************************************************************
   function ionc4_fnljmax(nom_fichier)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   integer :: ionc4_fnljmax
 
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id
   
! ******** FIN DES DECLARATIONS **************
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,'ionc4_fnljmax', ' ' ,nom_fichier)
      stop
   else
      ionc4_fnljmax = ionc4_fljmax(nc_id)
   end if
   end function ionc4_fnljmax
   
! ***************************************************************
! * fonction ionc4_flkmax                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    22/02/01                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                  *
! *      *
! * Role : lit l'indice max des niveaux                         *
! *      *
! * Parametres :      *
! *       entree:- nc_id :  id du fichier        *
! * Resultat :      *
! *       indice niveau max*
! ***************************************************************
   function ionc4_flkmax(nc_id)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   integer :: nc_id
   integer :: ionc4_flkmax
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: dim_k_id,lkma,nc_err
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   
   nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomz),dim_k_id)
   call ionc4_err(nc_err,'ionc4_flkmax', 'nf90_inq_dimid',trim(ionc_nomz))
   
   nc_err = nf90_inquire_dimension(nc_id,dim_k_id,len=lkma)
   call ionc4_err(nc_err,'ionc4_flkmax','nf90_inquire_dimension',trim(ionc_nomz))
   
   ionc4_flkmax = lkma
   return
   end function ionc4_flkmax
   
! ***************************************************************
! * fonction ionc4_fnlkmax                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    28/02/02                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                   *
! *      *
! * Role : lit l'indice max des niveaux                         *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! * Resultat :      *
! *       indice niveau max*
! ***************************************************************
   function ionc4_fnlkmax(nom_fichier)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   integer :: ionc4_fnlkmax
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: dim_k_id,lkma,nc_err,nc_id
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,'ionc4_fnlkmax', ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomz),dim_k_id)
      call ionc4_err(nc_err,'ionc4_fnlkmax', 'nf90_inq_dimid',trim(ionc_nomz))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_k_id,len=lkma)
      call ionc4_err(nc_err,'ionc4_fnlkmax','nf90_inquire_dimension',trim(ionc_nomz))
      
      ionc4_fnlkmax = lkma
   end if
   end function ionc4_fnlkmax
   
! ***************************************************************
! * function ionc4_fltmax                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    22/02/01                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                  *
! *      *
! * Role : lit l'indice max des temps                           *
! *      *
! * Parametres :      *
! *       entree:- nc_id :  id du fichier        *
! * Resultat :      *
! *       indice temps max*
! ***************************************************************
   function ionc4_fltmax(nc_id)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   integer :: nc_id
   integer :: ionc4_fltmax
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: dim_time_id,nc_err,ltma
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   
   nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
   call ionc4_err(nc_err,'ionc4_fltma', 'nf90_inq_dimid',trim(ionc_nomtime))
   
   nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=ltma)
   call ionc4_err(nc_err,'ionc4_fltma','nf90_inquire_dimension',trim(ionc_nomtime))
   
   ionc4_fltmax = ltma
   
   end function ionc4_fltmax
   
! ***************************************************************
! * function ionc4_fnltmax                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    28/02/02                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                   *
! *      *
! * Role : lit l'indice max des temps                           *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! * Resultat :      *
! *       indice temps max*
! ***************************************************************
   function ionc4_fnltmax(nom_fichier)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   integer :: ionc4_fnltmax
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: dim_time_id,nc_err,ltma,nc_id
   
! ******** FONCTIONS **************
   integer :: trim
   
! ******** FIN DES DECLARATIONS **************
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,'ionc4_fnltmax', ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomtime),dim_time_id)
      call ionc4_err(nc_err,'ionc4_fnltmax', 'nf90_inq_dimid',trim(ionc_nomtime))
      
      nc_err = nf90_inquire_dimension(nc_id,dim_time_id,len=ltma)
      call ionc4_err(nc_err,'ionc4_fnltmax','nf90_inquire_dimension',trim(ionc_nomtime))
      
      ionc4_fnltmax = ltma
   end if
   end function ionc4_fnltmax
 
! ***************************************************************
! * subroutine ionc4_write_z_real                   *
! *      *
! * auteur         :    jfleroux,vgarnier                       *
! * org            :    IFREMER                                 *
! * date creation  :    26/08/02                                *
! *      *
! * Role : ecrit une variable reelle 1D =  z(kmax)              *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - var         :  variable                *
! *       - kmax        :  borne max dim k*
! ***************************************************************
   subroutine ionc4_write_z_real(nom_fichier,var_name,var,kmax,fill_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   
   character(len=*) :: nom_fichier,var_name
   integer :: kmax
   real(kind=4),dimension(kmax) :: var
   real(kind=4),optional,intent(in) :: fill_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,dimk_id,nc_err,dimk,k,indk
   real(kind=4) :: add_offset,scale_factor
   integer,dimension(1) :: start,count
   real(kind=4),dimension(:),allocatable :: bid
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_write_z_real"

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else

      nc_err = nf90_inq_dimid(nc_id,trim(ionc_nomz),dimk_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_dimid',trim(ionc_nomz))

      nc_err = nf90_inquire_dimension(nc_id,dimk_id,len=dimk)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_dimension',trim(ionc_nomz))

      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif

      allocate(bid(dimk))
      start  = (/1/)
      count = (/dimk/)

      if (scale_factor /= 1.e0 .and. present(fill_value)) then
         indk = 1
         !do k=ionc_lkmi,ionc_lkma, ionc_pask
         do k=1,dimk
            if (var(k) /= fill_value) then
              bid(indk)=real(nint((var(k)- add_offset)/scale_factor))
            else
               ! Promotion to integer2 when writing since the variable is saved in short if l_out_pack
               bid(indk)=-HUGE(1_2)-1.0_4
            endif
            indk = indk + 1
         enddo
      else
         indk = 1
         !do k=ionc_lkmi,ionc_lkma, ionc_pask
         do k=1,dimk
            bid(indk)=var(k)
            indk = indk + 1
         enddo
      endif

      nc_err = nf90_put_var(nc_id,var_id,bid,start=start,count=count)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_var',trim(var_name))
      deallocate(bid)
   endif
   
   return
   end subroutine ionc4_write_z_real
   
! ***************************************************************
! * subroutine ionc4_read_z_real                         *
! *      *
! * auteur         :    jfleroux,vgarnier                       *
! * org            :    IFREMER                                 *
! * date creation  :    26/08/02                                *
! *      *
! * Role : lit une variable reelle 1D =  z                       *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - kmax        :  borne max dim k*
! *     sortie: - var         :  variable                      *
! ***************************************************************
   subroutine ionc4_read_z_real(nom_fichier,var_name,var,kmin,kmax,fillvar)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier,var_name
   integer :: kmin,kmax
   real(kind=4) :: var(kmin:kmax)
   real(kind=4),optional,intent(in) :: fillvar
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,k,nc_err_fill
   real(kind=4) :: add_offset,scale_factor,fillvalue_read
   integer ,dimension(1):: start,count
   real(kind=4) :: bid(kmin:kmax)
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_read_z_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor', scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset', add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      nc_err_fill = nf90_get_att(nc_id,var_id,'_FillValue', fillvalue_read)
      
      start(1)= 1
      count(1) = kmax-kmin+1
      
      nc_err = nf90_get_var(nc_id,var_id,bid,start=start,count=count)
      
      do k = kmin,kmax
         var(k)=bid(k)*scale_factor+add_offset
      end do

      if (present(fillvar) .and. nc_err_fill.eq.nf90_noerr) then
         where (var(:)==fillvalue_read*scale_factor+add_offset ) var(:)=fillvar
      endif

   endif
   
   return
   end subroutine ionc4_read_z_real

! ***************************************************************
! * subroutine ionc4_showvar                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    26/10/02                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                    *
! *      *
! * Role : affichage des variables du fichier                   *
! *      sans les variables lat, lon, z et time*
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! * Resultat :      *
! ***************************************************************
   subroutine ionc4_showvar(nom_fichier)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: nom_fichier
   
! ******** VARIABLES DE TRAVAIL **************
   character(len=30) :: var_name
   integer :: nvars, ndims, ngatts,xdimid
   integer :: var_id, var_type, var_ndims, var_natts
   integer,dimension(10) :: var_ndim
   integer :: nc_err,nc_id
   
! ******** FIN DES DECLARATIONS **************
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,'ionc4_showvar', ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_inquire(nc_id, ndims, nvars, ngatts, xdimid)
      
      if (nvars .gt. 0) then
         write (*,*) "variables:"
         do var_id = 1, nvars
            nc_err = nf90_inquire_variable(nc_id, var_id, var_name,var_type, var_ndims,var_ndim, var_natts)
            if ((var_name.NE.ionc_pnomlat).AND.   &
                (var_name.NE.ionc_pnomlon).AND.   &
                (var_name.NE.ionc_pnomz)  .AND.   &
                (var_name.NE.ionc_pnomtime)) then
               write (*,*) var_name
            endif
         end do
      end if
   end if
   end subroutine ionc4_showvar
   
! ***************************************************************
! * subroutine ionc4_version                          *
! *      *
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    25/06/07                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)     *
! *      *
! * Role : donne le numero de version de la librairie IONETCDF  *
! *      *
! ***************************************************************
   subroutine ionc4_version(version)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) :: version
   
! ******** VARIABLES DE TRAVAIL **************
   
! ******** FIN DES DECLARATIONS **************
!tina       version = ionc4_pversion
   version = '01/01/2011'
   return
   end subroutine ionc4_version
   
! ***************************************************************
! * subroutine ionc4_vatt_axis                   *
! *      *
! * auteur         :    ggaliber (ALTRAN)                 *
! * org            :    ALTRAN                                 *
! * date creation  :    22/09/08                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)    *
! *      *
! * Role : ajoute l'attribut axis a une variable          *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - axis_value   :  valeur de axis          *
! ***************************************************************
   subroutine ionc4_vatt_axis(nom_fichier,var_name,axis_value)
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name,axis_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_axis"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_put_att(nc_id, var_id,'axis', axis_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:axis',var_name)
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_axis
   
! ***************************************************************
! * subroutine ionc4_vatt_positive                   *
! *      *
! * auteur         :    ggaliber (ALTRAN)                 *
! * org            :    ALTRAN                                 *
! * date creation  :    22/09/08                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)    *
! *      *
! * Role : ajoute l'attribut positive a une variable (depth par ex.)          *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - positive_value   :  valeur de conventions          *
! ***************************************************************
   subroutine ionc4_vatt_positive(nom_fichier,var_name,positive_value)
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name,positive_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_positive"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_put_att(nc_id, var_id,'positive', positive_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:positive',var_name)
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_positive
   
! ***************************************************************
! * subroutine ionc4_vatt_formula                                *
! *                                                             *
! * auteur         :    V. Garnier                              *
! * org            :    IFREMER                                 *
! * date creation  :    18/12/08                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                 *
! *                                                             *
! * Role : ajoute l'attribut formula a une variable             *
! *                                                             *
! * Parametres :                                                *
! *       entree:- nom_fichier :  nom du fichier                *
! *       - var_name    :  nom de la variable                   *
! *       - formula_terms   :  valeur de conventions            *
! ***************************************************************
   subroutine ionc4_vatt_formula(nom_fichier,var_name,formula_terms)
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name,formula_terms
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_formula"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_put_att(nc_id, var_id,'formula_terms', formula_terms)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:formula_terms',var_name)
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_formula

! ***************************************************************
! * subroutine ionc4_vatt_chain                                 *
! *                                                             *
! * auteur         :    V. Garnier                              *
! * org            :    IFREMER                                 *
! * date creation  :    22/06/12                                *
! * derniere modif :                                            *
! *                                                             *
! * Role : ajoute l'attribut comment a une variable             *
! *                                                             *
! * Parametres :                                                *
! *       entree:- nom_fichier :  nom du fichier                *
! *       - var_name    :  nom de la variable                   *
! *       - att_name    :  nom de l attribut a introduire       *
! *       - att_terms   :  chaine de caractere a introduire     *
! ***************************************************************
   subroutine ionc4_vatt_chain(nom_fichier,var_name,att_name,att_terms)
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name,att_name,att_terms

! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err

! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_chain"

   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')

      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)

      nc_err = nf90_put_att(nc_id, var_id,trim(att_name), trim(att_terms))
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:comment',var_name)

      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif

   return
   end subroutine ionc4_vatt_chain
 
! ***************************************************************
! * subroutine ionc4_vatt_coordinate                             *
! *                                                             *
! * auteur         :    V. Garnier                              *
! * org            :    IFREMER                                 *
! * date creation  :    18/12/08                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)                *
! *                                                             *
! * Role : ajoute l'attribut standard_name a une variable       *
! *        pour specifier le type de coordonnees verticles      *
! *                                                             *
! * Parametres :                                                *
! *       entree:- nom_fichier :  nom du fichier                *
! *       - var_name    :  nom de la variable                   *
! *       - coordinate_type   :  valeur de conventions            *
! ***************************************************************
   subroutine ionc4_vatt_coordinate(nom_fichier,var_name,coordinate_type)
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name,coordinate_type
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_coordinate"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_put_att(nc_id, var_id,'standard_name',coordinate_type)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:standard_name',var_name)
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_coordinate
   
! ***************************************************************
! * subroutine ionc4_vatt_conventions                   *
! *      *
! * auteur         :    ggaliber (ALTRAN)                 *
! * org            :    ALTRAN                                 *
! * date creation  :    22/09/08                                *
! * derniere modif :    22/09/08 par ggaliber (ALTRAN)    *
! *                     2 5/03/09 vgarnier (Maj a Conventions)     *
! *                     22/12/10 par rramel (ALYOTECH)   (F90)* 
! * Role : ajoute l'attribut conventions a une variable (time par ex.)          *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - conventions_value   :  valeur de conventions          *
! ***************************************************************
   subroutine ionc4_vatt_conventions(nom_fichier,var_name,conventions_value)
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name,conventions_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_conventions"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_put_att(nc_id, var_id,'conventions', conventions_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:conventions',var_name)
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_conventions
   
! ***************************************************************
! * subroutine ionc4_vatt_valid_min_double                   *
! *      *
! * auteur         :    ggaliber (ALTRAN)                 *
! * org            :    IFREMER                                 *
! * date creation  :    22/09/08                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)    *
! *      *
! * Role : ajoute l attribut valid_min a une variable          *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - valid_min_value   :  valeur de valid_min          *
! ***************************************************************
   subroutine ionc4_vatt_valid_min_double(nom_fichier,var_name,valid_min_value)
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   real(kind=8) :: valid_min_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,nc_type
   real(kind=4) :: scale_factor,add_offset
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_valid_min_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
     
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      nc_err = nf90_inquire_variable(nc_id,var_id,xtype=nc_type)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_variable',var_name)      

      if (nc_type .eq. nf90_short) then
         if (scale_factor /= 1.e0 .and. add_offset /= 0.e0) then
            nc_err = nf90_put_att(nc_id,var_id,'valid_min', nint((valid_min_value-add_offset)/scale_factor,2))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
         else
            nc_err = nf90_put_att(nc_id,var_id,'valid_min', nint(valid_min_value,2))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
         endif   
      elseif (nc_type .eq. nf90_int) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_min', nint(valid_min_value))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
      elseif (nc_type .eq. nf90_real) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_min', real(valid_min_value))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
      elseif (nc_type .eq. nf90_double) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_min', valid_min_value)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
      endif
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_valid_min_double
   
! ***************************************************************
! * subroutine ionc4_vatt_valid_min_real                   *
! *      *
! * auteur         :    ggaliber (ALTRAN)                 *
! * org            :    IFREMER                                 *
! * date creation  :    22/09/08                                *
! * derniere modif :   22/12/10 par rramel (ALYOTECH)     *
! *      *
! * Role : ajoute l attribut valid_min a une variable          *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - valid_min_value   :  valeur de valid_min          *
! ***************************************************************
   subroutine ionc4_vatt_valid_min_real(nom_fichier,var_name,valid_min_value)
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   real(kind=4) :: valid_min_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,nc_type
   real(kind=4) :: scale_factor,add_offset
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_valid_min_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      nc_err = nf90_inquire_variable(nc_id,var_id,xtype=nc_type)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_variable',var_name)      

      if (nc_type .eq. nf90_short) then
         if (scale_factor /= 1.e0 .and. add_offset /= 0.e0) then
            nc_err = nf90_put_att(nc_id,var_id,'valid_min', nint((valid_min_value-add_offset)/scale_factor,2))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
         else
            nc_err = nf90_put_att(nc_id,var_id,'valid_min', nint(valid_min_value,2))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
         endif   
      elseif (nc_type .eq. nf90_int) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_min', nint(valid_min_value))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
      elseif (nc_type .eq. nf90_real) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_min', valid_min_value)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
      elseif (nc_type .eq. nf90_double) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_min', REAL(valid_min_value,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_min',var_name)
      endif
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_valid_min_real
   
! ***************************************************************
! * subroutine ionc4_vatt_valid_max_double                   *
! *      *
! * auteur         :    ggaliber (ALTRAN)                 *
! * org            :    IFREMER                                 *
! * date creation  :    22/09/08                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)    *
! *      *
! * Role : ajoute l attribut valid_max a une variable          *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - valid_min_value   :  valeur de valid_max          *
! ***************************************************************
   subroutine ionc4_vatt_valid_max_double(nom_fichier,var_name,valid_max_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   real(kind=8) :: valid_max_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,nc_type
   real(kind=4) :: scale_factor,add_offset
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_valid_max_double"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      nc_err = nf90_inquire_variable(nc_id,var_id,xtype=nc_type)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_variable',var_name)      

      if (nc_type .eq. nf90_short) then
         if (scale_factor /= 1.e0 .and. add_offset /= 0.e0) then
            nc_err = nf90_put_att(nc_id,var_id,'valid_max', nint((valid_max_value-add_offset)/scale_factor,2))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
         else
            nc_err = nf90_put_att(nc_id,var_id,'valid_max', nint(valid_max_value,2))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
         endif   
      elseif (nc_type .eq. nf90_int) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_max', nint(valid_max_value))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
      elseif (nc_type .eq. nf90_real) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_max', real(valid_max_value))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
      elseif (nc_type .eq. nf90_double) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_max', valid_max_value)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
      endif
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_valid_max_double
   
! ***************************************************************
! * subroutine ionc4_vatt_valid_max_real                   *
! *      *
! * auteur         :    ggaliber (ALTRAN)                 *
! * org            :    IFREMER                                 *
! * date creation  :    22/09/08                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)  *
! *      *
! * Role : ajoute l attribut valid_max a une variable          *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - valid_max_value   :  valeur de valid_max          *
! ***************************************************************
   subroutine ionc4_vatt_valid_max_real(nom_fichier,var_name,valid_max_value)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name
   real(kind=4) :: valid_max_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err,nc_type
   real(kind=4) :: scale_factor,add_offset
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_valid_max_real"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_get_att(nc_id,var_id,'scale_factor',scale_factor)
      if (nc_err.ne.nf90_noerr) then
         scale_factor = 1.e0
      endif
      nc_err = nf90_get_att(nc_id,var_id,'add_offset',add_offset)
      if (nc_err.ne.nf90_noerr) then
         add_offset = 0.e0
      endif
      nc_err = nf90_inquire_variable(nc_id,var_id,xtype=nc_type)
      call ionc4_err(nc_err,ionc_rout,'nf90_inquire_variable',var_name)      

      if (nc_type .eq. nf90_short) then
         if (scale_factor /= 1.e0 .and. add_offset /= 0.e0) then
            nc_err = nf90_put_att(nc_id,var_id,'valid_max', nint((valid_max_value-add_offset)/scale_factor,2))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
         else
            nc_err = nf90_put_att(nc_id,var_id,'valid_max', nint(valid_max_value,2))
            call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
         endif   
      elseif (nc_type .eq. nf90_int) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_max', nint(valid_max_value))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
      elseif (nc_type .eq. nf90_real) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_max', valid_max_value)
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
      elseif (nc_type .eq. nf90_double) then
         nc_err = nf90_put_att(nc_id,var_id,'valid_max',REAL(valid_max_value,8))
         call ionc4_err(nc_err,ionc_rout,'nf90_put_att:valid_max',var_name)
      endif
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_valid_max_real
   
! ***************************************************************
! * subroutine ionc4_vatt_standard_name                   *
! *      *
! * auteur         :    ggaliber (ALTRAN)                 *
! * org            :    ALTRAN                                 *
! * date creation  :    22/09/08                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)  *
! *      *
! * Role : ajoute l'attribut standard_name a une variable          *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - standard_name_value   :  valeur de standard_name          *
! ***************************************************************
   subroutine ionc4_vatt_standard_name(nom_fichier,var_name,standard_name_value)
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name,standard_name_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_standard_name"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_put_att(nc_id, var_id,'standard_name', standard_name_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:standard_name',var_name)
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_standard_name
   
! ***************************************************************
! * subroutine ionc4_vatt_long_name                   *
! *      *
! * auteur         :    ggaliber (ALTRAN)                 *
! * org            :    ALTRAN                                 *
! * date creation  :    22/09/08                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)  *
! *      *
! * Role : ajoute l'attribut long_name a une variable          *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - var_name    :  nom de la variable        *
! *       - long_name_value   :  valeur de long_name          *
! ***************************************************************
   subroutine ionc4_vatt_long_name(nom_fichier,var_name,long_name_value)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,var_name,long_name_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_vatt_long_name"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_inq_varid',var_name)
      
      nc_err = nf90_put_att(nc_id, var_id,'long_name', long_name_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att:long_name',var_name)
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',var_name)
   endif
   
   return
   end subroutine ionc4_vatt_long_name
   
! ***************************************************************
! * subroutine ionc4_gatt_char                   *
! *      *
! * auteur         :    ggaliber (ALTRAN)                 *
! * org            :    ALTRAN                                 *
! * date creation  :    22/09/08                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)  *
! *      *
! * Role : ajoute un attribut global au fichier netCDF           *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - gatt_name    :  nom de l'attribut global        *
! *       - gatt_value   :  valeur de l'attribut global          *
! ***************************************************************
   subroutine ionc4_gatt_char(nom_fichier,gatt_name,gatt_value)
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,gatt_name,gatt_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_gatt_char"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_put_att(nc_id, nf90_global,gatt_name, gatt_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',gatt_name)
   endif
   
   return
   end subroutine ionc4_gatt_char
   
! ***************************************************************
! * subroutine ionc4_gatt_int                   *
! *      *
! * auteur         :    jfleroux                 *
! * org            :    IFREMER                                 *
! * date creation  :    06/09/09                                *
! * derniere modif :    22/12/10 par rramel (ALYOTECH)    *
! *      *
! * Role : ajoute un attribut global de type int           *
! *      *
! * Parametres :      *
! *       entree:- nom_fichier :  nom du fichier        *
! *       - gatt_name    :  nom de l'attribut global        *
! *       - gatt_value   :  valeur de l'attribut global          *
! ***************************************************************
   subroutine ionc4_gatt_int(nom_fichier,gatt_name,gatt_value)
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier,gatt_name
   integer ::  gatt_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_gatt_int"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      nc_err = nf90_put_att(nc_id, nf90_global,gatt_name,gatt_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',gatt_name)
   endif
   
   return
   end subroutine ionc4_gatt_int
   
! ***************************************************************
! * subroutine ionc4_gatt_conv                                   *
! *                                                             *
! * auteur         :    V. Garnier                              *
! * org            :    IFREMER                                 *
! * date creation  :    18/12/08                                *
! * derniere modif :    18/12/08 par V. Garnier                 *
! *                                                             *
! * Role : ajoute un attribut global au fichier netCDF          *
! *                                                             *
! * Parametres :                                                *
! *       entree:- nom_fichier :  nom du fichier                *
! *       - gatt_name    :  nom de l'attribut global            *
! *       - gatt_value   :  valeur de l'attribut global         *
! ***************************************************************
   subroutine ionc4_gatt_conv(nom_fichier)

! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*) ::   nom_fichier
   character(len=100)  ::gatt_name,gatt_value
   
! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_gatt_conv"
   
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then
      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop
   else
      nc_err = nf90_redef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_redef',' ')
      
      gatt_name='data_type'
      gatt_value=ionc_conv_data_type
      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
      gatt_name='format_version'
      gatt_value=ionc_conv_format_version
      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
      gatt_name='Conventions'
      gatt_value=ionc_conv_conventions
      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
      gatt_name='netcdf_version'
      gatt_value=ionc_conv_netcdf_version
      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
      gatt_name='product_version'
      gatt_value=ionc_conv_product_version
      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
      gatt_name='software_version'
      gatt_value=ionc_conv_software_version
      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
      gatt_name='references'
      gatt_value=ionc_conv_references
      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
      gatt_name='easting'
      gatt_value=ionc_conv_easting
      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
      gatt_name='northing'
      gatt_value=ionc_conv_northing
      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
      gatt_name='grid_projection'
      gatt_value=ionc_conv_grid_projection
      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
      gatt_name='contact'
      gatt_value=ionc_conv_contact
      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
!      gatt_name='minimum_altitude'
!      gatt_value=ionc_conv_minimum_altitude
!      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
!      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
!      gatt_name='maximum_altitude'
!      gatt_value=ionc_conv_maximum_altitude
!      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
!      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
!      gatt_name='altitude_resolution'
!      gatt_value=ionc_conv_altitude_resolution
!      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
!      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
      gatt_name='distribution_statement'
      gatt_value=ionc_conv_distribution
      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
      gatt_name='quality_index'
      gatt_value=ionc_conv_quality_index
      nc_err = nf90_put_att(nc_id,nf90_global,gatt_name,gatt_value)
      call ionc4_err(nc_err,ionc_rout,'nf90_put_att',gatt_name)
      
      nc_err = nf90_enddef(nc_id)
      call ionc4_err(nc_err,ionc_rout,'nf90_enddef',gatt_name)
   endif
   
   return
   end subroutine ionc4_gatt_conv

! ***************************************************************
! * subroutine ionc4_var_exists                                 *
! *                                                             *
! * auteur         :    S. PETTON                               *
! * org            :    IFREMER                                 *
! * date creation  :    14/09/11                                *
! *                                                             *
! * Role : Verifie l'existence d'une variable                   *
! *                                                             *
! * Parametres :                                                *
! *       entree:- nom_fichier : nom du fichier NETCDF          *
! *              - var_name    : nom de la variable NETCDF      *
! *       sortie:- exists : booleen .true. si existe            *
! ***************************************************************
   subroutine ionc4_var_exists(nom_fichier,var_name,exists)
   
! ******** PARAMETRES DE LA SUBROUTINE *******
   character(len=*),intent(in) :: nom_fichier,var_name
   logical, intent(out)        :: exists

! ******** VARIABLES DE TRAVAIL **************
   integer :: nc_id,var_id,nc_err
   
! ******** FIN DES DECLARATIONS **************
   ionc_rout = "ionc4_var_exists"
   
   exists=.false.
   call ionc4_corres(nom_fichier,nc_id)
   if (nc_id .eq. 0) then

      call ionc4_err(ionc_errfich,ionc_rout, ' ' ,nom_fichier)
      stop

   else

      nc_err = nf90_inq_varid(nc_id,var_name,var_id)
      if(nc_err .eq. nf90_noerr) then
          exists=.true.
      endif  

   endif

   return
   end subroutine ionc4_var_exists
   
   end module ionc4
