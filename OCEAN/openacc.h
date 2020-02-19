      integer, parameter :: sync_ruv_nbq_avg1     = 1
      integer, parameter :: sync_ruv_int_nbq      = 2
      integer, parameter :: sync_rho_rufrc_z_w    = 3
      integer, parameter :: sync_Hz               = 4
      integer, parameter :: sync_ruv_nbq_avg2     = 5
      integer, parameter :: sync_DU_avg1_avg2     = 6

      integer            :: ngpus
      integer            :: devicenum
      integer            :: device_type


#define ACC_TILE_X  Lm     
#define ACC_TILE_Y  Mm    
#define ACC_TILE_Z  N     
