
# define WORK UFx
      call step3d_fbijk_nbq (Istr,Iend,Jstr,Jend, WORK
     &     ,Hzw_half_nbq_inv,Hzr_half_nbq_inv
     &     ,Hzw_half_nbq_inv_u,Hzw_half_nbq_inv_v
     &     ,work_nbq(PRIVATE_2D_SCRATCH_ARRAY,1)
     &     ,work_nbq(PRIVATE_2D_SCRATCH_ARRAY,3)
     &     ,work_nbq(PRIVATE_2D_SCRATCH_ARRAY,5)
     &     ,work_nbq(PRIVATE_2D_SCRATCH_ARRAY,7)
     &     ,work_nbq(PRIVATE_2D_SCRATCH_ARRAY,9)
     &     ,work_nbq(PRIVATE_2D_SCRATCH_ARRAY,10)
     &     ,work_nbq(PRIVATE_2D_SCRATCH_ARRAY,11)
!     &   ,Hzu_half_qdmu, Hzv_half_qdmv                   
     &     ,work3d_nbq(PRIVATE_2D_SCRATCH_ARRAY,1,1)
     &     ,work3d_nbq(PRIVATE_2D_SCRATCH_ARRAY,1,2)
!    &     ,work_nbq(START_2D_ARRAY,1), work_nbq(START_2D_ARRAY,3)
!    &     ,work_nbq(START_2D_ARRAY,5), work_nbq(START_2D_ARRAY,7)
!    &     ,work_nbq(START_2D_ARRAY,9), work_nbq(START_2D_ARRAY,10)
!    &     ,work_nbq(START_2D_ARRAY,11)
!    &     ,work3d_nbq(START_2D_ARRAY,1,1)
!    &     ,work3d_nbq(START_2D_ARRAY,1,2)
     &     )
# undef WORK       
