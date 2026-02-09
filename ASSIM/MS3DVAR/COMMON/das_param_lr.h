! Dimensions of Physical Grid and array dimensions: 
! =========== == ======== ==== === ===== ============
! LLm_lr,MMm_lr  Number of the internal points of the PHYSICAL grid.
!          in the XI- and ETA-directions [physical side boundary
!          points and peroodic ghost points (if any) are excluded].
!
! Lm_lr,Mm_lr    Number of the internal points [see above] of array
!          covering a Message Passing subdomain. In the case when
!          no Message Passing partitioning is used, these two are
!          the same as LLm,MMm. 
!
      integer  nratio,nhalf ,  LLm_lr,Lm_lr,  MMm_lr,Mm_lr
      parameter (nratio=1, nhalf=nratio/2+1) 
      parameter (LLm_lr=(LLm+2-nhalf)/nratio-1,
     &           MMm_lr=(MMm+2-nhalf)/nratio-1)
!
      parameter ( Lm_lr=LLm_lr, Mm_lr=MMm_lr)

#define GLOBAL_2D_ARRAY_LR 0:Lm_lr+1+padd_X,0:Mm_lr+1+padd_E
