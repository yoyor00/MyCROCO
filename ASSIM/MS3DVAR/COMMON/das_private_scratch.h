c     The following private scratch parameters are 
c     defined in Croco param.h, with different 
c     size => modifying namespace
      integer NNSA, NN2d, NN3d
      parameter (NNSA=8)
c     The following parameters are consistently
c     defined in Croco param.h
c     integer size_XI,size_ETA
c---#define ALLOW_SINGLE_BLOCK_MODE
c-#ifdef  ALLOW_SINGLE_BLOCK_MODE
c      parameter (size_XI=6+Lm, size_ETA=6+Mm)
c-#else
c      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X,
c     &           size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
c-#endif
c     The following private scratch parameters are 
c     defined in param.h with different size 
c     modifying namespace
      integer see,ssee, szz,sszz
      parameter (ssee=size_ETA/NDASp, sszz=NDASp/size_ETA)
      parameter (see=ssee/(ssee+sszz),   szz=1-see)

      parameter (NN2d=size_XI*(see*size_ETA+szz*NDASp),
     &                   NN3d=size_XI*size_ETA*NDASp)
#ifdef SGI
      real A2d(NN2d,NNSA,0:NPP-1), A3d(NN3d,2,0:NPP-1)
      common /das_private_scratch_A2d/A2d
     &       /das_private_scratch_A3d/A3d
#elif defined CRAY
      real A2d(NN2d,NNSA,0:0), A3d(NN3d,2,0:0)
      task common /das_private_scratch/ A2d,A3d
#else
      real A2d(NN2d,NNSA,0:NPP-1), A3d(NN3d,2,0:NPP-1)
      common /das_private_scratch/A2d,A3d
#endif
