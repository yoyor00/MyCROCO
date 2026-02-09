#include "cppdefs.h"
                                      ! correlation
      subroutine das_write_corr  ! to binary file
      implicit none
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
      integer i,j,k,ios,lenstr,len_corr,len_corr_vert
!
! compute the entire model domain
! but the lower triangular matrices
!
# ifdef EW_PERIODIC
#  define IR_ENTIRE 1,LLm
#  define IU_ENTIRE 1,LLm
# else
#  define IR_ENTIRE 0,LLm+1
#  define IU_ENTIRE 1,LLm+1
# endif

# ifdef NS_PERIODIC
#  define JR_ENTIRE 1,MMm
#  define JV_ENTIRE 1,MMm
# else
#  define JR_ENTIRE 0,MMm+1
#  define JV_ENTIRE 1,MMm+1
# endif

!      
      len_corr=lenstr(file_corr)
      len_corr_vert=lenstr(file_corr_vert)
!
! 1D horizontal correlations
!
      open(81, file=file_corr(1:len_corr),form='unformatted',
     &              STATUS='unknown',ACTION='write',IOSTAT=ios)
      if(ios.ne.0) then
        write(*,'(6x,A,(/6x,A))') '??? WARNING: can not find ',
     &           file_corr(1:len_corr)
        return
      endif
!
! temp
!
      write(81) ((ctX_das(i,j,itemp),
     &             i=IR_ENTIRE),j=IR_ENTIRE)
      write(81) ((ctE_das(i,j,itemp),
     &             i=JR_ENTIRE),j=JR_ENTIRE)
!
! salt
!
      write(81) ((ctX_das(i,j,isalt),
     &             i=IR_ENTIRE),j=IR_ENTIRE)
      write(81) ((ctE_das(i,j,isalt),
     &             i=JR_ENTIRE),j=JR_ENTIRE)
!
! psi
!
      write(81) ((cpsiX_das(i,j),
     &             i=IR_ENTIRE),j=IR_ENTIRE)
      write(81) ((cpsiE_das(i,j),
     &             i=JR_ENTIRE),j=JR_ENTIRE)
!
! chi
!
      write(81) ((cchiX_das(i,j),
     &             i=IR_ENTIRE),j=IR_ENTIRE)
      write(81) ((cchiE_das(i,j),
     &             i=JR_ENTIRE),j=JR_ENTIRE)
!
! zeta
!
      write(81) ((czX_das(i,j),
     &             i=IR_ENTIRE),j=IR_ENTIRE)
      write(81) ((czE_das(i,j),
     &             i=JR_ENTIRE),j=JR_ENTIRE)
!
! cross T/S
!
!      write(81) (((corr_ts(i,j,k),
!     &          i=IR_ENTIRE),j=JR_ENTIRE),k=1,ndas)
!
!      close(81)
!
! vertical correlation, cholesky factorization
!
      open(81, file=file_corr_vert(1:len_corr_vert),
     &              STATUS='unknown',ACTION='write',IOSTAT=ios)
      if(ios.ne.0) then
        write(*,'(6x,A,(/6x,A))') '??? WARNING: can not find ',
     &           file_corr_vert(1:len_corr_vert)
        return
      endif
!
! temp
!
       write(81) ((ctZ_das(i,j,itemp),
     &           i=1,ndas),j=1,ndas)
!
! salt
!
      write(81) ((ctZ_das(i,j,isalt),
     &           i=1,ndas),j=1,ndas)
!
! psi & chi
!
      write(81) ((cpsiZ_das(i,j),
     &           i=1,ndas),j=1,ndas)
      write(81) ((cchiZ_das(i,j),
     &           i=1,ndas),j=1,ndas)
!
! u & v
!
      write(81) ((cuZ_das(i,j),
     &           i=1,ndas),j=1,ndas)
      write(81) ((cvZ_das(i,j),
     &           i=1,ndas),j=1,ndas)
!
      close(81)
      return
      end
