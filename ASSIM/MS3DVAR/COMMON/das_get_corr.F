#include "cppdefs.h"
#if defined DAS_BVAR_CORR 
                                      ! correlation
      subroutine das_get_corr  ! from binary file
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

#if defined DAS_COUNT_PSICHI
      integer :: count_non_null_psi
      integer :: count_non_null_chi
#endif

!      
      len_corr=lenstr(file_corr)
      len_corr_vert=lenstr(file_corr_vert)
!
! 1D horizontal correlations
!
      open(81, file=file_corr(1:len_corr),form='unformatted',
     &              STATUS='old',IOSTAT=ios)
      if(ios.ne.0) then
        write(*,'(6x,A,(/6x,A))') '??? WARNING: can not find ',
     &           file_corr(1:len_corr)
        return
      endif
!
! temp
!
      read(81) ((ctX_das(i,j,itemp),
     &             i=IR_ENTIRE),j=IR_ENTIRE)
      read(81) ((ctE_das(i,j,itemp),
     &             i=JR_ENTIRE),j=JR_ENTIRE)
!
! salt
!
      read(81) ((ctX_das(i,j,isalt),
     &             i=IR_ENTIRE),j=IR_ENTIRE)
      read(81) ((ctE_das(i,j,isalt),
     &             i=JR_ENTIRE),j=JR_ENTIRE)
!
! psi
!
      read(81) ((cpsiX_das(i,j),
     &             i=IR_ENTIRE),j=IR_ENTIRE)
      read(81) ((cpsiE_das(i,j),
     &             i=JR_ENTIRE),j=JR_ENTIRE)
!
! chi
!
      read(81) ((cchiX_das(i,j),
     &             i=IR_ENTIRE),j=IR_ENTIRE)
      read(81) ((cchiE_das(i,j),
     &             i=JR_ENTIRE),j=JR_ENTIRE)
!
! zeta
!
      read(81) ((czX_das(i,j),
     &             i=IR_ENTIRE),j=IR_ENTIRE)
      read(81) ((czE_das(i,j),
     &             i=JR_ENTIRE),j=JR_ENTIRE)
!
! cross T/S
!
      read(81) (((corr_ts(i,j,k),
     &          i=IR_ENTIRE),j=JR_ENTIRE),k=1,ndas)

      close(81)

#if defined DAS_COUNT_PSICHI
      count_non_null_psi=0
      count_non_null_chi=0
      do j=IR_ENTIRE
        do i=IR_ENTIRE
          if (abs(cpsiX_das(i,j)) .gt. 0. ) then
            count_non_null_psi=count_non_null_psi+1
          endif
        enddo
      enddo
      do j=IR_ENTIRE
        do i=IR_ENTIRE
          if (abs(cchiX_das(i,j)) .gt. 0. ) then
            count_non_null_chi=count_non_null_chi+1
          endif
        enddo
      enddo
      if ( count_non_null_psi .gt. 0. ) then
         write(*, '(6x,A,i)') 
     &     'YES non null cpsiX_das = ',count_non_null_psi
      else
         write(*, '(6x,A,i)') 
     &     'Ouch all null cpsiX_das = ',count_non_null_chi
      endif 
      if ( count_non_null_chi .gt. 0. ) then
         write(*, '(6x,A,i)') 
     &     'YES non null cchiX_das = ',count_non_null_chi
      else
         write(*, '(6x,A,i)') 
     &     'Ouch all null cchiX_das = ',count_non_null_chi
      endif 
#endif
#if defined DAS_COUNT_PSICHI
      count_non_null_psi=0
      count_non_null_chi=0
      do j=JR_ENTIRE
        do i=JR_ENTIRE
          if (abs(cpsiE_das(i,j)) .gt. 0. ) then
            count_non_null_psi=count_non_null_psi+1
          endif
        enddo
      enddo
      do j=JR_ENTIRE
        do i=JR_ENTIRE
          if (abs(cchiE_das(i,j)) .gt. 0. ) then
            count_non_null_chi=count_non_null_chi+1
          endif
        enddo
      enddo
      if ( count_non_null_psi .gt. 0. ) then
         write(*, '(6x,A,i)') 
     &     'YES non null cpsiE_das = ',count_non_null_psi
      else
         write(*, '(6x,A,i)') 
     &     'Ouch all null cpsiE_das = ',count_non_null_chi
      endif 
      if ( count_non_null_chi .gt. 0. ) then
         write(*, '(6x,A,i)') 
     &     'YES non null cchiE_das = ',count_non_null_chi
      else
         write(*, '(6x,A,i)') 
     &     'Ouch all null cchiE_das = ',count_non_null_chi
      endif 
#endif


!
! vertical correlation, cholesky factorization
!
      open(81, file=file_corr_vert(1:len_corr_vert),
     &            form='unformatted',STATUS='old',IOSTAT=ios)
      if(ios.ne.0) then
        write(*,'(6x,A,(/6x,A))') '??? WARNING: can not find ',
     &           file_corr_vert(1:len_corr_vert)
        return
      endif
!
! temp
!
       read(81) ((ctZ_das(i,j,itemp),
     &           i=1,ndas),j=1,ndas)
!
! salt
!
      read(81) ((ctZ_das(i,j,isalt),
     &           i=1,ndas),j=1,ndas)
!
! psi & chi
!
      read(81) ((cpsiZ_das(i,j),
     &           i=1,ndas),j=1,ndas)
      read(81) ((cchiZ_das(i,j),
     &           i=1,ndas),j=1,ndas)
!
! u & v
!
      read(81) ((cuZ_das(i,j),
     &           i=1,ndas),j=1,ndas)
      read(81) ((cvZ_das(i,j),
     &           i=1,ndas),j=1,ndas)
!
      close(81)

#if defined DAS_COUNT_PSICHI
      count_non_null_psi=0
      count_non_null_chi=0
      do j=ndas,ndas-5,-1
        do i=ndas,ndas-5,-1
          if (abs(cpsiZ_das(i,j)) .gt. 0. ) then
            count_non_null_psi=count_non_null_psi+1
          endif
        enddo
      enddo
      do j=ndas,ndas-5,-1
        do i=ndas,ndas-5,-1
          if (abs(cchiZ_das(i,j)) .gt. 0. ) then
            count_non_null_chi=count_non_null_chi+1
          endif
        enddo
      enddo
      if ( count_non_null_psi .gt. 0. ) then
         write(*, '(6x,A,i)') 
     &     'YES non null cpsiZ_das = ',count_non_null_psi
      else
         write(*, '(6x,A,i)') 
     &     'Ouch all null cpsiZ_das = ',count_non_null_chi
      endif 
      if ( count_non_null_chi .gt. 0. ) then
         write(*, '(6x,A,i)') 
     &     'YES non null cchiZ_das = ',count_non_null_chi
      else
         write(*, '(6x,A,i)') 
     &     'Ouch all null cchiZ_das = ',count_non_null_chi
      endif 
#endif

      return
      end
#else
      subroutine das_get_corr_empty
      return
      end
#endif


      
      

      

      


