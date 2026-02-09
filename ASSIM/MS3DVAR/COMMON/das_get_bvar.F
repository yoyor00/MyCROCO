#include "cppdefs.h"
#if !defined DAS_ANA_BVAR

      subroutine das_get_bvar  ! from binary file
      implicit none
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
!
      integer i,j,k,ios,lenstr,len_bvar
!
! compute the entire model domain
! but the lower triangular matrices
!
      len_bvar=lenstr(file_bvar)
!
! depth dependent error variance
!
      open(81, file=file_bvar(1:len_bvar),
     &            form='unformatted',STATUS='old',IOSTAT=ios)
      if(ios.ne.0) then
        write(*,'(6x,A,(/6x,A))') '??? WARNING: can not find ',
     &           file_bvar(1:len_bvar)
        return
      endif
!
! t
!
      read(81)
     &      (((bt_das(i,j,k,itemp),i=0,Lm+1),j=0,Mm+1),k=1,ndas)
!
! s
!
      read(81) 
     &      (((bt_das(i,j,k,isalt),i=0,Lm+1),j=0,Mm+1),k=1,ndas)
!
! psi
!
      read(81)
     &      (((bpsi_das(i,j,k),i=1,Lm+1),j=1,Mm+1),k=1,ndas)
!
! chi
!
      read(81)
     &      (((bchi_das(i,j,k),i=0,Lm+1),j=0,Mm+1),k=1,ndas)
!
! u
!
      read(81) 
     &      (((bu_das(i,j,k),i=1,Lm+1),j=0,Mm+1),k=1,ndas)
!
! v
!
      read(81)
     &      (((bv_das(i,j,k),i=0,Lm+1),j=1,Mm+1),k=1,ndas)
!
! zeta
!
      read(81) 
     &      ((bz_das(i,j),i=0,Lm+1),j=0,Mm+1)
!
      close(81)
!
      return
      end
#else
      subroutine das_get_bvar_empty
      return
      end
#endif
