#include "cppdefs.h"

      subroutine das_innov_print
      implicit none
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_innov.h"
!
      integer jj,k,lenstr,len_file
!
! print sio data
!

      len_file=lenstr(file_innov_sio)
      open (81,file=file_innov_sio(1:len_file),form='formatted',
     &              STATUS='unknown')
      write(81,11) 'INNOVATIONS FOR SIO, NUM OF PROFILES: '
     &             ,prf_num_sio
      do jj=1,prf_num_sio
        write(81,12) 'XI & ETA', Isio(jj), Jsio(jj)
     &      ,lonr(Isio(jj), Jsio(jj)),latr(Isio(jj), Jsio(jj))               
        write(81,9) 'DEPTH','TEMP','SALT','MASK','FCST T','FCST S'
        write(81,10) 
     &        (z_das(ndas-k+1), t_sio(jj,k),
     &         s_sio(jj,k), mask_sio(jj,k), 
     &         t_das(Isio(jj),Jsio(jj),ndas-k+1,itemp),
     &         t_das(Isio(jj),Jsio(jj),ndas-k+1,isalt),
     &                                     k=1,max_prf)
      enddo
!
!mooring
!
      write(81,11) 'INNOVATIONS FOR MOORINGS, NUM OF PROFILES: '
     &             ,prf_num_moor
      do jj=1,prf_num_moor
        write(81,12) 'XI & ETA', Imoor(jj), Jmoor(jj)
     &      ,lonr(Imoor(jj), Jmoor(jj)),latr(Imoor(jj), Jmoor(jj))
        write(81,9) 'DEPTH','TEMP','SALT','MASK','FCST T','FCST S'
        write(81,10)
     &        (z_das(ndas-k+1), t_moor(jj,k),
     &         s_moor(jj,k), smask_moor(jj,k),
     &         t_das(Imoor(jj),Jmoor(jj),ndas-k+1,itemp),
     &         t_das(Imoor(jj),Jmoor(jj),ndas-k+1,isalt),
     &                                     k=1,max_prf)
      enddo

      close(81)
!
! print whoi data
!
      len_file=lenstr(file_innov_whoi)
      open (81,file=file_innov_whoi(1:len_file),form='formatted',
     &              STATUS='unknown')
      write(81,*) 'INNOVATIONS FOR WHOI, NUM OF PROFILES:'
     &             ,prf_num_whoi
      do jj=1,prf_num_whoi
        write(81,12) 'XI & ETA', Iwhoi(jj), Jwhoi(jj)
     &      ,lonr(Iwhoi(jj), Jwhoi(jj)),latr(Iwhoi(jj), Jwhoi(jj))               
        write(81,9) 'DEPTH','TEMP','SALT','MASK','FCST T','FCST S'
        write(81,10)
     &        (z_das(ndas-k+1),t_whoi(jj,k),
     &         s_whoi(jj,k),mask_whoi(jj,k), 
     &         t_das(Iwhoi(jj),Jwhoi(jj),ndas-k+1,itemp),
     &         t_das(Iwhoi(jj),Jwhoi(jj),ndas-k+1,isalt),
     &                                     k=1,max_prf)

      enddo
      close(81)

10    format (6f11.5)
11    format (A38,I6)
12    format (A8,2I6,2f11.5)
9     format (6A11)
      return
      end
      

