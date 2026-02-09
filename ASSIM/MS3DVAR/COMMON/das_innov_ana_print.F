#include "cppdefs.h"

      subroutine das_innov_ana_print
      implicit none
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_innov.h"
# include "das_ocean9.h"
!
      integer jj,k,lenstr,len_file
      real ht_ana(max_prf), hs_ana(max_prf)
!
      write(81,11) 'INNOVATIONS FOR MOORING, NUM OF PROFILES: '
     &             ,prf_num_moor
      do jj=1,prf_num_moor
        write(81,12) 'XI & ETA', Imoor(jj), Jmoor(jj)
     &      ,lonr(Imoor(jj), Jmoor(jj)),latr(Imoor(jj), Jmoor(jj))
        write(81,9) 'DEPTH','TEMP','SALT','MASK','DAS T','DAS S'
     &             ,'FCST T','FCST S'
        do k=1,max_prf
          ht_ana(k)=tmask_moor(jj,k)*(t_moor(jj,k)
     &             -t_s(Imoor(jj),Jmoor(jj),ndas-k+1,itemp))
          hs_ana(k)=smask_moor(jj,k)*(s_moor(jj,k)
     &             -t_s(Imoor(jj),Jmoor(jj),ndas-k+1,isalt))
        enddo
        write(81,10)
     &        (z_das(ndas-k+1), t_moor(jj,k),
     &        s_moor(jj,k), smask_moor(jj,k), ht_ana(k),hs_ana(k),
     &         t_das9(Imoor(jj),Jmoor(jj),ndas-k+1,itemp),
     &         t_das9(Imoor(jj),Jmoor(jj),ndas-k+1,isalt),
     &                                     k=1,max_prf)
      enddo
      close(81)
!
      close(81)

10    format (8f10.5)
11    format (A38,I6)
12    format (A8,2I6,2f11.5)
9     format (8A10)
      return
      end
      

