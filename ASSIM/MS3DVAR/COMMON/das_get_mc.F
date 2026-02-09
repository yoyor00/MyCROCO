#include "cppdefs.h"
#if defined DAS_MCSST

      subroutine das_get_mc

      implicit none

# include "param.h"
# include "das_param.h"
# include "grid.h"
# include "das_innov.h"
#if defined DAS_SSTMCMSK 
# include "das_buffer.h"
#endif

      integer status, nc_id, v_id, i, j
      logical, parameter :: flag_unit_kelvin=.true.    
#if defined DAS_FDN_MCSST
      integer ifndlev      
#endif
#if defined DAS_SSTMCMSK 
      logical :: flag_mc_mask
#endif

#include "netcdf.inc"
!
! READ mc sst
!
      write(*,'(6x,A,(/8x,A))') 'NOW READING MC SST'

      flag_mc=.true.

      status = NF_OPEN(file_mc,0,nc_id)
      if (status.ne.NF_NOERR) then
        write(*,'(6x,A,(/8x,A))') '??? WARNING: NOT FOUND OR OPENED',
     &          file_mc
        flag_mc=.false.
        goto 999
      endif
!
! sst
!
! replacing 'sst' by analysed_sst SST-L4
!           'adjusted_sea_surface_temperature' SST-L3S
! if SST-L3S WARNING requires to get the cloudy pixel mask !!!
!     status = NF_INQ_VARID(nc_id,'analysed_sst',v_id)
      status = NF_INQ_VARID(nc_id,'sst',v_id)
      if (status.ne.NF_NOERR) then
          print*, '??? WARNING: error in inquiring mc sst'
          print*, '??? trying analysed_sst SST-L4'
          status = NF_INQ_VARID(nc_id,'analysed_sst',v_id)
      endif
      if (status.ne.NF_NOERR) then
          print*, '??? WARNING: error in inquiring mc sst'
          print*, '??? trying adjusted_sea_surface_temperatureanalysed_sst'
          status = 
     & NF_INQ_VARID(nc_id,'adjusted_sea_surface_temperature',v_id)
      endif   
      if (status.ne.NF_NOERR) then
        print*, '??? WARNING: error in inquiring mc sst'
        flag_mc=.false.
        goto 999
      endif
      status = NF_GET_VAR_DOUBLE(nc_id, v_id, anc)
      if (status.ne.NF_NOERR) then
        print*, '??? WARNING: error in reading mc sst'
        flag_mc=.false.
        goto 999
      endif
                
      if (flag_unit_kelvin) then
        print*, '... Kelvin to Celsius conversion'
        !write(*,*) '    sst = ',sst_mc(Lm/2,Mm/2)
        do j=0,Mm+1
          do i=0,Lm+1
#if ! defined DAS_FDN_MCSST
              sst_mc(i,j)=anc(i+1,j+1) - 273.15
#else
              sst_mc(i,j,:)=anc(i+1,j+1) - 273.15
#endif

          enddo
        enddo
      else
        do j=0,Mm+1
          do i=0,Lm+1
#if ! defined DAS_FDN_MCSST
               sst_mc(i,j)=anc(i+1,j+1)
#else
              sst_mc(i,j,:)=anc(i+1,j+1)
#endif
          enddo
        enddo
      endif

#if defined DAS_SSTMCMSK 
      flag_mc_mask=.true.
      if(flag_mc_mask)then
       status = NF_INQ_VARID(nc_id, 'mask_src',v_id)
       if (status.ne.NF_NOERR) then
         print*, '??? WARNING: error in inquiring sst mc mask '
         flag_mc_mask=.false.
       endif
      endif
      if(flag_mc_mask)then
          status = NF_GET_VAR_DOUBLE(nc_id, v_id, anc)
          if (status.ne.NF_NOERR) then
              print*, '??? WARNING: error in reading sst mask'
              flag_mc_mask=.false.
          endif
      endif
      if(flag_mc_mask)then
           do j=0,Mm+1
            do i=0,Lm+1
            ! BLXD shift indices :
            ! anc        array declared as 1:Lm+2,1:Mm+2
            ! das_buff2  array declared as 0:Lm+1,0:Mm+1
              das_buff2(i,j)=anc(i+1,j+1)
            enddo
           enddo
      endif
#endif

!
! mask
!
#if ! defined DAS_FDN_MCSST
      do j=0,Mm+1
        do i=0,Lm+1
          if (sst_mc(i,j) .gt. 0.5 ) then
            mc_mask(i,j)=1.0
          else
            mc_mask(i,j)=0.0
          endif
        enddo
      enddo
#else
      do ifndlev=1,nfndlev
      do j=0,Mm+1
        do i=0,Lm+1
            if (sst_mc(i,j,ifndlev) .gt. 0.5 ) then
              mc_mask(i,j,ifndlev)=1.0
            else
              mc_mask(i,j,ifndlev)=0.0
            endif
        enddo
      enddo
      enddo
#endif

#if defined DAS_SSTMCMSK 
# if ! defined DAS_FDN_MCSST
      if(flag_mc_mask)then
           do j=0,Mm+1
            do i=0,Lm+1
              mc_mask(i,j)=mc_mask(i,j)*das_buff2(i,j)
            enddo
           enddo
      endif
# else
      if(flag_mc_mask)then
        do ifndlev=1,nfndlev
           do j=0,Mm+1
            do i=0,Lm+1
              mc_mask(i,j,ifndlev)=mc_mask(i,j,ifndlev)
     &                             *das_buff2(i,j)
            enddo
           enddo
        enddo
      endif
# endif
#endif

      status = NF_CLOSE(nc_id)

999   continue

      if (flag_mc) then
        write(*, '(6x, A,(/8x,A))') 'OK: MC SST READ IN',
     &         file_mc
      else
        write(*, '(6x, A,(/8x,A))') 'ISSUES: MC SST NOT READ',
     &         file_mc
      endif
      return
      end
#else
      subroutine das_get_mc_empty
      return
      end
#endif /* DAS_MCSST */

