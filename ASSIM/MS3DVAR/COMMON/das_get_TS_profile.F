#include "cppdefs.h"
      subroutine das_get_TS_obs
      implicit none
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_innov.h"
      integer i,j,k,ios,lenstr,len_file,it,num
!
!NOTE
!
!     T/S vertical profiles
!     one data with one DA layer ( not at DA levels as before)
!     Each profiles thus have NDAS-1 data points
!     missing valuses as -999.0
!
! initilaizzation
!
      flag_TS_obs=.true.
      prf_num_TS_obs=0
!
!      max_TS_obs=1200;     ! in das_param.h
!
! T/S Obs profles
!
      do i=1,max_TS_obs
        lon_obs(i)=-9999.0
        lat_obs(i)=-9999.0
      enddo
      do k=1,NDAS-1
        do i=1,max_TS_obs
          t_obs_raw(i,k)=-9999.0
          s_obs_raw(i,k)=-9999.0
        enddo
      enddo

!
! profiles
!
      len_file=lenstr(file_prf_TS_obs)
      open(81, file=file_prf_TS_obs(1:len_file),form='formatted',
     &         STATUS='old',IOSTAT=ios)
      if(ios.ne.0) then
        flag_TS_obs=.false.
        write(*,'(6x,A,(/8x,A))') '??? WARNING: TS Profile NOT FOUND ',
     &          file_prf_TS_obs
        goto 10
      endif
!
      do it=1,max_TS_obs
        do k=1,NDAS-1
          read(81,900,end=100)
     &     I_obs(it), J_obs(it), K_obs(it),
     &     t_obs_raw(it,k),dep_tobs_raw(it,k),
     &     s_obs_raw(it,k),dep_sobs_raw(it,k),
     &     time_obs(it),err_tobs_raw(it,k),err_sobs_raw(it,k),
     &     lon_obs(it),lat_obs(it),inst_obs(it)
        enddo
      enddo
100   continue
      prf_num_TS_obs=it-1
      write(*,'(6x,A,I8,(/8x,A))')  'OK: TS Profile OBS  READ IN',
     &         prf_num_TS_obs, file_prf_TS_obs
10    continue
      close(81)
!
      do it=1,prf_num_TS_obs
        do k=1,NDAS-1
           write(*,900)
     &     I_obs(it), J_obs(it), K_obs(it),
     &     t_obs_raw(it,k),dep_tobs_raw(it,k),
     &     s_obs_raw(it,k),dep_sobs_raw(it,k),
     &     time_obs(it),err_tobs_raw(it,k),err_sobs_raw(it,k),
     &     lon_obs(it),lat_obs(it),inst_obs(it)
        enddo
      enddo
900   format(3I8,9F11.4,I8)
!901   format(3I8,9F11.4,I8)

      return
      end
