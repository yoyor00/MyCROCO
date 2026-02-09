#include "cppdefs.h"
      subroutine das_get_hfradar
      implicit none
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_innov.h"
      integer i, ios, lenstr,len_file,it, num0
      real hh, spd,uerr,verr,amperr,cr_amp, cr_fit
      real cr_amp6, cr_fit6
!NOTE
!     SSH 00-06, 06-12, 12-18, 18-19 UTC
!     we pick 12-18

!
! initilaizzation
!
      flag_hfradar=.true.
      flag_hfradar6=.true.
      num_hfradar=0
      num_hfradar6=0
      cr_amp=8.5
      cr_amp6=8.5
      cr_fit=8.5
      cr_fit6=8.5

#if defined DAS_HUGEHRADAR
!      Parameter in das_param.h fhuge = 1000.
!      Not used :
!      cr_amphuge =   fhuge * cr_amp
!      fvmaskhuge = - fhuge * 900. 
#endif
!
      do i=1,max_hfradar
        lon_hf(i)=-999.0
        lat_hf(i)=-999.0
        u_hf_raw(i)=-999.0
        v_hf_raw(i)=-999.0
      enddo

!
! NOTE: u/v in cm/s ------
!
!
! HF radar
!
      len_file=lenstr(file_hfradar_uv)
      open(81, file=file_hfradar_uv(1:len_file),form='formatted',
     &         STATUS='old',IOSTAT=ios)
      if(ios.ne.0) then
        flag_hfradar=.false.
        write(*,'(6x,A,(/8x,A))') '??? WARNING: NOT FOUND',
     &          file_hfradar_uv
        goto 10
      endif
!skip the 11 head lines
      flag_hfradar=.true.
#if defined DAS_HFHEADER
      do it=1,2
        read(81,*)
      enddo
#endif
      do it=1,max_hfradar
!        read(81,800,end=100)
         read(81,*,end=100)
     &        lat_hf(it),lon_hf(it),u_hf_raw(it),v_hf_raw(it)
     &       ,spd, uerr,verr
        amperr = sqrt(uerr*uerr + verr*verr)
        amperr_hf(it)=amperr
        fitdif_hf(it)=0.0
!        if (amperr .gt. cr_amp .or. fitdif .gt. cr_fit ) then
        if (amperr .gt. cr_amp ) then
          u_hf_raw(it)=-999.0
          v_hf_raw(it)=-999.0
        endif
      enddo
100   continue
      num_hfradar=it-1
      write(*,'(6x,A,I6, (/8x,A))')
     &   'OK: HFRADAR UV READ IN TOTAL: ', num_hfradar, file_hfradar_uv
10    continue
      close(81)
800   format(f9.4,6f10.4)

#if !defined DAS_HUGEHRADAR
      num0=0
      do it=1,num_hfradar
        if (u_hf_raw(it).gt. -900.0 .and. v_hf_raw(it) .gt. -900.0) then
          num0=num0+1
          lon_hf(num0)=lon_hf(it)
          lat_hf(num0)=lat_hf(it)
          amperr_hf(num0)=amperr_hf(it)
          fitdif_hf(num0)=fitdif_hf(it)
# if !defined DAS_SYNTHRADAR
          u_hf_raw(num0)=u_hf_raw(it)*0.01
          v_hf_raw(num0)=v_hf_raw(it)*0.01
# else
          u_hf_raw(num0)=u_hf_raw(it)
          v_hf_raw(num0)=v_hf_raw(it)
# endif

        endif
      enddo
      num_hfradar=num0
#else
      num0=0
      do it=1,num_hfradar
        if (u_hf_raw(it).gt. -900.0 .and. v_hf_raw(it) .gt. -900.0) then
          num0=num0+1
          lon_hf(num0)=lon_hf(it)
          lat_hf(num0)=lat_hf(it)
          amperr_hf(num0)= fhuge * amperr_hf(it)
          fitdif_hf(num0)= fhuge * fitdif_hf(it)
! BLXD cm/s to m/s conversion removed for "huge" test 
          u_hf_raw(num0) = fhuge * u_hf_raw(it)
          v_hf_raw(num0) = fhuge * v_hf_raw(it)
        endif
      enddo
#endif

      num_hfradar=num0
      write(*,'(6x,A,I6, (/8x,A))')
     &   'OK: HFRADAR UV SCREENED IN TOTAL: ', num_hfradar

      write(*,801) (lon_hf(it),lat_hf(it),u_hf_raw(it),v_hf_raw(it)
     &          , amperr_hf(it),fitdif_hf(it),it=1,num_hfradar)
801   format(2(1x,f10.4),1x,2f8.5,2f8.5/)

!
! 6KM HG radar
!
      do i=1,max_hfradar6
        lon_hf6(i)=-999.0
        lat_hf6(i)=-999.0
        u_hf_raw6(i)=-999.0
        v_hf_raw6(i)=-999.0
      enddo

!
! NOTE: u/v in cm/s ------
!
!
! HF radar
!
      len_file=lenstr(file_hfradar_uv6)
      open(81, file=file_hfradar_uv6(1:len_file),form='formatted',
     &         STATUS='old',IOSTAT=ios)
      if(ios.ne.0) then
        flag_hfradar=.false.
        write(*,'(6x,A,(/8x,A))') '??? WARNING: NOT FOUND',
     &          file_hfradar_uv6(1:len_file)
        goto 20
      endif
!skip 2 head lines
      flag_hfradar=.true.
#if defined DAS_HFHEADER
      do it=1,2
        read(81,*)
      enddo
#endif
      do it=1,max_hfradar6
        read(81,700,end=200)
     &        lat_hf6(it),lon_hf6(it),u_hf_raw6(it),v_hf_raw6(it)
     &       ,spd, uerr,verr
        amperr = sqrt(uerr*uerr + verr*verr)
        amperr_hf6(it)=amperr
        fitdif_hf6(it)=0.0
!        if (amperr .gt. cr_amp6 .or. fitdif .gt. cr_fit6 ) then
        if (amperr .gt. cr_amp6 ) then
          u_hf_raw6(it)=-999.0
          v_hf_raw6(it)=-999.0
        endif
      enddo
200   continue
      num_hfradar6=it-1
      write(*,'(6x,A,I6, (/8x,A))')
     &   'OK: HFRADAR 6KM READ IN TOTAL: ', num_hfradar6,
     &    file_hfradar_uv6(1:len_file)

20    continue
      close(81)
700   format(f9.4,6f10.4)

      num0=0
      do it=1,num_hfradar6
        if (u_hf_raw6(it).gt. -900.0 .and. v_hf_raw6(it) .gt. -900.0) then
          num0=num0+1
          lon_hf6(num0)=lon_hf6(it)
          lat_hf6(num0)=lat_hf6(it)
          amperr_hf6(num0)=amperr_hf6(it)
          fitdif_hf6(num0)=fitdif_hf6(it)
          u_hf_raw6(num0)=u_hf_raw6(it)*0.01
          v_hf_raw6(num0)=v_hf_raw6(it)*0.01
        endif
      enddo
      num_hfradar6=num0


!      write(*,802) (lon_hf6(it),lat_hf6(it),u_hf_raw6(it),v_hf_raw6(it)
!     &          , amperr_hf6(it),fitdif_hf6(it),it=1,100)
802   format(2f10.4,2f8.2,2f8.2)

      return
      end
