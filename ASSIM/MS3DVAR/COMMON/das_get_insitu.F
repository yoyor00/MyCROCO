#include "cppdefs.h"
      subroutine das_get_insitu
      implicit none
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_innov.h"
      integer i,j,k,ios,ios1,ios2,lenstr,len_file,it,num,itp
      real hh
!NOTE
!     SSH 00-06, 06-12, 12-18, 18-19 UTC
!     we pick 12-18

!
! initilaizzation
!
      flag_sio=.true.
      flag_whoi=.true.
      flag_ptsur=.true.
      flag_martn=.true.
      flag_map=.true.
      flag_prof=.true.
      flag_flight=.true.
      flag_ship=.true.
      prf_num_sio=0
      prf_num_whoi=0
      prf_num_ptsur=0
      prf_num_martn=0
      prf_num_map=0
      num_flight=0
!
! sio
!
      do i=1,max_sio 
        lon_sio(i)=-9999.0
        lat_sio(i)=-9999.0
      enddo
      do j=1,max_prf
        do i=1,max_sio
          t_sio_raw(i,j)=-9999.0
          s_sio_raw(i,j)=-9999.0
        enddo
      enddo
!
! whoi
!
      do i=1,max_whoi
        lon_whoi(i)=-9999.0
        lat_whoi(i)=-9999.0
      enddo
      do j=1,max_prf
        do i=1,max_whoi
          t_whoi_raw(i,j)=-9999.0
          s_whoi_raw(i,j)=-9999.0
        enddo
      enddo
!     
! mapped
!
      do i=1,max_map
        lon_map(i)=-9999.0
        lat_map(i)=-9999.0
      enddo
      do j=1,max_prf
        do i=1,max_map
          t_map_raw(i,j)=-9999.0
          s_map_raw(i,j)=-9999.0
        enddo
      enddo

!
!
! ptsur
!
      do i=1,max_ptsur
        lon_ptsur(i)=-9999.0
        lat_ptsur(i)=-9999.0
      enddo
      do j=1,max_prf
        do i=1,max_ptsur
          t_ptsur_raw(i,j)=-9999.0
          s_ptsur_raw(i,j)=-9999.0
        enddo
      enddo
!
!martn
!
      do i=1,max_martn
        lon_martn(i)=-9999.0
        lat_martn(i)=-9999.0
        Imartn(i)=-9999
        Jmartn(i)=-9999
      enddo
      do j=1,max_prf
        do i=1,max_martn
          mask_martn(i,j)=0.0
          t_martn(i,j)=-9999.0
          s_martn(i,j)=-9999.0
        enddo
      enddo
!
      do i=1,max_flight
        lon_fl(i)=-9999.0
        lat_fl(i)=-9999.0
        sst_fl_raw(i)=-9999.0
      enddo
!
      do i=1,max_flight
        lon_fl(i)=-9999.0
        lat_fl(i)=-9999.0
        sst_fl_raw(i)=-9999.0
      enddo

!
! sio_glider
!
      len_file=lenstr(file_sio_glider)
      open(81, file=file_sio_glider(1:len_file),form='formatted',
     &         STATUS='old',IOSTAT=ios)
      if(ios.ne.0) then
        flag_sio=.false.
        write(*,'(6x,A,(/8x,A))') '??? WARNING: SIO NOT FOUND ',
     &          file_sio_glider
        goto 10
      endif
!
      do it=1,max_sio
        do k=1,max_prf
          read(81,900,end=100)
     &     lon_sio(it),lat_sio(it),hh,
     &     t_sio_raw(it,k),s_sio_raw(it,k)      
        enddo
      enddo
100   continue
      prf_num_sio=it-1
      write(*,'(6x,A,I5, (/8x,A))')  'OK: SIO GLIDER READ IN',
     &         prf_num_sio,file_sio_glider
10    continue
      close(81)
!
! whoi_glider
!
      len_file=lenstr(file_whoi_glider)
      open(81, file=file_whoi_glider(1:len_file),form='formatted',
     &         STATUS='old',IOSTAT=ios)
      if(ios.ne.0) then
        flag_whoi=.false.
        write(*,'(6x,A,(/8x,A))') '??? WARNING: WHOI NOT FOUND ',
     &          file_whoi_glider
        goto 20
      endif
!
      do it=1,max_whoi
        do k=1,max_prf
          read(81,900,end=200)
     &     lon_whoi(it),lat_whoi(it),hh,
     &     t_whoi_raw(it,k),s_whoi_raw(it,k)
        enddo
      enddo
200   continue
      prf_num_whoi=it-1
      write(*,'(6x,A,I5, (/8x,A))')  'OK: WHOI GLIDER READ IN',
     &         prf_num_whoi,file_whoi_glider
20    continue
      close(81)

!                                                                                                           
! maping ctd                                                                                               
!                                                                                                           
      len_file=lenstr(file_map_ctd)                              
      open(81, file=file_map_ctd(1:len_file),form='formatted',   
     &         STATUS='old',IOSTAT=ios)                           
      if(ios.ne.0) then                                           
        flag_map=.false.                                         
        write(*,'(6x,A,(/8x,A))') '??? WARNING: NOT FOUND ',      
     &         file_map_ctd                                      
        goto 150                                                  
      endif                                                       
      do it=1,max_map                                            
        do k=1,max_prf                                            
          read(81,900,end=1500)                                   
     &     lon_map(it),lat_map(it),hh,
     &     t_map_raw(it,k),s_map_raw(it,k)
        enddo                                                    
      enddo                                                      
1500  continue                                                   
      prf_num_map=it-1                                          
      write(*,'(6x,A,I5,(/8x,A))')                               
     &         'OK: MAPPED TS READ IN',                       
     &          prf_num_map,file_map_ctd                      
150    continue                                                 
      close(81)
!    
! profile ctd
!     
      len_file=lenstr(file_prof_ctd)
      open(81, file=file_prof_ctd(1:len_file),form='formatted',
     &         STATUS='old',IOSTAT=ios)
      if(ios.ne.0) then
        flag_prof=.false.
        write(*,'(6x,A,(/8x,A))') '??? WARNING: NOT FOUND ',
     &         file_prof_ctd
        goto 160
      endif
      do it=1,max_map-prf_num_map
          itp=prf_num_map+it
        do k=1,max_prf
          read(81,900,end=1600)
     &     lon_map(itp),lat_map(itp),hh,t_map_raw(itp,k),
     &     s_map_raw(itp,k)
        enddo
      enddo
1600  continue
      prf_num_map=itp-1
      write(*,'(6x,A,(/8x,A))')
     &         'OK: PROFILE CDT READ IN',file_prof_ctd
160    continue
      close(81)

!
! ptsur ctd
!
      len_file=lenstr(file_ptsur_ctd)
      open(81, file=file_ptsur_ctd(1:len_file),form='formatted',
     &         STATUS='old',IOSTAT=ios)
      if(ios.ne.0) then
        flag_ptsur=.false.
        write(*,'(6x,A,(/8x,A))') '??? WARNING: PT NOT FOUND ',
     &         file_ptsur_ctd
        goto 40
      endif
      do it=1,max_ptsur
        do k=1,max_prf
          read(81,900,end=400)
     &     lon_ptsur(it),lat_ptsur(it),hh,
     &     t_ptsur_raw(it,k),s_ptsur_raw(it,k)
        enddo
      enddo
400   continue
      prf_num_ptsur=it-1
      write(*,'(6x,A,I5,(/8x,A))')
     &         'OK: PTSUR CDT READ IN',
     &          prf_num_ptsur,file_ptsur_ctd
40    continue
      close(81)
!
! martn ctd
!
      len_file=lenstr(file_martn_ctd)
      open(81, file=file_martn_ctd(1:len_file),form='formatted',
     &         STATUS='old',IOSTAT=ios)
      if(ios.ne.0) then
        flag_martn=.false.
        write(*,'(6x,A,(/8x,A))') '??? WARNING: NOT FOUND ',
     &         file_martn_ctd
        goto 50
      endif
      do it=1,max_martn
        do k=1,max_prf
          read(81,900,end=500)
     &     lon_martn(it),lat_martn(it),hh,t_martn(it,k),
     &     s_martn(it,k)
        enddo
      enddo
500   continue
      prf_num_martn=it-1
      write(*,'(6x,A,(/8x,A))')
     &         'OK: MARTIN CDT READ IN',file_martn_ctd
50    continue
      close(81)
!
      do it=1,num_flight
        if ( lon_fl(it) .lt. 0.0)
     &       lon_fl(it)=360.0+lon_fl(it)
      enddo
!
! ship
!
      len_file=lenstr(file_ship_sst)
      open(81, file=file_ship_sst(1:len_file),form='formatted',
     &         STATUS='old',IOSTAT=ios)
      if(ios.ne.0) then
        flag_ship=.false.
        write(*,'(6x,A,(/8x,A))') '??? WARNING: NOT FOUND ',
     &          file_ship_sst
        goto 130
      endif
      do it=1,max_flight
        read(81,800,end=1300)
     &   lon_fl(it),lat_fl(it),sst_fl_raw(it)
      enddo
1300  continue
      num_flight=it-1
      write(*,'(6x,A,(/8x,A))')
     &         'OK: SHIP SST READ IN',file_ship_sst
130   continue
      close(81)
!
! flight
!
      len_file=lenstr(file_flight_sst)
      open(81, file=file_flight_sst(1:len_file),form='formatted',
     &         STATUS='old',IOSTAT=ios)
      if(ios.ne.0) then
        flag_flight=.false.
        write(*,'(6x,A,(/8x,A))') '??? WARNING: NOT FOUND ',
     &          file_flight_sst
        goto 30
      endif
      do it=1,max_flight-num_flight
        itp=num_flight+it
        read(81,800,end=300)
     &   lon_fl(itp),lat_fl(itp),sst_fl_raw(itp)
      enddo
300   continue
      num_flight=itp-1
      write(*,'(6x,A,(/8x,A))')
     &         'OK: FLIGHT SST READ IN',file_flight_sst
30    continue
      close(81)
!
800   format(3f10.4)
900   format(5f11.4)
      return
      end
