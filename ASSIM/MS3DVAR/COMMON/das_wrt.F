#include "cppdefs.h"

      subroutine das_wrt

      implicit none
#include "netcdf.inc"
#include "param.h"
#include "das_param.h"
#include "das_ocean.h"
#include "das_ocean_diag.h"
#include "grid.h"
#include "ncscrum.h"
#include "das_ocean9.h"

      character das_start_day*20, vnts*4
      integer i,j,k,lstr, lenstr, das_nf_fwrite, indx,varid,itrc
      integer inclstr
      character*180 incname

* error status return
      integer  iret
* netCDF id
      integer  ncid
* dimension ids
      integer  xi_rho_dim
      integer  xi_u_dim
      integer  eta_rho_dim
      integer  eta_v_dim
      integer  s_rho_dim
      integer  time_dim
* dimension lengths
      integer  xi_rho_len
      integer  xi_u_len
      integer  eta_rho_len
      integer  eta_v_len
      integer  s_rho_len
      integer  time_len
      parameter (time_len = NF_UNLIMITED)
* variable ids
      integer  depth_id
      integer  u_da_id
      integer  v_da_id
#if defined DAS_UVAGEO_DIAG
      integer  u_ap_da_id
      integer  v_ap_da_id
      integer  u_ac_da_id
      integer  v_ac_da_id
#endif
#if defined DAS_BKGZETA
      integer  zetabssm_da_id
      integer  zetadas9_da_id
      integer  rhodas9_da_id
#endif
      integer  t_da_id
      integer  s_da_id
      integer  zeta_da_id
      integer  zeta_h_da_id
      integer  mask_da_id
      integer  lat_da_id
      integer  lon_da_id
* rank (number of dimensions) for each variable
      integer  depth_rank
      integer  u_da_rank
      integer  v_da_rank
      integer  t_da_rank
      integer  s_da_rank
      integer  mask_da_rank
      integer  lat_da_rank
      integer  lon_da_rank
      integer  zeta_da_rank
      parameter (depth_rank = 1)
      parameter (u_da_rank = 3)
      parameter (v_da_rank = 3)
      parameter (t_da_rank = 3)
      parameter (s_da_rank = 3)
      parameter (mask_da_rank = 3)
      parameter (lat_da_rank = 2)
      parameter (lon_da_rank = 2)
      parameter (zeta_da_rank = 2)
* variable shapes
      integer  depth_dims(depth_rank)
      integer  u_da_dims(u_da_rank)
      integer  v_da_dims(v_da_rank)
      integer  t_da_dims(t_da_rank)
      integer  s_da_dims(s_da_rank)
      integer  mask_da_dims(mask_da_rank)
      integer  lat_da_dims(lat_da_rank)
      integer  lon_da_dims(lon_da_rank)
      integer  zeta_da_dims(zeta_da_rank)

      real doubleval(1)

      integer nc_start(3) 
      integer nc_count(3) 
      integer nc_start2(2) 
      integer nc_count2(2) 


 
      lstr=lenstr(rstname)
!     remove .nc from rstname
      inclstr=lstr-3
      incname=rstname(1:inclstr)//'_inc.nc'
      inclstr=lenstr(incname)
* enter define mode
!     iret = nf_create(rstname(1:lstr), NF_CLOBBER, ncid)
      iret = nf_create(incname(1:inclstr), NF_CLOBBER, ncid)
      print *, '... Before creating increment of analysis file'
      print *,incname(1:inclstr)
      call check_err(iret)
* define dimensions

      iret = nf_def_dim(ncid, 'xi_rho', xi_rho, xi_rho_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'xi_u', xi_u, xi_u_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'eta_rho', eta_rho, eta_rho_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'eta_v', eta_v, eta_v_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'ndas', ndas, s_rho_dim)
      call check_err(iret)
* define variables

!      depth_dims(1) = s_rho_dim
!      iret = nf_def_var(ncid, 'depth', NF_DOUBLE, 1, depth_dims
!     1, depth_id)
!      call check_err(iret)

      u_da_dims(3) =  s_rho_dim
      u_da_dims(2) =  eta_rho_dim
      u_da_dims(1) = xi_u_dim
      iret = nf_def_var(ncid, 'u_da', NF_DOUBLE, u_da_rank, 
     &       u_da_dims, u_da_id)
      call check_err(iret)

      v_da_dims(3) =  s_rho_dim
      v_da_dims(2) =  eta_v_dim
      v_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'v_da', NF_DOUBLE, v_da_rank,
     &       v_da_dims, v_da_id)
      call check_err(iret)

#if defined DAS_UVAGEO_DIAG
      write(*,*) '... Defining ageostrophic velocities'
      u_da_dims(3) =  s_rho_dim
      u_da_dims(2) =  eta_rho_dim
      u_da_dims(1) = xi_u_dim
      iret = nf_def_var(ncid, 'u_ap_da', NF_DOUBLE, u_da_rank, 
     &       u_da_dims, u_ap_da_id)
      call check_err(iret)
      iret = nf_def_var(ncid, 'u_ac_da', NF_DOUBLE, u_da_rank, 
     &       u_da_dims, u_ac_da_id)
      call check_err(iret)

      v_da_dims(3) =  s_rho_dim
      v_da_dims(2) =  eta_v_dim
      v_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'v_ap_da', NF_DOUBLE, v_da_rank,
     &       v_da_dims, v_ap_da_id)
      call check_err(iret)
      iret = nf_def_var(ncid, 'v_ac_da', NF_DOUBLE, v_da_rank,
     &       v_da_dims, v_ac_da_id)
      call check_err(iret)
      write(*,*) '... End defining ageostrophic velocities'
#endif

      t_da_dims(3) =  s_rho_dim
      t_da_dims(2) =  eta_rho_dim
      t_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 't_da', NF_DOUBLE, t_da_rank,
     &       t_da_dims, t_da_id)
      call check_err(iret)

      s_da_dims(3) =  s_rho_dim
      s_da_dims(2) =  eta_rho_dim
      s_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 's_da', NF_DOUBLE, s_da_rank,
     &       s_da_dims, s_da_id)
      call check_err(iret)

      mask_da_dims(3) =  s_rho_dim                                                                                                 
      mask_da_dims(2) =  eta_rho_dim                                                                                               
      mask_da_dims(1) = xi_rho_dim                                                                                                 
      iret = nf_def_var(ncid, 'mask_da', NF_DOUBLE, mask_da_rank,                                                                     
     &       mask_da_dims, mask_da_id)                                                                                                
      call check_err(iret)

#if defined DAS_BKGZETA
      write(*,*) '... Defining DA modified bkg SSH'
      zeta_da_dims(2) = eta_rho_dim
      zeta_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'zetab_h_da', NF_DOUBLE, zeta_da_rank,
     &       zeta_da_dims, zetabssm_da_id)
      call check_err(iret)
      write(*,*) '... Defining original bkg SSH'
      iret = nf_def_var(ncid, 'zetadas9_da', NF_DOUBLE, zeta_da_rank,
     &       zeta_da_dims, zetadas9_da_id)
      call check_err(iret)
      write(*,*) '... Defining DA bkg density'
      t_da_dims(3) =  s_rho_dim
      t_da_dims(2) =  eta_rho_dim
      t_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'rhodas9_da', NF_DOUBLE, t_da_rank,
     &       t_da_dims, rhodas9_da_id)
      call check_err(iret)
#endif
      write(*,*) '... Defining DA total SSH increment'
      zeta_da_dims(2) =  eta_rho_dim
      zeta_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'zeta_da', NF_DOUBLE, zeta_da_rank,
     &       zeta_da_dims, zeta_da_id)
      call check_err(iret)

      write(*,*) '... Defining DA steric SSH increment'
      zeta_da_dims(2) =  eta_rho_dim
      zeta_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'zeta_h_da', NF_DOUBLE, zeta_da_rank,
     &       zeta_da_dims, zeta_h_da_id)
      call check_err(iret)

      lat_da_dims(2) =  eta_rho_dim
      lat_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'lat_da', NF_DOUBLE, lat_da_rank,
     &       lat_da_dims, lat_da_id)
      call check_err(iret)

      lon_da_dims(2) =  eta_rho_dim
      lon_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'lon_da', NF_DOUBLE, lon_da_rank,
     &       lon_da_dims, lon_da_id)
      call check_err(iret)

* assign attributes

      iret = nf_put_att_text(ncid, u_da_id, 'long_name', 1, 'u') 
      call check_err(iret)
!      doubleval(1) = -9999
      doubleval(1) = 0.0
!      iret = nf_put_att_double(ncid, u_da_id, 'missing_value', NF_DOUBLE
!     1, 1, doubleval)
!     call check_err(iret)
      iret = nf_put_att_text(ncid, v_da_id, 'long_name', 1, 'v') 
      call check_err(iret)
!      iret = nf_put_att_double(ncid, v_da_id, 'missing_value', NF_DOUBLE
!     1, 1, doubleval)
#if defined DAS_UVAGEO_DIAG
      write(*,*) '... Writing DA ageo vel increment attributes'
      iret = nf_put_att_text(ncid, u_ap_da_id, 'long_name', 1, 'u_ap') 
      call check_err(iret)
      iret = nf_put_att_text(ncid, v_ap_da_id, 'long_name', 1, 'v_ap') 
      call check_err(iret)
      iret = nf_put_att_text(ncid, u_ac_da_id, 'long_name', 1, 'u_ac') 
      call check_err(iret)
      iret = nf_put_att_text(ncid, v_ac_da_id, 'long_name', 1, 'v_ac') 
      call check_err(iret)
      write(*,*) '... End writing DA ageo vel increment attributes'
#endif

      iret = nf_put_att_text(ncid, t_da_id, 'long_name', 4, 'temp') 
      call check_err(iret)
!      iret = nf_put_att_double(ncid, t_da_id, 'missing_value', NF_DOUBLE
!     1, 1, doubleval)
      iret = nf_put_att_text(ncid, s_da_id, 'long_name', 4, 'salt') 
      call check_err(iret)
!      iret = nf_put_att_double(ncid, s_da_id, 'missing_value', NF_DOUBLE
!     1, 1, doubleval)
      iret = nf_put_att_text(ncid, mask_da_id, 'long_name', 4, 'salt') 
      call check_err(iret)
!      iret = nf_put_att_double(ncid, mask_da_id, 'missing_value', NF_DOUBLE
!     1, 1, doubleval)
      iret = nf_put_att_text(ncid, zeta_da_id, 'long_name', 4, 'zeta') 
      call check_err(iret)
!      iret = nf_put_att_double(ncid, zeta_da_id, 'missing_value', NF_DOUBLE
!     1, 1, doubleval)
      write(*,*) '... Writing DA steric SSH increment attributes'
      iret = nf_put_att_text(ncid, zeta_h_da_id, 'long_name', 4, 'zeta_h') 
      call check_err(iret)
      write(*,*) '... End writing DA steric SSH increment attributes'
      iret = nf_put_att_text(ncid, lat_da_id, 'long_name', 8, 'latitude') 
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_da_id, 'long_name', 9, 'longitude') 
      call check_err(iret)

      iret = nf_put_att_text(ncid, NF_GLOBAL, 'type', 7, 'ROMS LR')
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'data', 16, 'ROMS vert interp')
      call check_err(iret)

* leave define mode
      iret = nf_enddef(ncid)
      call check_err(iret)
       
      iret = nf_close(ncid)
      call check_err(iret)


!
! write
      indx=1              ! record
      iret = NF_OPEN(incname(1:inclstr),NF_WRITE,ncid) 

#if defined DAS_BKGZETA
      write(*,*) '... Writing DA modified bkg SSH'
      iret=nf_inq_varid (ncid, 'zetab_h_da', varid)                     
      if (iret .eq. nf_noerr) then                                  
        iret=das_nf_fwrite (zetab_das9(START_2D_ARRAY), ncid, varid,   
     &                                         indx, r2dvar)      
        if (iret .ne. nf_noerr) then                             
          write(stdout,2) 'zetab_h_da',indx,incname(1:inclstr)        
          goto 99                                         !--> ERROR                                                             
        endif                                                        
      else                                                          
        write(stdout,1) 'zetab_h_da', incname(1:inclstr)
        goto 99                                           !--> ERROR                                                             
      endif                                                           

      write(*,*) '... Writing DA original bkg SSH'
      iret=nf_inq_varid (ncid, 'zetadas9_da', varid)                     
      if (iret .eq. nf_noerr) then                                  
        iret=das_nf_fwrite (zeta_das9(START_2D_ARRAY), ncid, varid,   
     &                                         indx, r2dvar)      
        if (iret .ne. nf_noerr) then                             
          write(stdout,2) 'zetadas9_da',indx,incname(1:inclstr)        
          goto 99                                         !--> ERROR                                                             
        endif                                                        
      else                                                          
        write(stdout,1) 'zetadas9_da', incname(1:inclstr)
        goto 99                                           !--> ERROR                                                             
      endif                                                           

#endif

      write(*,*) '... Writing DA total SSH increment'
      iret=nf_inq_varid (ncid, 'zeta_da', varid)                     
      if (iret .eq. nf_noerr) then                                  
        iret=das_nf_fwrite (zeta_s(START_2D_ARRAY), ncid, varid,   
     &                                         indx, r2dvar)      
        if (iret .ne. nf_noerr) then                             
          write(stdout,2) 'zeta_da',indx,incname(1:inclstr)        
          goto 99                                         !--> ERROR                                                             
        endif                                                        
      else                                                          
        write(stdout,1) 'zeta_da', incname(1:inclstr)
        goto 99                                           !--> ERROR                                                             
      endif                                                           

      write(*,*) '... Writing DA steric SSH increment'
      iret=nf_inq_varid (ncid, 'zeta_h_da', varid)                     
      if (iret .eq. nf_noerr) then                                  
        iret=das_nf_fwrite (zeta_h(START_2D_ARRAY), ncid, varid,   
     &                                         indx, r2dvar)      
        if (iret .ne. nf_noerr) then                             
          write(stdout,2) 'zeta_h_da',indx,incname(1:inclstr)        
          goto 99                                         !--> ERROR                                                             
        endif                                                        
      else                                                          
        write(stdout,1) 'zeta_h_da', incname(1:inclstr)
        goto 99                                           !--> ERROR                                                             
      endif                                                           

!
! lat
!
      iret=nf_inq_varid (ncid, 'lat_da', varid)
      if (iret .eq. nf_noerr) then
        iret=das_nf_fwrite (latr(START_2D_ARRAY), ncid, varid,
     &                                         indx, r2dvar)
        if (iret .ne. nf_noerr) then
          write(stdout,2) 'lat_da',indx,incname(1:inclstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'lat_da', incname(1:inclstr)
        goto 99                                           !--> ERROR
      endif
!
!lon
!
      iret=nf_inq_varid (ncid, 'lon_da', varid)
      if (iret .eq. nf_noerr) then
        iret=das_nf_fwrite (lonr(START_2D_ARRAY), ncid, varid,
     &                                         indx, r2dvar)
        if (iret .ne. nf_noerr) then
          write(stdout,2) 'lon_da',indx,incname(1:inclstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'lon_da', incname(1:inclstr)
        goto 99                                           !--> ERROR
      endif
!
! 3D momentum component in the XI-direction.
!
      iret=nf_inq_varid (ncid, 'u_da', varid)
      if (iret .eq. nf_noerr) then
        iret=das_nf_fwrite (u_s(START_2D_ARRAY,1), ncid, varid,
     &                                        indx, u3dvar)
        if (iret .ne. nf_noerr) then
          write(stdout,2) 'u_da', indx,
     &                                     incname(1:inclstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'u_da', incname(1:inclstr)
        goto 99                                           !--> ERROR
      endif
!
! 3D momentum component in the ETA-direction.
!
      iret=nf_inq_varid (ncid, 'v_da', varid)
      if (iret .eq. nf_noerr) then
        iret=das_nf_fwrite (v_s(START_2D_ARRAY,1), ncid, varid,
     &                                        indx, v3dvar)
        if (iret .ne. nf_noerr) then
          write(stdout,2) 'v_da', indx,
     &                                       incname(1:inclstr)
          goto 99                                        !--> ERROR
        endif
      else
        write(stdout,1) 'v_da', incname(1:inclstr)
        goto 99                                           !--> ERROR
      endif

#if defined DAS_UVAGEO_DIAG
      write(*,*) '... Writing ageostrophic u_ap term'
      iret=nf_inq_varid (ncid, 'u_ap_da', varid)
      if (iret .eq. nf_noerr) then
        iret=das_nf_fwrite (u_ap(START_2D_ARRAY,1), ncid, varid,
     &                                        indx, u3dvar)
        if (iret .ne. nf_noerr) then
          write(stdout,2) 'u_ap_da', indx,
     &                                     incname(1:inclstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'u_ap_da', incname(1:inclstr)
        goto 99                                           !--> ERROR
      endif
      write(*,*) '... Writing ageostrophic u_ac term'
      iret=nf_inq_varid (ncid, 'u_ac_da', varid)
      if (iret .eq. nf_noerr) then
        iret=das_nf_fwrite (u_ac(START_2D_ARRAY,1), ncid, varid,
     &                                        indx, u3dvar)
        if (iret .ne. nf_noerr) then
          write(stdout,2) 'u_ac_da', indx,
     &                                     incname(1:inclstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'u_ac_da', incname(1:inclstr)
        goto 99                                           !--> ERROR
      endif

      write(*,*) '... Writing ageostrophic v_ap term'
      iret=nf_inq_varid (ncid, 'v_ap_da', varid)
      if (iret .eq. nf_noerr) then
        iret=das_nf_fwrite (v_ap(START_2D_ARRAY,1), ncid, varid,
     &                                        indx, v3dvar)
        if (iret .ne. nf_noerr) then
          write(stdout,2) 'v_ap_da', indx,
     &                                     incname(1:inclstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'v_ap_da', incname(1:inclstr)
        goto 99                                           !--> ERROR
      endif
      write(*,*) '... Writing ageostrophic v_ac term'
      iret=nf_inq_varid (ncid, 'v_ac_da', varid)
      if (iret .eq. nf_noerr) then
        iret=das_nf_fwrite (v_ac(START_2D_ARRAY,1), ncid, varid,
     &                                        indx, v3dvar)
        if (iret .ne. nf_noerr) then
          write(stdout,2) 'v_ac_da', indx,
     &                                     incname(1:inclstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'v_ac_da', incname(1:inclstr)
        goto 99                                           !--> ERROR
      endif
#endif

!
! Tracer variables.
!
      do itrc=1,NT
        if (itrc .eq. 1) then
           vnts='t_da'
        else if (itrc .eq. 2) then
           vnts='s_da'
        else
           vnts='BIO-VARIABLES'
        endif

        write(*,*) '... Writing tracers'
        write(*,*) vnts
        iret=nf_inq_varid (ncid, vnts, varid)
        if (iret .eq. nf_noerr) then           
          iret=das_nf_fwrite (t_s(START_2D_ARRAY,1,itrc), ncid,  varid,
     &                                               indx, r3dvar)
          if (iret .ne. nf_noerr) then                    
            write(stdout,2) vnts, indx,
     &                                              incname(1:inclstr)
            goto 99                                       !--> ERROR
          endif                                           
        else
          write(stdout,1) vnts, incname(1:inclstr)
          goto 99                                         !--> ERROR
        endif
      enddo

#if defined DAS_BKGZETA
      write(*,*) '... Inquiring DA bkg density var ID'
      iret=nf_inq_varid (ncid, 'rhodas9_da', varid)                     
      if (iret .eq. nf_noerr) then                                  
          write(*,*) '... Writing DA bkg density'
         iret=das_nf_fwrite (rho_das9(START_2D_ARRAY,1), ncid,  varid,
     &                                               indx, r3dvar)
        if (iret .ne. nf_noerr) then                             
          write(stdout,2) 'rhodas9_da',indx,incname(1:inclstr)        
          goto 99                                         !--> ERROR                                                             
        endif                                                        
      else                                                          
        write(stdout,1) 'rhodas9_da', incname(1:inclstr)
        goto 99                                           !--> ERROR                                                             
      endif                                                           
#endif

!    
! mask  
!    
      iret=nf_inq_varid (ncid, 'mask_da', varid)
      if (iret .eq. nf_noerr) then     
        iret=das_nf_fwrite (rmask_das(START_2D_ARRAY,1), ncid,  varid,
     &                                               indx, r3dvar)
        if (iret .ne. nf_noerr) then
          write(stdout,2) 'mask_da', indx,
     &                                              incname(1:inclstr)
          goto 99                                       !--> ERROR
        endif
      else
        write(stdout,1) 'mask_da', incname(1:inclstr)
        goto 99                                         !--> ERROR
      endif
      
  1   format(/1x,'WRT_INC - unable to find variable:',    1x,A,
     &                            /15x,'in input NetCDF file:',1x,A/)
  2   format(/1x,'WRT_INC - error while writing variable:',1x, A,
     &    2x,'at time record =',i4/15x,'in input NetCDF file:',1x,A/)
!    
!  Close input NetCDF file.                              
!       
      iret=nf_close(ncid)
      return
  99  stop 819

      return
!!!!!!!!!!!!!!!!!!!!!!
      end
       
      subroutine check_err(iret)
      integer iret
#include "netcdf.inc"
      if (iret .ne. NF_NOERR) then
      print *, nf_strerror(iret)
      stop
      endif
      end
