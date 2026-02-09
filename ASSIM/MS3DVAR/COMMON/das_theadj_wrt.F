#include "cppdefs.h"

      subroutine das_theadj_wrt(iflag)

      implicit none
#include "netcdf.inc"
#include "param.h"
#include "das_param.h"
#include "das_ocean.h"
#include "grid.h"
#include "ncscrum.h"

      integer :: iflag
      character(len=6) :: sflag
      character das_start_day*20, vnts*5
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
      integer  psi_da_id
      integer  chi_da_id
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
      integer  psi_da_rank
      integer  chi_da_rank
      integer  mask_da_rank
      integer  lat_da_rank
      integer  lon_da_rank
      integer  zeta_da_rank
      parameter (depth_rank = 1)
      parameter (u_da_rank = 3)
      parameter (v_da_rank = 3)
      parameter (t_da_rank = 3)
      parameter (s_da_rank = 3)
      parameter (psi_da_rank = 3)
      parameter (chi_da_rank = 3)
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
      integer  psi_da_dims(psi_da_rank)
      integer  chi_da_dims(chi_da_rank)
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
      write(sflag,"(i6.6)") iflag
      incname=rstname(1:inclstr)//'_theadj'//sflag//'.nc'
      inclstr=lenstr(incname)
* enter define mode
!     iret = nf_create(rstname(1:lstr), NF_CLOBBER, ncid)
      iret = nf_create(incname(1:inclstr), NF_CLOBBER, ncid)
      print *, '... Before creating adjoint of analysis file'
      print *,incname(1:inclstr)
      call check_err3(iret)
* define dimensions

      iret = nf_def_dim(ncid, 'xi_rho', xi_rho, xi_rho_dim)
      call check_err3(iret)
      iret = nf_def_dim(ncid, 'xi_u', xi_u, xi_u_dim)
      call check_err3(iret)
      iret = nf_def_dim(ncid, 'eta_rho', eta_rho, eta_rho_dim)
      call check_err3(iret)
      iret = nf_def_dim(ncid, 'eta_v', eta_v, eta_v_dim)
      call check_err3(iret)
      iret = nf_def_dim(ncid, 'ndas', ndas, s_rho_dim)
      call check_err3(iret)
* define variables

!      depth_dims(1) = s_rho_dim
!      iret = nf_def_var(ncid, 'depth', NF_DOUBLE, 1, depth_dims
!     1, depth_id)
!      call check_err3(iret)

      write(*,*) '... Defining DA psi adjoint'
!     psi_da_dims(3) = s_rho_dim
!     psi_da_dims(2) = eta_v_dim
!     psi_da_dims(1) = xi_u_dim
      psi_da_dims(3) = s_rho_dim
      psi_da_dims(2) = eta_rho_dim
      psi_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'psi_adj', NF_DOUBLE, psi_da_rank, 
     &       psi_da_dims, psi_da_id)
      call check_err3(iret)

      write(*,*) '... Defining DA chi adjoint'
      chi_da_dims(3) =  s_rho_dim
      chi_da_dims(2) =  eta_rho_dim
      chi_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'chi_adj', NF_DOUBLE, chi_da_rank,
     &       chi_da_dims, chi_da_id)
      call check_err3(iret)

      write(*,*) '... Defining DA T adjoint'
      t_da_dims(3) =  s_rho_dim
      t_da_dims(2) =  eta_rho_dim
      t_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 't_adj', NF_DOUBLE, t_da_rank,
     &       t_da_dims, t_da_id)
      call check_err3(iret)

      write(*,*) '... Defining DA S adjoint'
      s_da_dims(3) =  s_rho_dim
      s_da_dims(2) =  eta_rho_dim
      s_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 's_adj', NF_DOUBLE, s_da_rank,
     &       s_da_dims, s_da_id)
      call check_err3(iret)

      mask_da_dims(3) =  s_rho_dim                                                                                                 
      mask_da_dims(2) =  eta_rho_dim                                                                                               
      mask_da_dims(1) = xi_rho_dim                                                                                                 
      iret = nf_def_var(ncid, 'mask_da', NF_DOUBLE, mask_da_rank,                                                                     
     &       mask_da_dims, mask_da_id)                                                                                                
      call check_err3(iret)

      write(*,*) '... Defining DA total SSH adjoint'
      zeta_da_dims(2) =  eta_rho_dim
      zeta_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'zeta_adj', NF_DOUBLE, zeta_da_rank,
     &       zeta_da_dims, zeta_da_id)
      call check_err3(iret)

      write(*,*) '... Defining DA steric SSH adjoint'
      zeta_da_dims(2) =  eta_rho_dim
      zeta_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'zeta_h_adj', NF_DOUBLE, zeta_da_rank,
     &       zeta_da_dims, zeta_h_da_id)
      call check_err3(iret)

      lat_da_dims(2) =  eta_rho_dim
      lat_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'lat_da', NF_DOUBLE, lat_da_rank,
     &       lat_da_dims, lat_da_id)
      call check_err3(iret)

      lon_da_dims(2) =  eta_rho_dim
      lon_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'lon_da', NF_DOUBLE, lon_da_rank,
     &       lon_da_dims, lon_da_id)
      call check_err3(iret)

* assign attributes

      iret = nf_put_att_text(ncid, t_da_id, 'long_name', 4, 'temp') 
      call check_err3(iret)
!      iret = nf_put_att_double(ncid, t_da_id, 'missing_value', NF_DOUBLE
!     1, 1, doubleval)
      iret = nf_put_att_text(ncid, s_da_id, 'long_name', 4, 'salt') 
      call check_err3(iret)
!      iret = nf_put_att_double(ncid, s_da_id, 'missing_value', NF_DOUBLE
!     1, 1, doubleval)
      iret = nf_put_att_text(ncid, mask_da_id, 'long_name', 4, 'salt') 
      call check_err3(iret)
!      iret = nf_put_att_double(ncid, mask_da_id, 'missing_value', NF_DOUBLE
!     1, 1, doubleval)
      iret = nf_put_att_text(ncid, zeta_da_id, 'long_name', 4, 'zeta') 
      call check_err3(iret)
!      iret = nf_put_att_double(ncid, zeta_da_id, 'missing_value', NF_DOUBLE
!     1, 1, doubleval)
      write(*,*) '... Writing DA steric SSH adjoint attributes'
      iret = nf_put_att_text(ncid, zeta_h_da_id, 'long_name', 4, 'zeta_h') 
      call check_err3(iret)
      write(*,*) '... End writing DA steric SSH adjoint attributes'
      iret = nf_put_att_text(ncid, lat_da_id, 'long_name', 8, 'latitude') 
      call check_err3(iret)
      iret = nf_put_att_text(ncid, lon_da_id, 'long_name', 9, 'longitude') 
      call check_err3(iret)

      iret = nf_put_att_text(ncid, NF_GLOBAL, 'type', 7, 'ROMS LR')
      call check_err3(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'data', 16, 'ROMS vert interp')
      call check_err3(iret)

* leave define mode
      iret = nf_enddef(ncid)
      call check_err3(iret)
       
      iret = nf_close(ncid)
      call check_err3(iret)


!
! write
      indx=1              ! record
      iret = NF_OPEN(incname(1:inclstr),NF_WRITE,ncid) 

      write(*,*) '... Writing DA total SSH adjoint'
      iret=nf_inq_varid (ncid, 'zeta_adj', varid)                     
      if (iret .eq. nf_noerr) then                                  
        iret=das_nf_fwrite (zeta_adj(START_2D_ARRAY), ncid, varid,   
     &                                         indx, r2dvar)      
        if (iret .ne. nf_noerr) then                             
          write(stdout,2) 'zeta_adj',indx,incname(1:inclstr)        
          goto 99                                         !--> ERROR                                                             
        endif                                                        
      else                                                          
        write(stdout,1) 'zeta_adj', incname(1:inclstr)
        goto 99                                           !--> ERROR                                                             
      endif                                                           

      write(*,*) '... Writing DA steric SSH adjoint'
      iret=nf_inq_varid (ncid, 'zeta_h_adj', varid)                     
      if (iret .eq. nf_noerr) then                                  
        iret=das_nf_fwrite (zeta_h_adj(START_2D_ARRAY), ncid, varid,   
     &                                         indx, r2dvar)      
        if (iret .ne. nf_noerr) then                             
          write(stdout,2) 'zeta_h_adj',indx,incname(1:inclstr)        
          goto 99                                         !--> ERROR                                                             
        endif                                                        
      else                                                          
        write(stdout,1) 'zeta_h_adj', incname(1:inclstr)
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

      write(*,*) '... Writing chi'
      iret=nf_inq_varid (ncid, 'chi_adj', varid)
      if (iret .eq. nf_noerr) then
        iret=das_nf_fwrite (chi_adj(START_2D_ARRAY,1), ncid, varid,
     &                                        indx, r3dvar)
        if (iret .ne. nf_noerr) then
          write(stdout,2) 'chi_adj', indx,
     &                                     incname(1:inclstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'chi_adj', incname(1:inclstr)
        goto 99                                           !--> ERROR
      endif

      write(*,*) '... Writing psi'
      iret=nf_inq_varid (ncid, 'psi_adj', varid)
      if (iret .eq. nf_noerr) then
        iret=das_nf_fwrite (psi_adj(START_2D_ARRAY,1), ncid, varid,
     &                                        indx, r3dvar)
        if (iret .ne. nf_noerr) then
          write(stdout,2) 'psi_adj', indx,
     &                                     incname(1:inclstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'psi_adj', incname(1:inclstr)
        goto 99                                           !--> ERROR
      endif

!
! Tracer variables.
!
      do itrc=1,NT
        if (itrc .eq. 1) then
           vnts='t_adj'
        else if (itrc .eq. 2) then
           vnts='s_adj'
        else
           vnts='BIO-VARIABLES'
        endif

        write(*,*) '... Writing tracers'
        write(*,*) vnts
        iret=nf_inq_varid (ncid, vnts, varid)
        if (iret .eq. nf_noerr) then           
          iret=das_nf_fwrite(t_adj(START_2D_ARRAY,1,itrc),ncid,varid,
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
      
  1   format(/1x,'WRT_ADJ - unable to find variable:',    1x,A,
     &                            /15x,'in input NetCDF file:',1x,A/)
  2   format(/1x,'WRT_ADJ - error while writing variable:',1x, A,
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
       
      subroutine check_err3(iret)
      integer iret
#include "netcdf.inc"
      if (iret .ne. NF_NOERR) then
      print *, nf_strerror(iret)
      stop
      endif
      end
