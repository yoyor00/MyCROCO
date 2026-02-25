#include "cppdefs.h"

      subroutine netcdf_filtered(file_nc)

      implicit none

#include "netcdf.inc"
#include "param.h"
#include "das_param.h"
#include "das_ocean.h"
#include "grid.h"

!!!!!!!!!!!!!!!!!!
      integer xi_rho,eta_rho,xi_u,eta_v,ntime,ns
      parameter(ns=nratio/2+1)               ! nratio=1
      parameter(xi_rho=(Lm+2-ns)/nratio+1,   ! xi_rho  Lm+2
     &          eta_rho=(Mm-ns+2)/nratio+1,  ! eta_rho Mm+2
     &  xi_u=xi_rho-1, eta_v=eta_rho-1, ntime=1)

      character file_nc*96
      character start_day*20
      integer i,j,k,i0,j0, lstr, lenstr


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
      integer  t_da_id
      integer  s_da_id
      integer  zeta_da_id
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

* attribute vectors
      real tt(xi_rho,eta_rho,NDAS)                                                                                                                   
      real uu(xi_u,eta_rho,NDAS)
      real vv(xi_rho,eta_v,NDAS)
      real zz(xi_rho,eta_rho)                                                                                                                        
      common /TTT/tt,uu,vv,zz

      real doubleval(1)                                                                                                                              

      integer nc_start(3)
      integer nc_count(3)                                                                                                                            
      integer nc_start2(2)
      integer nc_count2(2)

!
      lstr=lenstr(file_nc)
* enter define mode
      iret = nf_create(file_nc(1:lstr), NF_CLOBBER, ncid)
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
      
      zeta_da_dims(2) =  eta_rho_dim
      zeta_da_dims(1) = xi_rho_dim
      iret = nf_def_var(ncid, 'zeta_da', NF_DOUBLE, zeta_da_rank,
     &       zeta_da_dims, zeta_da_id)
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
      call check_err(iret)
      iret = nf_put_att_text(ncid, v_da_id, 'long_name', 1, 'v')
      call check_err(iret)
!      iret = nf_put_att_double(ncid, v_da_id, 'missing_value', NF_DOUBLE
!     1, 1, doubleval)
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
      iret = nf_put_att_text(ncid, lat_da_id, 'long_name', 8, 'latitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_da_id, 'long_name', 9, 'longitude')
      call check_err(iret)
      
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'type', 11, 'ROMS SCB LR')
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'data', 16, 'ROMS vert interp')
      call check_err(iret)
      
* leave define mode
      iret = nf_enddef(ncid)
      call check_err(iret)

      iret = nf_close(ncid)
      call check_err(iret)


!        print*,'xi_rho', xi_rho
!        print*,'eta_rho', eta_rho
!
! write
      iret = NF_OPEN(file_nc(1:lstr),NF_WRITE,ncid) 

      nc_start(1)=1
      nc_start(2)=1
      nc_start(3)=1
      nc_count(3)=ndas
      nc_count(2)=eta_rho   
      nc_count(1)=xi_u  

      i0=0 
      j0=0 
      do k=1,ndas
        do i=1,xi_u ! das_filtered fills u-point field using IU_RANGE 1:Lm+1
          do j=1,eta_rho
            i0=(i-1)*nratio+1  ! 1:Lm+1 no i-shift
            j0=(j-1)*nratio    ! 0:Mm+1 no j-shift
            uu(i,j,k)=u_s(i0,j0,k)*umask_das(i0,j0,k)
          enddo 
        enddo 
      enddo 
      i0=0 
      j0=0 
      do k=1,ndas
        do i=1,xi_rho
          do j=1,eta_v ! das_filtered fills v-point field using IV_RANGE 1:Mm+1
            i0=(i-1)*nratio   ! 0:Lm+1 no i-shift
            j0=(j-1)*nratio+1 ! 1:Mm+1 no j-shift
            vv(i,j,k)=v_s(i0,j0,k)*vmask_das(i0,j0,k)
          enddo
        enddo
      enddo

      iret = NF_INQ_VARID(ncid, 'u_da',u_da_id)
      if (iret.ne.NF_NOERR) then 
        print*, '   error inquiring u_da'
        stop
      endif
      iret = NF_PUT_VARA_DOUBLE(ncid, u_da_id,
     &                            nc_start,nc_count,uu)
      if (iret.ne.NF_NOERR) then
        print*, '   error write u_da'
        stop 
      endif


      nc_count(3)=ndas
      nc_count(2)=eta_v   
      nc_count(1)=xi_rho
      iret = NF_INQ_VARID(ncid, 'v_da',v_da_id)
      if (iret.ne.NF_NOERR) then
        print*, '   error inquiring v_da'
        stop 
      endif
      iret = NF_PUT_VARA_DOUBLE(ncid, v_da_id,
     &                            nc_start,nc_count,vv)
      if (iret.ne.NF_NOERR) then
        print*, '   error write v_da'
        stop 
      endif

      i0=0 
      j0=0 
      do k=1,ndas
        do i=1,xi_rho
          do j=1,eta_rho
            i0=(i-1)*nratio   ! 0:Lm+1 => no i-shift
            j0=(j-1)*nratio   ! 0:Mm+1 => no j-shift
            tt(i,j,k)=t_s(i0,j0,k,1)*rmask_das(i0,j0,k)
          enddo
        enddo
      enddo
      nc_count(3)=ndas
      nc_count(2)=eta_rho
      nc_count(1)=xi_rho
      iret = NF_INQ_VARID(ncid, 't_da',t_da_id)
      if (iret.ne.NF_NOERR) then
        print*, '   error inquiring t_da'
        stop
      endif
      iret = NF_PUT_VARA_DOUBLE(ncid, t_da_id,
     &                            nc_start,nc_count,tt)
      if (iret.ne.NF_NOERR) then
        print*, '   error write t_da'
        stop
      endif


      i0=0 
      j0=0 
      do k=1,ndas
        do i=1,xi_rho
          do j=1,eta_rho
            i0=(i-1)*nratio
            j0=(j-1)*nratio
            tt(i,j,k)=t_s(i0,j0,k,2)*rmask_das(i0,j0,k)
          enddo
        enddo
      enddo
      nc_count(3)=ndas
      nc_count(2)=eta_rho
      nc_count(1)=xi_rho
      iret = NF_INQ_VARID(ncid, 's_da',s_da_id)
      if (iret.ne.NF_NOERR) then
        print*, '   error inquiring s_da'
        stop
      endif
      iret = NF_PUT_VARA_DOUBLE(ncid, s_da_id,
     &                            nc_start,nc_count,tt)
      if (iret.ne.NF_NOERR) then
        print*, '   error write s_da'
        stop
      endif
!
! mask
!
      i0=0 
      j0=0 
      do k=1,ndas                                                                                                               
        do i=1,xi_rho
          do j=1,eta_rho
            i0=(i-1)*nratio
            j0=(j-1)*nratio
            tt(i,j,k)=rmask_das(i0,j0,k)                                                                       
          enddo                                                                                                                 
        enddo                                                                                                                   
      enddo                                                                                                                     
      nc_count(3)=ndas                                                                                                          
      nc_count(2)=eta_rho                                                                                                       
      nc_count(1)=xi_rho                                                                                                        
      iret = NF_INQ_VARID(ncid, 'mask_da',mask_da_id)                                                                                 
      if (iret.ne.NF_NOERR) then                                                                                                
        print*, '   error inquiring mask_da'                                                                                       
        stop                                                                                                                    
      endif                                                                                                                     
      iret = NF_PUT_VARA_DOUBLE(ncid, mask_da_id,                                                                                  
     &                            nc_start,nc_count,tt)                                                                         
      if (iret.ne.NF_NOERR) then                                                                                                
        print*, '   error write mask_da'                                                                                           
        stop                                                                                                                    
      endif


      i0=0 
      j0=0 
        do i=1,xi_rho
          do j=1,eta_rho
            i0=(i-1)*nratio
            j0=(j-1)*nratio
            zz(i,j)=zeta_s(i0,j0)*rmask_das(i0,j0,ndas)
          enddo
        enddo
      nc_start2(1)=1
      nc_start2(2)=1
      nc_count2(2)=eta_rho    
      nc_count2(1)=xi_rho    
      iret = NF_INQ_VARID(ncid, 'zeta_da',zeta_da_id)
      if (iret.ne.NF_NOERR) then
        print*, '   error inquiring zeta_da'
        stop
      endif
      iret = NF_PUT_VARA_DOUBLE(ncid, zeta_da_id,
     &                            nc_start2,nc_count2,zz)

      i0=0                                                                                                                                           
      j0=0                                                                                                                                           
        do i=1,xi_rho                                                                                                                                
          do j=1,eta_rho                                                                                                                             
            i0=(i-1)*nratio                                                                                                                        
            j0=(j-1)*nratio                                                                                                                        
            zz(i,j)=latr(i0,j0)
          enddo                                                                                                                                      
        enddo                                                                                                                                        
      nc_start2(1)=1                                                                                                                                 
      nc_start2(2)=1                                                                                                                                 
      nc_count2(2)=eta_rho                                                                                                                           
      nc_count2(1)=xi_rho                                                                                                                            
      iret = NF_INQ_VARID(ncid, 'lat_da',lat_da_id)                                                                                                
      if (iret.ne.NF_NOERR) then                                                                                                                     
        print*, '   error inquiring lat_da'                                                                                                         
        stop                                                                                                                                         
      endif                                                                                                                                          
      iret = NF_PUT_VARA_DOUBLE(ncid, lat_da_id,                                                                                                    
     &                            nc_start2,nc_count2,zz)                                                                                            

      i0=0
      j0=0                                                                                                                                           
        do i=1,xi_rho                                                                                                                                
          do j=1,eta_rho                                                                                                                             
            i0=(i-1)*nratio                                                                                                                        
            j0=(j-1)*nratio                                                                                                                        
            zz(i,j)=lonr(i0,j0)
          enddo                                                                                                                                      
        enddo                                                                                                                                        
      nc_start2(1)=1                                                                                                                                 
      nc_start2(2)=1                                                                                                                                 
      nc_count2(2)=eta_rho                                                                                                                           
      nc_count2(1)=xi_rho                                                                                                                            
      iret = NF_INQ_VARID(ncid, 'lon_da',lon_da_id)                                                                                                  
      if (iret.ne.NF_NOERR) then                                                                                                                     
        print*, '   error inquiring lon_da' 
        stop
      endif
      iret = NF_PUT_VARA_DOUBLE(ncid, lon_da_id, 
     &                            nc_start2,nc_count2,zz)
!
       iret = NF_CLOSE(ncid)

!!!!!!!!!!!!!!!!!!!!!!
      end
