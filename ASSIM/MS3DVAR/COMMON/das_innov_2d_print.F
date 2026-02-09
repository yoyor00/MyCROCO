#include "cppdefs.h"

      subroutine das_innov_2D_print
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
      include 'netcdf.inc'

      real aa(Lm+2,Mm+2)
!
      if ( flag_flight .and. flag_mc ) then 

!
! enter define mode
!
        len_file=lenstr(file_innov_2d)
        iret = nf_create(file_innov_2d(1:len_file), NF_CLOBBER, ncid)

        xi_rho=Lm+2
        eta_rho=Mm+2
!
! define dimensions
!
        iret = nf_def_dim(ncid, 'xi_rho', xi_rho, xi_dim)
        iret = nf_def_dim(ncid, 'eta_rho', eta_rho, eta_dim)
!
! define variables
!
        dim2(2)=eta_dim
        dim2(1)=xi_dim

        iret = nf_def_var(ncid, 'fcst_sst', NF_DOUBLE, 2, dim2, fcst_sst_id)
        if ( flag_flight) then
          iret = nf_def_var(ncid, 'fl_sst', NF_DOUBLE, 2, dim2, fl_sst_id)
          iret = nf_def_var(ncid, 'fl_sst_mask', NF_DOUBLE, 2, dim2, fl_mask_id)
        endif
          iret = nf_def_var(ncid, 'mc_sst', NF_DOUBLE, 2, dim2, fl_sst_id)
          iret = nf_def_var(ncid, 'mc_sst_mask', NF_DOUBLE, 2, dim2, mc_mask_id)
! fcst
        iret = nf_put_att_text(ncid, fcst_sst_id, 'long_name', 11, 
     &                        'Forecat SST')
        iret = nf_put_att_text(ncid, fcst_sst_id, 'units', 9, 'degrees C')
!
!innov fl sst
!
        iret = nf_put_att_text(ncid, fl_sst_id, 'long_name', 20, 
     &                        'Innovation Fight SST')
        iret = nf_put_att_text(ncid, fl_sst_id, 'units', 9,'degrees C')

        iret = nf_put_att_text(ncid, fl_mask_id, 'long_name', 24, 
     &                        'Innovation Fight SST mask')
        iret = nf_put_att_text(ncid, fl_mask_id, 'value', 38, 
     &                        'mask=1, with SST; mask=0, without SST ')
!
!innov mc sst
!
        iret = nf_put_att_text(ncid, mc_sst_id, 'long_name', 17, 
     &                        'Innovation MC SST')
        iret = nf_put_att_text(ncid, mc_sst_id, 'units', 9, 'degrees C')

        iret = nf_put_att_text(ncid, mc_mask_id, 'long_name', 21,
     &                        'Innovation MC SST mask')
        iret = nf_put_att_text(ncid, mc_mask_id, 'value', 38,    
     &                        'mask=1, with SST; mask=0, without SST ')

!
! global attributes
!
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'type',11,'INNOVATIONS')
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'fileid', 21, 'MONTEREY ZO
     1OM LEVEL 3')
! leave define mode
      iret = nf_enddef(ncid)
       
      iret =nf_put_vara_double(ncid, fl_sst_id, 2,dim2, sst_fl)
      iret =nf_put_vara_double(ncid, fl_mask_id, 2,dim2, fl_mask)

      iret =nf_put_vara_double(ncid, mc_sst_id, 2,dim2, sst_mc)
      iret =nf_put_vara_double(ncid, mc_mask_id, 2,dim2, mc_mask)
       

      iret = nf_close(ncid)
      
      return
      end

      

