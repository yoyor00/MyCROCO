#include "cppdefs.h"

      subroutine das_copy_steric (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_copy_steric_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_copy_steric_tile (Istr,Iend,Jstr,Jend) 
!
! compute ssh increments by integration of the
! linearized hydrostatic equation
!
      implicit none
# include "param.h"
# include "das_param.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_ocean9.h"
# include "scalars.h"

      integer Istr,Iend,Jstr,Jend, i,j,k 
      integer is,ie,js,je,i0,j0
      real rr,rc,ro
      integer count,num,n0
!
      integer srad,sdim
      parameter (srad=sm_rad, sdim=(2*srad+1)*(2*srad+1))
      real zz(sdim),dis(sdim)
!

# include "compute_extended_bounds.h"
!
! rho is the deviation from the Boussinesque approximation mean density
! total density is rho+rho0 (rho0=1000kg/m^3)
!
#define DISTANCES
#ifdef DISTANCES                                                                                                                         
      ro = sqrt ((lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))                                                                                   
     &          *(lonr(Lm/2,Mm/2)-lonr(Lm/2+1,Mm/2+1))                                                                                   
     &          +(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1))                                                                                   
     &          *(latr(Lm/2,Mm/2)-latr(Lm/2+1,Mm/2+1)))                                                                                  
#endif                                                                                                                                   
!                                                                                                                                        
      do j=JstrR,JendR
        do i=IstrR,IendR
          if (rmask_das(i,j,NDAS) .gt. 0.5) then                                                                                         
            count=0                
            js=max(j-srad, 0)     
            je=min(j+srad, Mm+1) 
            is=max(i-srad,0)    
            ie=min(i+srad, Lm+1)      
            do j0=js,je              
              do i0=is,ie           
                if (rmask_das(i0,j0,NDAS) .gt. 0.5) then                                                                                 
#ifdef DISTANCES                                                                                                                         
                  rr=sqrt( (lonr(i,j)-lonr(i0,j0)) 
     &                    *(lonr(i,j)-lonr(i0,j0))
     &                    +(latr(i,j)-latr(i0,j0))
     &                    *(latr(i,j)-latr(i0,j0)) )                                                                                     
#endif                                 
                  count=count+1         
                  zz(count)=zeta_h(i0,j0)
#ifdef DISTANCES         
                  dis(count)=max(ro,rr)                                                                                                  
#endif                    
                endif       
              enddo        
            enddo                                                                                                                        
            rc=0.0        
            zeta_s(i,j)=0.0                                                                                                             
            do n0=1,count  
#ifdef DISTANCES                                                                                                                         
              rc=rc+1./dis(n0) 
              zeta_s(i,j)=zeta_s(i,j)+zz(n0)/dis(n0)
#else                        
              rc=rc+1.      
              zeta_s(i,j)=zeta_s(i,j)+zz(n0)       
#endif                                            
            enddo                                                                                                                        
            zeta_s(i,j)=zeta_s(i,j)/rc   
          else                          
            zeta_s(i,j)=0.0                                                                                                             
          endif               
        enddo                                                                                                                            
      enddo
      return
      end

