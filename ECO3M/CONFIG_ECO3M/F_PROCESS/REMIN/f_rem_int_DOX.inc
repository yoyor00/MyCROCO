!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!m.baklouti@univmed.fr; vincent.faure@univmed.fr
!
!This software (Eco3M) is a computer program whose purpose is to perform 
!biogeochemical or coupled physical-biogeochemical modelling.
!
!This software is governed by the CeCILL license under French law and
!abiding by the rules of distribution of free software. You can  use, 
!modify and/ or redistribute the software under the terms of the CeCILL
!license as circulated by CEA, CNRS and INRIA at the following URL
!"http://www.cecill.info". 
!
!As a counterpart to the access to the source code and  rights to copy,
!modify and redistribute granted by the license, users are provided only
!with a limited warranty  and the software''s author,  the holder of the
!economic rights,  and the successive licensors  have only  limited
!liability. 
!
!In this respect, the user''s attention is drawn to the risks associated
!with loading,  using,  modifying and/or developing or reproducing the
!software by the user in light of its specific status of free software,
!that may mean  that it is complicated to manipulate,  and  that  also
!therefore means  that it is reserved for developers  and  experienced
!professionals having in-depth computer knowledge. Users are therefore
!encouraged to load and test the software''s suitability as regards their
!requirements in conditions enabling the security of their systems and/or 
!data to be ensured and,  more generally, to use and operate it in the 
!same conditions as regards security. 
!
!The fact that you are presently reading this means that you have had
!knowledge of the CeCILL license and that you accept its terms.
!***************************************************************************
!***************************************************************************
!--------------------------------------------------------------------
!
  function f_rem_int_DOX(ili,jcol,indBX,iDOX,iDOC,KDOC,VDOC,QminC,QmaxC)
!
!--------------------------------------------------------------------
! A function for calculating the rate of remineralisation 
! either inside or outside bacterial cells as a function of carbon quota
!
!
! Melika Baklouti Avril 2026
!-------------------------------------------------------------------------

    Implicit None

    !-- fonction:
    Real(8):: f_rem_int_DOX(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: indBX,QminC,QmaxC
    Real(8) :: iDOX,iDOC,VDOC,KDOC
    Integer :: ili,jcol

    !--variables locales:
    Integer :: identite
    Integer :: iiX,iiC,i,j,k,iiCell,iiDOC,iiDOX
    Real(8) :: QC,fQC,hQC,fluxC,fluxX
    REAL    :: nn_C
    Character(L_VAR)::chain

    ! Index of bacterial concentration of element X 
    iiX = nint(indBX)
   
    iiDOX  = nint(iDOX)
    iiDOC  = nint(iDOC)

    ! Organism number for bacteria 
    identite = var(iiX)%idorg

    ! Indices of bacterial concentrations expressed in terms of biomass and abundance
    chain = 'C'
    iiC = f_idorg2id(identite,chain)
    chain = 'cell'
    iiCell = f_idorg2id(identite,chain)
       

    nn_C = 2 

 
    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max  

             ! bacterial carbon cell quota :
                QC =  var(iiC)%conc(i,j,k) / (var(iiCell)%conc(i,j,k)+1.d-80)
                fQC = (QmaxC - QC) / (QmaxC -QminC)
     
                fQC = max(0.d0,fQC)
                fQC = min(1.d0,fQC)
	        fQC = fQC**nn_C
                  
             ! mineralization rate (intracellular)             
		fluxC =  VDOC * var(iiCell)%conc(i,j,k) * var(iiDOC)%conc(i,j,k)/(var(iiDOC)%conc(i,j,k)+KDOC) * fQC 
		fluxX =  fluxC * var(iiDOX)%conc(i,j,k)/(var(iiDOC)%conc(i,j,k)+1d-50) 

                 if ((var(iiDOX)%conc(i,j,k) > fluxX * dt_bio) .and. (var(iiDOC)%conc(i,j,k) > fluxC * dt_bio)) then
                     f_rem_int_DOX(i,j,k) =  fluxX
		 else
                     f_rem_int_DOX(i,j,k) =0.0d0
	         endif
            enddo
        enddo
    enddo
                

end function f_rem_int_DOX
!--------------------------------------------------------------------
