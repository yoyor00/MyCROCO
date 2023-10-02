!***************************************************************************
!***************************************************************************
!Copyright or � or Copr. CNRS/IRD/Universit� de la M�diterran�e
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/n006)
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
function f_rem_DOM(ili,jcol,indBX,indDOM,Qmin,Qmax,Vmax,Kdemi)
    !
    !--------------------------------------------------------------------
    ! This function will control  DON and DOP mineralization by bacteria
    ! 
    ! La fonction de quota est donn�e par le mod�le de Geider et al. (1998) 
    !
    ! Variables et parametres et leurs unites:
    !-----------------------------------------
    ! indBX: index of X concentration in Bacteria
    ! indDOM: index of DOM to be mineralized
    !
    ! Ref: Geider et al.(1998) Limnol. & Oceanogr. 43(4):679-694
    !
    ! Melika Baklouti 7/06/2016
    !-------------------------------------------------------------------------

    Implicit None

    !-- fonction:
    Real(8):: f_rem_DOM(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: indBX,indDOM,Qmin,Qmax,Vmax
    Real(8),optional::Kdemi
    Integer :: ili,jcol,iindX

    !--variables locales:
    Integer :: identite
    Integer :: i,j,k,iiCell,iiBX,iiDOM,iiN,iiP
    Real(8) :: Q,fQ,nn
    Character(L_VAR)::chain,chain2,chain3

    ! Indice de l'element a contraindre
       iiBX = nint(indBX)
       iiDOM = nint(indDOM)
    ! Numero de l'Organismes
       identite = var(iiBX)%idorg
       chain = 'cell'
    chain2 = 'N'
    chain3 = 'P'
       iiCell = f_idorg2id(identite,chain)
       iiN =  f_idorg2id(identite,chain2)
       iiP =  f_idorg2id(identite,chain3)

      if (iiBX == iiP) then 
          nn = 0.06
      elseif(iiBX==iiN) then
          nn = 0.6
     else
        write(*,*) 'pb ds f_rem_DOM'
     endif


    !-- Calcul de la fonction de limitation :
if (present(Kdemi)) then
    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max ! good position
                Q = var(iiBX)%conc(i,j,k) /  (var(iiCell)%conc(i,j,k)+1.d-80)

                ! Calcul de fQ d'uptake
                if (Q .LE. Qmax) then
                    fQ = ((Qmax - Q) / (Qmax - Qmin))**nn      
                else
                    fQ = 0.d0;
                endif
 
!-- Is there sufficient DOM to be mineralized during 20 minutes ?
                if (VAR(iiDOM)%conc(i,j,k) > Vmax * VAR(iiDOM)%conc(i,j,k)/ (VAR(iiDOM)%conc(i,j,k) + Kdemi) &
                                                                         * var(iiCell)%conc(i,j,k) * dt_bio) then
                    f_rem_DOM(i,j,k) = Vmax *  VAR(iiDOM)%conc(i,j,k)/ (VAR(iiDOM)%conc(i,j,k) + Kdemi) &
                                                                         * var(iiCell)%conc(i,j,k)  * fQ
                else
                    f_rem_DOM(i,j,k) = 0.d0
                endif
              
            enddo
        enddo
    enddo

else
    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max ! good position
                Q = var(iiBX)%conc(i,j,k) /  (var(iiCell)%conc(i,j,k)+1.d-80)

                ! Calcul de fQ d'uptake
                if (Q .LE. Qmax) then
                    fQ = ((Qmax - Q) / (Qmax - Qmin))**nn      
                else
                    fQ = 0.d0;
                endif
 
!-- Is there sufficient DOM to be mineralized during 20 minutes ?
                if (VAR(iiDOM)%conc(i,j,k) > Vmax * var(iiCell)%conc(i,j,k) * dt_bio) then
                    f_rem_DOM(i,j,k) = Vmax * var(iiCell)%conc(i,j,k)  * fQ
                else
                    f_rem_DOM(i,j,k) = 0.d0
                endif
              
            enddo
        enddo
    enddo
endif
    !WRITE(*,*) 'FUNCTION f_rem_DOM '
end function f_rem_DOM
!--------------------------------------------------------------------