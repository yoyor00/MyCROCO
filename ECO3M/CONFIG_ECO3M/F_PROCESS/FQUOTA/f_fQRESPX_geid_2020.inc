!***************************************************************************
!***************************************************************************
!Copyright or � or Copr. CNRS/IRD/Universit� de la M�diterran�e
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
function f_fQRESPX_geid_2020(ili,jcol,indC,Qmin,Qmax,CP_min,CP_max,CN_min,CN_max)
    !
    !--------------------------------------------------------------------
    ! Quota function for net X uptake (used to calculate the energetic costs)
    !
    ! Melika Baklouti 17/07/2019; last version May 2020
    !-------------------------------------------------------------------------

    Implicit None

    !-- function:
    Real(8):: f_fQRESPX_geid_2020(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: indC,Qmin,Qmax
    Real(8), optional :: CP_min,CP_max,CN_min,CN_max
    Integer :: ili,jcol


    !-- Calcul of the quota function for net uptake :

    f_fQRESPX_geid_2020 = 1.d0 - f_fQPPX_geid_2020(ili,jcol,indC,Qmin,Qmax,CP_min,CP_max,CN_min,CN_max)

end function f_fQRESPX_geid_2020
!--------------------------------------------------------------------