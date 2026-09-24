!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée (now
!Aix-Marseille University)
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!melika.baklouti@univ-amu.fr; 
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
!---------------------------------------------------------------------------
                     Module mod_eco3m_irrad

! Module specfic to light irradiance and extinction
!--------------------------------------------------------------------------------
    use mod_eco3m_vartypes  ! Module for variable types

implicit none
    real(8) :: dt_irrad  !< Time-step for the irradiance reading/calculation 
    real(8) :: irr2par   !< Conversion factor from irradiance to par (set to 1 if PAR is used instead of irradiance) 
    real(8) :: irr_param  !<  Parameter used to convert PAR(0+) into PAR(0-)
    Character(15) :: fichirrad  !< Definition of irradiance calculation.
    !! ( among IRR_FONCTION, IRR_filename, IRR_CODEPHYS or PAR_CODEPHYS)
    Character(40) :: fichirrad_long  !< Name of the irradiance file (if
    !! irradiance is read from a file)
    REAL(8) :: irrad_MAX  !< Maximum irradience (used if fichirrad=IRR_FONCTION)
    TYPE(PARAM) :: IRR_PAR  !< Parameters of the irradiance function
    REAL(8), pointer :: irrad(:,:)  !< 2D irradiance
    real(8) :: albedo  !< Albedo

CONTAINS

End Module mod_eco3m_irrad
