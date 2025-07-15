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
                   Module mod_eco3m_outputs

! Module specific to Eco3M outputs
!--------------------------------------------------------------------------------
    
implicit none

    character(124) :: output_dir  !< Output directory (namelist)
    integer :: num_outphy  !< Number of Physical variables that will be saved (namelist)
    integer :: num_outflux !< Number of points at which biogeochemical fluxes will be saved (namelist)
    
    type outflux_type  !< Type associated with the saving of fluxes (namelist)
        integer :: indi  !< i-index 
        integer :: indj  !< j-index
        integer :: indk  !< k-index
    end type
    
    type outphy_type  !< Type associated with the saving of physical variables (namelist)
        character(124) :: outvarname    !< Name of the variable to output 
    end type

    type(outflux_type), allocatable :: outflux(:) !< Array of coordinates where to save fluxes
    type(outphy_type), allocatable :: outphy(:)   !< Array of physical variables to save

    integer :: filecount !< File counter (iterated after each writing in a file)

CONTAINS

end Module mod_eco3m_outputs
