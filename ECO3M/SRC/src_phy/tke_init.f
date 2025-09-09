c!-----------------------------------------------------------------------------
      subroutine init_phys_model


c!> Initialises the Tke Model 
c!!
c!! - Call of the *RUNDEF* subroutine (reading of the *runparam.i* file)
c!! - Call of the *WRUNDEF* subroutine
c!! - Call of the *RUNGRID* subroutine (reading of the *cgrille* file)
c!! - Call of the *WRUNGRID* subroutine
c!! - Call of the *DYNDEF* subroutine (reading of the *dynparam.i* file)
c!! - Call of the *WDYNDEF* subroutine
c!! - Dynamical allocation of the *TENEUR*, *DIFFBIO*, *SEDBIO* and *SEDFLUX* arrays
c!! - Call of the *BIODEF* subroutine
c!! - Initialisation of the TKE time-stepping variables 
c!! - Call of the *DYNINIT* subroutine
c!! - Call of the *WDYNINIT* subroutine 
c!! - Call of the *BIOINIT* subroutine
c!! - The Eco3M value of albedo is replaced by the TKE albedo value
c!! - The Eco3M value of *bathymetry* is initialised by the *DEPW(nz_max)* values of the TKE model
c!! - The Eco3M value of *depth_z* is initialised by the TKE *DEPT* variable
c!! - The Eco3M value of *dz* is initialised by the TKE *E3T* variable
c!!
c!! \author Melika Baklouti
c!!
c!! \warning In the *biodef_eco3m.f* file, one and only one sedimentation 
c!! speed is read in the file. This value is then used for each "det" compartment.
c!! If your detritic department has another name, then modify the 
c!! *biodef_eco3m.f* file accordingly
      use mod_eco3m
      use comrunmod
      use mod_eco3m_outputs
      use mod_varphy_coupl
      use comdynmod, only: albedo_phy,DE3T
      use mod_eco3m_irrad, only: albedo

      implicit none


      integer :: m, n, l  ! Loop indexes
      integer:: istat
        
c! Initialising the jptract variable with the Eco3M variable
        jptract = nbvar
        Allocate(CNMTRA(jptract))
        Allocate(MSTOCK(jptract))
        Allocate(MBIODIF(jptract))
        Allocate(BIOAVER(jpzt,jptract),BIOAIR(jptract),
     &FLUXAIR(jptemps,jptract))
        Allocate(VITSED0(jptract),VITSED(jptract))
        Allocate(MNEG(jpzt,jptract),MNEGB(jpzt,jptract))

        if (jpzt .ne. nz_max) then 
            write(*,*) "***********************************************"
            write(*,*) "The number of vertical levels is inconsistent"
            write(*,*) "TKE model: ", jpzt
            write(*,*) "Eco3M model: ", nz_max
            write(*,*) "This program will stop"
            write(*,*) "***********************************************"
            stop
        endif

        if (pos_nzmax .ne. "BOTT") then
            write(*,*) "**************************************"
            write(*,*) "In the TKE model, the bottom is at the
     & nz_max level. It has therefore been switched from ", 
     &       pos_nzmax, " to BOTT"
            write(*,*) "**************************************"
            pos_nzmax = "BOTT"
        end if
        
c! Allocation of dynamical table
      Allocate(teneur(1:jpzt,1:nbvar),stat=istat)
      if (istat /=0) write(*,*) "problem in the allocation of teneur"
      Allocate(diffbio(1:jpzt,1:nbvar), stat=istat)
      if (istat /=0) write(*,*) "problem in the allocation of diffbio"
      Allocate(sedbio(1:jpzt,1:nbvar),stat=istat) 
      if (istat /=0) write(*,*) "problem in the allocation of sedbio"
      Allocate(sedflux(1:jpzt,1:nbvar),stat=istat)
      if (istat /=0) write(*,*) "problem in the allocation of sedflux"


c        ! Output file where the outputs of the W* subroutines will be written
        open(99, file='biotke.check')

c        ! Reading of the run definition (../CONFIG_PHY/tke/inputs/runparam.i file)
c        ! Contains general parameters (config name, time-step, etc)
        call RUNDEF
c        ! Writting out the run definition 
        call WRUNDEF

c        ! Reading of the grid variable (../CONFIG_PHY/tke/inputs/grille300.i file)
c        ! Depth for each grid cell)
        call RUNGRID
c        ! Writting out the grid definition
        call WRUNGRID

c        ! Reading of the physical parameters of the TKE model (../CONFIG_PHY/tke/inputs/dynparam.i file)
        call DYNDEF
        
c        ! Writting out the physical parameters
        call WDYNDEF

c        ! reading of some biogeochemical variables
c        ! (./PHY/inputs/bioparam.i file)
c        ! contains for instance the location of some
c        ! forcing inputs files 
        call BIODEF

        nsa = 0
        nstep = 0
        nstepneg = 0
        nstepneg0 = 0
        time = 0.0
        timep = tdebut
        nwrite0 = 0
        timeyr = float(int(timep/(xyear*day*hourm)))+1.
        timeday = timep/(hourm*day) - (timeyr-1.)*xyear+1.
        ntimeday = int(timeday+0.001)
        timeminu = (timeday - ntimeday)*day*hourm
        ntimeminu = int(timeminu + 0.001)

c        ! Initialisation of dynamical profiles (Temperature, salinity, velocities) 
        call DYNINIT
c        ! writting out the initialized dynamical profiles
        call WDYNINIT
c        ! Initialisation of the biological profiles
c        ! biological variables except those from Eco3M-MED.
c        ! These variables are used mainly for saving the outputs
        call BIOINIT
        
        albedo = albedo_phy
        
c        !-- profondeur (hauteur d''eau totale de la colonne d''eau):
        do m = ny_min, ny_max
            do n = nx_min, nx_max
                bathymetry(n, m) = DEPW(nz_max)
            enddo
        enddo
        
c        ! Cell altitudes (reference to the bottom)
        do m = nx_min, nx_max
            do n = ny_min, ny_max
                do l = 1, nz_max
                    depth_z(n, m, l) = DEPT(l)
                enddo
            enddo
        enddo
        
c        ! Cell dz (reference to the bottom)
        do m = nx_min, nx_max
            do n = ny_min, ny_max
                do l = 1, nz_max
                    dz(n, m, l) = DE3T(l)
                enddo
            enddo
        enddo

      if(dt_bio /= dts) then
         write(*,*) "****************************************"
         write(*,*) "The time-step of the TKE model, ", dts, 
     &   " is different from the time-step of the Eco3M model ",dt_bio
         write(*,*) "We set the Eco3M time-step equal to the TKE value"
         write(*,*) "****************************************"
         dt_bio = dts
      end if

      end subroutine init_phys_model

