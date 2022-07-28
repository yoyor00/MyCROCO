!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
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
!> This file contains subroutines dedicated to the saving of Eco3M variables. It also includes 
!! subroutine dedicated to the saving of Physical variables.
!! \author Melika Baklouti, Vincent Faure, Nicolas Barrier
!! \date 2016-02-16

    
!-------------------------------------------------------------------------------
    subroutine eco3m_outputs_init

    !> Initialization of the output variables. 
    !!
    !! - reads in the namelist the output parameters (output dir, number of physical variables
    !! and fluxes to be saved)
    !! - if fluxes are to be saved: dynamical allocation of the *outflux* array and filling 
    !! - if physical variables are to be saved: dynamical allocation of the *outphy* array
    !!  and reading its values in the table
    !! - opens all the input files
    !!
    !! \author Nicolas Barrier
    !! \date 2016-02-04
!-------------------------------------------------------------------------------
        use mod_eco3m_outputs
        use mod_eco3m_files, only: namelist_eco3m_id,eco3m_root_dir
        implicit none 

        ! local variables
        namelist /nam_eco3m_outpar/num_outphy,output_dir
        namelist /nam_eco3m_outflux_nb/num_outflux
        namelist /nam_eco3m_outflux_loc/outflux
        namelist /nam_eco3m_outphy/outphy

        ! reads the output parameters
        rewind(namelist_eco3m_id)
        read(namelist_eco3m_id, nam_eco3m_outpar)
        output_dir=trim(adjustl(eco3m_root_dir))//trim(adjustl(output_dir))

#ifdef SAVE_FLUX
        rewind(namelist_eco3m_id)
        read(namelist_eco3m_id, nam_eco3m_outflux_nb)
        ! if the number of biogeochemical fluxes to be saved is not 0,
        ! allocation of the list
        if (num_outflux > 0) then
            allocate(outflux(num_outflux))
            rewind(namelist_eco3m_id)
            read(namelist_eco3m_id, nam_eco3m_outflux_loc)
        end if
#else
        num_outflux = 0
#endif

        ! if the number of physical variables to be saved is not 0,
        ! allocation of the list
        if (num_outphy.gt.0) then
            allocate(outphy(num_outphy))
            rewind(namelist_eco3m_id)
            read(namelist_eco3m_id, nam_eco3m_outphy)
        end if
#if defined CALC
        call eco3m_open_outputfiles
#endif
        
    end subroutine eco3m_outputs_init
!-------------------------------------------------------------------------------
    subroutine eco3m_open_outputfiles

    !> Opens output files:
    !!
    !! - in NCOUPL mode, opening of the budget files
    !! - opening of the state variable files (i.e. concentration)
    !! - if *num_outphy* is not 0, opening of the physical variable files
    !! - if *num_outflux* is not 0, opening of the flux files
    !!
    !! \author Melika Baklouti, Vincent Faure, Nicolas Barrier
    !! \date 2016-02-04
!-------------------------------------------------------------------------------
        use mod_eco3m
        use mod_eco3m_outputs
        use mod_eco3m_files, only: file_budget_id,file_output_id

        implicit none

        ! local variables
        integer :: iflx  ! Loop index
        character(124) :: filename, varstring
        character(l_chain) :: stri, strj, strk
        integer :: longchain,file_flx_id

        ! opening of the output budget file
        filename = trim(output_dir) // trim(run_name) // "_" // "budgets.out" 
        open(file_budget_id, file=trim(filename))

        ! Initialization of the file counter
        filecount = file_output_id 

#if !defined ECO3M_SUB 
  ! opens the files for state variable outputs
  ! loops over all the state variables
        do iflx = 1, nbvar
            ! Definition of the var string (example: mid_NO3_N, etc)
            varstring = trim(VAR(iflx)%comp) // "_" // trim(VAR(iflx)%scomp) &
                // "_" // trim(var(iflx)%elmt)
            ! Filename is for example: output_dir/runname_var_mid_NO3_N.dat
            filename = trim(output_dir) // trim(run_name) // "_var_" // &
                trim(varstring) // ".out"
            ! opening of the file. If already existing, it is erased
            open(filecount, file=trim(filename), status="replace")
            ! writes  the header 
            write(filecount, *) "% Eco3M run: ", trim(run_name)  ! Writing the run name
            write(filecount, *) "% Eco3M state variable ", trim(varstring)  ! Variable name
            write(filecount, *) "% Time (h), C(k=1,nzmax)"  ! Write column headings
            filecount = filecount + 1  ! Iteration of the file counter
        end do 
    
  ! if variables other than the state variables of the biogeochemical Eco3M model (e.g. PAR value) are to be saved
        if(num_outphy > 0) then
            ! Looping over all the variables to save
            do iflx = 1, num_outphy
                ! definition of filename, same format as for state variables
                filename = trim(output_dir) // &
                    trim(run_name) // "_var_" // &
                    trim(outphy(iflx)%outvarname) // ".out"
                open(filecount, file=trim(filename))
                ! writing of the header
                write(filecount, *) "% Eco3M run: ", trim(run_name)  ! Writing the run name
                write(filecount, *) "% Eco3M variable ", trim(outphy(iflx)%outvarname)  ! Variable name
                write(filecount, *) "% Time (h), var(k=1,nzmax)"  ! Write column heading 
                filecount = filecount + 1 ! Iteration of the file counter
            end do
        end if
#endif /* ! defined ECO3M_SUB */

#ifdef SAVE_FLUX
  ! if Fluxes are to be saved
 if(num_outflux > 0) then
  ! Loop over all the elements of the outflux array
        do iflx = 1, num_outflux
  ! conversion of indi, indj and indk into strings (in order to use in file)
          stri = f_Int2chain(outflux(iflx)%indi, longchain)
          strj = f_Int2chain(outflux(iflx)%indj, longchain)
  ! construction of the file names
          strk = f_Int2chain(outflux(iflx)%indk, longchain)
          filename = trim(output_dir) // trim(run_name) // "_flux_" // &
                 trim(stri) // "_" // trim(strj) // "_" // &
                 trim(strk) // ".out"
  ! file number and opening
         file_flx_id = 5000 + outflux(iflx)%indi + outflux(iflx)%indj + outflux(iflx)%indk
          open(file_flx_id, file=trim(filename), status="replace")

  ! writing of the Header
         write(file_flx_id, '(2A)') " % Eco3M run: ", trim(run_name)  ! Writing the run name
         write(file_flx_id, "(A)") " % Eco3M flux variable :" 
         write(file_flx_id, "(A)") " % i-index j-index k-index :" 
         write(file_flx_id, "(A, 3I5.3)") " % ", outflux(iflx)%indi, &
                   outflux(iflx)%indj, outflux(iflx)%indk
     end do
 end if
#endif
  ! reinitialization of the filecounter to the first output file
        filecount = file_output_id

    end subroutine eco3m_open_outputfiles
!-------------------------------------------------------------------------------
    subroutine eco3m_write_varindex

    !> Writes into an output file the index and names of the state variables. 
    !! Usefull for plotting
    !! \author Nicolas Barrier
    !! \note The file is opened and closed in this function
!-------------------------------------------------------------------------------
        use mod_eco3m
        use mod_eco3m_outputs
        implicit none

        ! local variables
        integer :: ivar  ! Loop index
        character(124) :: filename  ! Name of the output file
        character(124) :: varname
        integer :: file_id !identification number for file

        file_id = 222
        ! Opening of a file containing the index of the different state variables 
        filename = trim(output_dir) // trim(run_name) // "_varindex.dat"

        ! Opening of the file, writing of the header
        open(file_id, file=filename, status="replace")
        write(file_id, *) "% Saving of the variable indexes"

        ! Looping over all the state variables and writing of their name and their index
        do ivar = 1, nbvar
            varname = trim(VAR(ivar)%comp) // "_" &
                // trim(VAR(ivar)%scomp) // "_" // trim(VAR(ivar)%elmt)
            write(file_id, "(I5, A,  A)") ivar, "    ", varname
        end do

        ! Closing the file
        close(file_id)

    end subroutine eco3m_write_varindex
#if defined M3D_NCOUPL  || (defined COUPL && (! defined ECO3M_SUB))
!-------------------------------------------------------------------------------
    subroutine eco3m_write_deptharray
    
    !> Writes into an output file the *depth_z* array, i.e 
    !! the distance from the water surface. This subroutine is
    !! called at the end of the run. It is used for plotting output data. 
    !! \author Nicolas Barrier
    !! \note The file is opened and closed in this function
!-------------------------------------------------------------------------------
        use mod_eco3m,only:nz_max,depth_z,run_name
        use mod_eco3m_outputs
        implicit none

        ! local variables
        integer :: ivar  ! Loop index
        character(124) :: filename  ! Name of the output file
        integer :: file_id !identification number for file

        file_id = 222
        ! Opening of a file containing the depth array
        filename = trim(output_dir) // trim(run_name) // "_depth.dat"

        ! Opening of the file, writing of the header
        open(file_id, file=filename)
        write(file_id, *) "% Saving of the depth_z array (reference to the surface)"

        ! Looping over all the depths intervals and writing of the depth_z array
        do ivar = 1, nz_max
            write(file_id, *) ivar, abs(depth_z(1, 1, ivar))
        end do

        ! Closing the file
        close(file_id)

    end subroutine eco3m_write_deptharray
#endif
!-------------------------------------------------------------------------------
#if defined CALC && !defined ECO3M_SUB
    subroutine eco3m_write_outputs

    !> Writes into all the input files. This function is called any time the current time step
    !! *tps* is  a multiple of the saving time step *dt_save_bio*
    !!
    !! - If NCOUPL, it writes the budget file
    !! - Writes the state variables concentrations
    !! - If *num_outphy* is not 0, it writes the values of the physical variables specified
    !! in the *outphy* array.
    !!
    !! \author Nicolas Barrier
    !! \date 2016-02-04
    !! \note The saving of fluxes is made in the *call_save_flux.inc* file.
!-------------------------------------------------------------------------------
        use mod_eco3m_outputs
        use mod_eco3m_files, only:file_output_id
        implicit none

      ! In UNCOUPLED mode, checking budgets   
        call eco3m_write_budgets
      
      ! Reinitialisation of the file counter
        filecount = file_output_id

      ! Calling of the saving of biological variables
        call eco3m_write_varoutputs

      ! If physical variables are to be saved
        if (num_outphy.gt.0) then
            call eco3m_write_glooutputs
        end if

    end subroutine  eco3m_write_outputs
#endif
!-------------------------------------------------------------------------------
#if defined CALC && ! defined ECO3M_SUB
    subroutine eco3m_write_varoutputs

    !> Writes the concentrations of each state variables into
    !! the corresponding output files. 
    !! \author Nicolas Barrier
    !! \date 2016-02-04
!-------------------------------------------------------------------------------
        use mod_eco3m, only:VAR, nbvar,nz_max
        use mod_eco3m_outputs
        implicit none

        ! local variables
        integer :: ivar  ! Loop index

        ! Looping over each state variable
        do ivar = 1, nbvar
            ! We write the vertical concentration profile into the file
            call eco3m_write_2d_var(var(ivar)%conc(1, 1, 1:nz_max))

            ! Iteration of the file counter
            filecount = filecount + 1  
        end do

    end subroutine eco3m_write_varoutputs
#endif
!-------------------------------------------------------------------------------
#ifdef CALC 

    subroutine eco3m_write_glooutputs

    !> Saves the Eco3M physical variables (i.e no state variables) into the 
    !! corresponding output file. Variables that can be saved are:
    !!
    !! - e_par (surface photo. available radiation)
    !! - irrad (surface irradiance)
    !! - e_parz (3D photo. available radiation)
    !! - kextinc (light extinction coefficients)
    !! - e_parz_r (3D photo. available radiation in RED) **Only if key RGB**
    !! - kextinc_r (light extinction coefficients in RED) **Only if key RGB**
    !! - e_parz_g (3D photo. available radiation in GREEN) **Only if key RGB**
    !! - kextinc_g (light extinction coefficients in GREEN) **Only if key RGB**
    !! - e_parz_b (3D photo. available radiation in BLUE) **Only if key RGB**
    !! - kextinc_b (light extinction coefficients in BLUE) **Only if key RGB**
    !! - temp_bio (water temperature) **Only of key COUPL**
    !! - sal_bio (water salinity) **Only of key COUPL**
    !! - raut (water density) **Only of key COUPL**
    !! - wind (wind speed) **Only of key COUPL**
    !! - wind (wind speed) **Only of key COUPL**
    !! - windx (zonal wind speed) **Only of key COUPL**
    !! - windy (meridional wind speed) **Only of key COUPL**
    !! 
    !! If the variable name (defined in the namelist) doesnt match any case, nothing is written in
    !! the output files. 
    !!
    !! \author Nicolas Barrier
    !! \date 2016-02-04
!-------------------------------------------------------------------------------
        use mod_eco3m
        use mod_eco3m_outputs
        use mod_eco3m_irrad
#if defined COUPL && defined key_tke_eco3m
        use comdynmod, only:raut,wind,windx,windy,AVM,UN,VN,DEPTHML,EN,nind,nsave,xdml
#endif
        implicit none 

        ! local variables
        integer :: ivar  ! Loop index
        Real(8):: Chl_tot(nx_min:nx_max,ny_min:ny_max,1:nz_max)
#if ! defined ECO3M_SUB
        do ivar = 1, num_outphy  
            select case (outphy(ivar)%outvarname)
            case ("e_par")
                call eco3m_write_1d_var(e_par(1,1))
            case ("irrad")
                call eco3m_write_1d_var(irrad(1,1))
            case ("chl_tot")
                call compute_chl_tot(Chl_tot)
                call eco3m_write_2d_var(Chl_tot(1, 1, 1:nz_max))
#if ! defined M0D_NCOUPL
            case ("e_parz")
                call eco3m_write_2d_var(E_PARZ(1, 1, 1:nz_max))
            case ("kextinc")
                call eco3m_write_2d_var(kextinc(1, 1, 1:nz_max))
#endif
#ifdef RGB
            case ("e_parz_r")
                call eco3m_write_2d_var(e_parz_rgb(3,1,1,1:nz_max))
            case ("kextinc_r")
                call eco3m_write_2d_var(kextinc_rgb(3,1,1,1:nz_max))
            case ("e_parz_g")
                call eco3m_write_2d_var(e_parz_rgb(2,1,1,1:nz_max))
            case ("kextinc_g")
                call eco3m_write_2d_var(kextinc_rgb(2,1,1,1:nz_max))
            case ("e_parz_b")
                call eco3m_write_2d_var(e_parz_rgb(1,1,1,1:nz_max))
            case ("kextinc_b")
                call eco3m_write_2d_var(kextinc_rgb(1,1,1,1:nz_max))
#endif

            case ("temp")
                call eco3m_write_2d_var(temp_bio(1,1,1:nz_max))
            case ("sal")
                call eco3m_write_2d_var(sal_bio(1,1,1:nz_max))
#endif
#if defined COUPL &&  defined key_tke_eco3m 
            case("raut")
                call eco3m_write_2d_var_r4(raut)
            case("wind")
                call eco3m_write_1d_var_r4(wind)
            case("windx")
                call eco3m_write_1d_var_r4(windx)
            case("windy")
                call eco3m_write_1d_var_r4(windy)
            case("Kz")
                call eco3m_write_2d_var_r4(AVM)
            case ("vitx")
                call eco3m_write_2d_var_r4(UN)
            case ("vity")
                call eco3m_write_2d_var_r4(VN)
            case ("mld")
!                call eco3m_write_1d_var_r4(DEPTHML(nind)/nsave)
                call eco3m_write_1d_var_r4(xdml)
            case ("tke")
                call eco3m_write_2d_var_r4(EN)
#endif
#if ! defined ECO3M_SUB
            end select

            filecount = filecount + 1  ! Iteration of the file counter

        end do  
#endif
    end subroutine eco3m_write_glooutputs
#endif /*CALC*/
!-------------------------------------------------------------------------------
#ifdef SAVE_FLUX
    subroutine eco3m_write_fluxoutputs(ss_flux, ili, jcol, nbfunc_loc)

    !> Saves biological fluxes of fluxes if *num_outflux* is not 0.
    !! This function is called in the *call.inc* filefor each flux defined in the *config.ini* file
    !!\author Melika Baklouti, Vincent Faure
    !! \author Nicolas Barrier
    !! \date 2016-02-04
!-------------------------------------------------------------------------------
        use mod_eco3m
        use mod_eco3m_outputs
        implicit none 

        ! local variables

        real(8),intent(in) :: ss_flux(nbfunc_loc+1,nx_min:nx_max, ny_min:ny_max, 1:nz_max) !< 
        !! Array containing, for each flux defined in the *config.ini* file, 
        !! the value of the sub-fluxes (columns 1 to *nbfunc_max*) and the 
        !! value of the flux itself (column *nbfunc_max+1*). 
        integer, intent(in) :: ili  !< Row in the flux_val array
        integer, intent(in) :: jcol !< Column in the flux_var array
        integer, intent(in) :: nbfunc_loc  !< Number of biological sub-processes
        ! within the flux calculation
        integer :: iii, jj  ! Loop indexes
        integer :: i_st, j_st, k_st 
        integer :: file_flx_id  


    ! If the number of fluxes to save is not 0
    ! Then we loop over all the elements of the outflux array
    if (num_outflux > 0) then
        do iii = 1, num_outflux   ! Index of the flux to save

     ! Recovering the i, j, and k coordinates of the flux to save
             i_st = outflux(iii)%indi
             j_st = outflux(iii)%indj
             k_st = outflux(iii)%indk
             file_flx_id = 5000 + i_st + j_st + k_st

     ! Writing the time (no line break)
             write(file_flx_id, "(E12.5,$)") tps/3600.

     ! Writing the coordinates of the flux in the FLUX_VAL array
             write(file_flx_id, "(2I4, $)")  ili, jcol

     ! Looping over the 1-nbfunc_loc column of the ss_flux array
     ! in order to save the sub-fluxes
            do jj=1, nbfunc_loc
                 write(file_flx_id,"(E12.4,$)") ss_flux(jj,i_st,j_st,k_st)
            enddo

     ! if nbfunc_loc < nbfunc_max, we fill the remaining with 0
            do jj=nbfunc_loc+1, nbfunc_max
                write(file_flx_id,"(A10,$)") "NaN" 
            enddo

     ! Writing of the value of the flux itsel (adding linebreak)
            write(file_flx_id, "(E12.4)") ss_flux(nbfunc_loc+1,i_st,j_st,k_st)

    enddo
 end if

    end subroutine eco3m_write_fluxoutputs
#endif
!-------------------------------------------------------------------------------
    subroutine eco3m_write_2d_var(vvar)

    !> Writing of 2-D (time, depth) output variables. **Only for real(8) variables**
    !! \author Nicolas Barrier
    !! \date 2016-02-04
!-------------------------------------------------------------------------------
        use mod_eco3m, only:tps,nz_max
        use mod_eco3m_outputs, only:filecount
        implicit none 

        ! local variables
        real(8), intent(in) :: vvar(1:nz_max)
        integer :: k ! Loop index

        ! Writing of the time-array in hours (no line break)
        write(filecount, "(E12.5,$)") tps/3600. ! time in hours
        ! Writing of the ndepth-1 first value (no line break)
        do k = 1, nz_max -1
            write(filecount, "(E12.4,$)") vvar(k)
        end do

        ! Writing of the last value with a line break
        write(filecount, "(E12.4)") vvar(nz_max)

    end subroutine  eco3m_write_2d_var
!-------------------------------------------------------------------------------
    subroutine eco3m_write_2d_var_r4(vvar)

    !> Writing of 2D (time, depth) output variables. **Only for real(4) variables**
    !! \author Nicolas Barrier
    !! \date 2016-02-04
    !! \note This function must be used when saving variables that come from the TKE model
!-------------------------------------------------------------------------------
        use mod_eco3m, only:tps,nz_max
        use mod_eco3m_outputs, only:filecount
        implicit none 

        ! local variables
        real, intent(in) :: vvar(1:nz_max)
        integer :: k ! Loop index

        write(filecount, "(E12.4,$)") tps/3600.! time in hours
        do k = 1, nz_max -1
            write(filecount, "(E12.4,$)") vvar(k)
        end do
        write(filecount, "(E12.4)") vvar(nz_max)

    end subroutine  eco3m_write_2d_var_r4
!-------------------------------------------------------------------------------

    !> Writing of scalar (time) output variables. **Only for real(8) variables**
    !! \author nicolas barrier
    !! \date 2016-02-04
    subroutine eco3m_write_1d_var(vvar)
!-------------------------------------------------------------------------------
        use mod_eco3m, only:tps
        use mod_eco3m_outputs, only:filecount
        implicit none 

        real(8), intent(in) :: vvar
        write(filecount, "(E12.4,$)") tps/3600. ! time in hours
        write(filecount, "(E12.4)") vvar
    end subroutine  eco3m_write_1d_var
!-------------------------------------------------------------------------------
    subroutine eco3m_write_1d_var_r4(vvar)

    !> Writing of scalar (time) output variables. **Only for real(8) variables**
    !! \author Nicolas Barrier
    !! \date 2016-02-04
    !! \note This function must be used when saving variables that come from the TKE model
!-------------------------------------------------------------------------------
        use mod_eco3m, only:tps
        use mod_eco3m_outputs, only:filecount
        implicit none 

        real, intent(in) :: vvar
        write(filecount, "(E12.4,$)") tps/3600. ! time in hours
        write(filecount, "(E12.4)") vvar
    end subroutine  eco3m_write_1d_var_r4
!-------------------------------------------------------------------------------
    subroutine eco3m_write_budgets


    !> **Only if NCOUPL.** Writing of total matter, C, N and P budgets into the budget file.
    !! Concentration are summed over all space indexes (x, y, z) and over all variables.
    !! \author Melika Baklouti
    !! \author Nicolas Barrier
!-------------------------------------------------------------------------------
        use mod_eco3m, only : VAR,nbvar
        use mod_eco3m_outputs
        use mod_eco3m_files, only:file_budget_id

        ! Local variables
        real(8) :: som,somC, somN ,somP
        integer :: i

        ! Conservativity of total matter
        som=0.d0
        do i = 1, nbvar
            som = som + sum(VAR(i)%conc)
        enddo
#ifdef key_nemo_eco3m
        if(lk_mpp) call mpp_sum(som)
#endif

        ! Sum of carbon concentrations (carbon budget)
        somC = 0.d0
        do i = 1, nbvar
            if (var(i)%elmt == 'C') then
                somC = somC + sum(VAR(i)%conc)
            endif
        enddo
#ifdef key_nemo_eco3m
        if(lk_mpp) call mpp_sum(somC)
#endif

        ! Sum of nitrogen concentrations (nitrogen budget)
        somN = 0.d0
        do i =1, nbvar
            if (var(i)%elmt == 'N') then
                somN = somN + sum(VAR(i)%conc)
            endif
        enddo
#ifdef key_nemo_eco3m
        if(lk_mpp) call mpp_sum(somN)
#endif

        ! Sum of phosphorus concentrations (phosphorus budget)
        somP = 0.d0
        do i = 1, nbvar
            if (var(i)%elmt == 'P') then  
                somP = somP + sum(VAR(i)%conc)
            endif
        enddo
#ifdef key_nemo_eco3m
        if(lk_mpp) call mpp_sum(somP)
#endif

        write(file_budget_id,'(g12.3,1X,g12.3,1X,g12.3,1X,g12.3)') &
            tps/3600.d0, somC, somN, somP

    end subroutine eco3m_write_budgets
!-------------------------------------------------------------------------------
#ifdef CALC
#ifdef SAVE_FLUX

subroutine eco3m_write_fluxes
! Subroutine allowing to save the biogeochemical fluxes when the SAVE_FLUX key
! is activated
!-------------------------------------------------------------------------------
  use mod_eco3m_vartypes
  use mod_eco3m
  use mod_eco3m_fprocess
  implicit None
  real(8),Allocatable:: ss_flux(:,:,:,:)
  integer :: nbproc_loc 
 
  include "call_save_flux.inc"

end subroutine eco3m_write_fluxes

#endif
#endif
!-------------------------------------------------------------------------------

