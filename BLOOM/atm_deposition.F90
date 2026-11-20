
#include "cppdefs.h"

MODULE atm_deposition
#if defined BLOOM

      USE croco_namelist
      USE comsubstance
      USE comBIOLink , ONLY: path_dust,dust_filename_prefix,dust_ts, &
              dust_vname, mod_vars, conv_fac, nb_dust_vars, &
              iv_nutr_NO3, iv_nutr_NH4, iv_nutr_PO4

# if defined XIOS
      USE iom
# endif

#  if defined MPI
      use scalars, ONLY: mynode, time
#  endif 

        IMPLICIT NONE

# include "param.h"
# include "netcdf.inc"

        

        REAL(KIND=rsh), DIMENSION(:,:,:,:), ALLOCATABLE ::  dust_N
        CHARACTER(len=19), EXTERNAL          :: tool_sectodat


        contains
                SUBROUTINE atm_deposition_init()


                        INTEGER :: iyear, var_id, ierr, ncid_dust, time_id
                        INTEGER :: start_year, start_month, end_year, end_month
                        INTEGER :: nrec_dust_r, nrec_dust, nrec, indx, i_rec, lstr, max_rec
                        INTEGER :: nf_fread, lenstr, dust_varid,checkdims
                        CHARACTER(len=200) :: dust_filename
                        CHARACTER(len=4) :: iyear_str
                        !CHARACTER(len=1) :: dust_ts

                        REAL(KIND=rsh), DIMENSION(:,:,:), ALLOCATABLE ::  dustmp



                        ! Table allocation


                        read(start_date(1:4), '(I4)') start_year
                        read(start_date(6:7), '(I2)') start_month
                        read(end_date(1:4)  , '(I4)') end_year
                        read(end_date(6:7)  , '(I2)') end_month

                        ! the number of time values in each file
                        IF(dust_ts == "M") THEN
                                nrec_dust_r = (end_year - start_year + 1) * 12

                                ALLOCATE(dustmp(GLOBAL_2D_ARRAY,12))

                        ELSEIF (dust_ts == "D") THEN
                                nrec_dust_r = 0
                                DO iyear = start_year, end_year
                                        nrec_dust_r = nrec_dust_r + 365
                                        IF(MOD(iyear,4) ==0) nrec_dust_r = nrec_dust_r + 1
                                ENDDO

                                ALLOCATE(dustmp(GLOBAL_2D_ARRAY, 366))
                        ELSE
                                MPI_master_only write(stdout,*) 'ERROR: Invalid value for `dust_ts` in the namelist. Expected values are `D` (daily) or `M` (monthly).'
                                error stop

                        ENDIF


                        ! Allocation 

                        ALLOCATE( dust_N(GLOBAL_2D_ARRAY,nrec_dust_r,nb_dust_vars))


                        ! Initialization:
                        dust_N(:,:,:,:) = 0.

                        nrec = 0


                        DO iyear = start_year, end_year
                                
                                write(iyear_str,'(I0)') iyear

                                ! Get the file name

                                dust_filename = TRIM(path_dust)//"/"//TRIM(dust_filename_prefix)//"_" //TRIM(iyear_str)//".nc"

                                ! Open the file
                                lstr = lenstr(dust_filename)

                                ierr=nf_open(dust_filename(1:lstr), nf_nowrite, ncid_dust)

                                IF (ierr /= nf_noerr) THEN
                                        MPI_master_only write(stdout,*) 'ERROR: Unable to open NetCDF file: ', dust_filename(1:lstr)
                                        MPI_master_only write(stdout,*) 'NetCDF error: ', nf_strerror(ierr)
                                        error stop
                                ELSE
                                        MPI_master_only write(stdout,*) 'Open dust file to read the atmospheric deposition :', dust_filename(1:lstr)
                                ENDIF

                                ! Check if the file dimensions are correct

                                ierr = checkdims(ncid_dust, dust_filename, lstr, max_rec)

                                ! Extract the number of records

                                ierr = nf_inq_dimid(ncid_dust,"time_counter",time_id)
                                IF (ierr /= nf_noerr) THEN
                                        MPI_master_only write(stdout,*) 'ERROR: Unable to find dimension "time_counter" in file: ', &
                                                                         dust_filename(1:lstr)
                                        MPI_master_only write(stdout,*) 'NetCDF error: ', nf_strerror(ierr)
                                        error stop
                                ENDIF

                                ierr = nf_inq_dimlen(ncid_dust,time_id,nrec_dust)
                                IF (ierr /= nf_noerr) THEN
                                        MPI_master_only write(stdout,*) 'ERROR: Unable to get the length of dimension "time_counter" in file: ', &
                                                                        dust_filename(1:lstr)
                                        MPI_master_only write(stdout,*) 'NetCDF error: ', nf_strerror(ierr)
                                        error stop
                                ELSEIF (nrec_dust /= nrec_dust_r) THEN
                                       MPI_master_only write(stdout,*) 'ERROR: The number of time records in file : ', &
                                               dust_filename(1:lstr), ' does not match the expected value: ' , nrec_dust_r

                                      error stop

                                ENDIF



                                DO var_id = 1,nb_dust_vars 

                                        dustmp(:,:,:) = 0               

                                        MPI_master_only write(stdout,*) 'Reading variable :', TRIM(dust_vname(var_id))
                                       
                                       ! Extract the variable index

                                        ierr=nf_inq_varid(ncid_dust, TRIM(dust_vname(var_id)), dust_varid)
                                        IF (ierr /= nf_noerr) THEN
                                                MPI_master_only write(stdout,*) 'ERROR: Unable to find variable "', TRIM(dust_vname(var_id)), &
                                                        '" in NetCDF file: ', dust_filename(1:lstr)
                                                MPI_master_only write(stdout,*) 'NetCDF error: ', nf_strerror(ierr)
                                                error stop
                                        ENDIF


                                        DO indx = 1,nrec_dust
                                           i_rec = nrec + indx
                                           
                                           ierr=nf_fread (dustmp(START_2D_ARRAY,indx), ncid_dust, dust_varid, indx, 0) 

                                           IF (ierr /= NF_NOERR) THEN
                                                   MPI_master_only write(stdout,*) 'ERROR: Unable to read data for variable "', &
                                                           TRIM(dust_vname(var_id)), '" from NetCDF file: ', &
                                                           dust_filename(1:lstr)
                                                   MPI_master_only write(stdout,*) 'Record index: ', indx
                                                   MPI_master_only write(stdout,*) 'NetCDF error: ', nf_strerror(ierr)
                                                   error stop
                                           END IF


                                        ENDDO

                                         dust_N(:,:,nrec+1:nrec+nrec_dust,var_id) = dustmp(:,:,1:nrec_dust)

                                 ENDDO

                                 nrec = nrec + nrec_dust


                        END DO ! iyear


                        DEALLOCATE(dustmp)


                END SUBROUTINE atm_deposition_init


                SUBROUTINE include_atm_deposition(i,j,epn,dc,flx_dep_N)

                        INTEGER, INTENT(IN) :: i,j
                        REAL(KIND=rsh)      :: epn
                        REAL(KIND=rsh),DIMENSION(nv_adv), intent(inout)  :: dc
                        REAL(KIND=rsh), DIMENSION(:),ALLOCATABLE, INTENT(INOUT) :: flx_dep_N

                        INTEGER,DIMENSION(12),PARAMETER :: nb_days_month=(/31,28,31,30,31,30,31,31,30,31,30,31/)
                        INTEGER, DIMENSION(3)           :: nut_idx
                        CHARACTER(len=25), PARAMETER :: nutrient_names(3) = &
                                [CHARACTER(len=25) :: 'nitrate','ammonium','dissolved_phosphate']

                        INTEGER :: start_year, current_year, current_month, mth, nb_sec_month, t_pos
                        INTEGER :: i_var, i_nut, iv_nut
                        CHARACTER(len=19) :: current_date


                        nut_idx=(/iv_nutr_NO3, iv_nutr_NH4, iv_nutr_PO4/)


                        read(start_date(1:4), '(I4)') start_year

                        current_date = TRIM(tool_sectodat(time))

                        read(current_date(1:4), '(I4)') current_year
                        read(current_date(6:7), '(I2)') current_month

                        IF(dust_ts == "M") THEN
                                mth = nb_days_month(current_month)
                                IF((current_month == 2) .AND. (MOD(current_year,4) == 0)) &
                                        mth = nb_days_month(current_month) + 1
                        ELSEIF(dust_ts == "D") THEN
                                mth = 1
                        ENDIF
       

                        ! the position to be read in the 3rd dimension of the table related to time

                        t_pos = (current_year - start_year) * 12 + current_month

                        ! the position to be read in the 4th dimension of the table related to variable
                        DO i_var=1,nb_dust_vars
    
                            ! Flux in mmole m-2 d-1

                             flx_dep_N(i_var) = dust_N(i,j,t_pos,i_var) * conv_fac(i_var) / mth
                             DO i_nut = 1, size(nutrient_names)
                                IF(TRIM(mod_vars(i_var)) == TRIM(nutrient_names(i_nut))) THEN
                                        iv_nut = nut_idx(i_nut)
                                        exit
                                END IF
                             END DO

                             ! update the dc of the corresponding state variable
                             dc(iv_nut) = dc(iv_nut) + flx_dep_N(i_var) /epn ! epn = Hz(i,j,N) the height of the shallowest layer

                        END DO



                END SUBROUTINE include_atm_deposition

#endif
END MODULE atm_deposition
