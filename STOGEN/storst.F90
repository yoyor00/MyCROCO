MODULE storst
   !!======================================================================
   !!                       ***  MODULE  storst  ***
   !! Purpose        : read and write restart files for stochastic fields
   !!=====================================================================
   !!   sto_rst_read  : read restart file
   !!   sto_rst_write : write restart file
   !!----------------------------------------------------------------------
   USE netcdf
   USE stoarray         ! stochastic arrays to store in restart
   USE stowhite         ! white noise generator
   USE storng_kiss      ! KISS random number generator
   USE storng_ziggurat  ! Ziggurat algorithm to generate normal numbers
   USE stoexternal, only : wp, lc, jpi, jpj, jpk, cn_mem, narea, ln_ensemble

   IMPLICIT NONE
   PRIVATE

   PUBLIC sto_rst_read, sto_rst_write

CONTAINS

   SUBROUTINE sto_rst_read
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_rst_read  ***
      !!
      !! ** Purpose :   read stochastic parameters from restart file
      !!----------------------------------------------------------------------
      ! Variables to describe NetCDF file
      CHARACTER(len=lc) :: rstfile
      INTEGER :: ierr, idf, idv0d, idv2d, idv3d, idvseed
      CHARACTER(len=4) :: cn_area
      ! Variables to describe seed of random number generator
      INTEGER(KIND=8) :: ziseed8(4) ! RNG seeds in 64-bit integer type
      INTEGER(KIND=4) :: ziseed4(4) ! RNG seeds in 32-bit integer type
      REAL(KIND=8)    :: zrseed8(4) ! RNG seeds in 64-bit real type (with same bits to save in restart)
      REAL(KIND=4)    :: zrseed4(4) ! RNG seeds in 32-bit real type (with same bits to save in restart)

      ! Set name of restart file (function of processor index and ensemble member index)
      rstfile = cn_storst_in
      WRITE(cn_area,'(i4.4)') narea
      IF ( narea > 1   ) rstfile = rstfile(1:len(rstfile)) //  cn_area
      IF ( ln_ensemble ) rstfile = cn_mem // rstfile(1:len(rstfile))

      ! Open NetCDF file
      ierr = NF90_OPEN(rstfile,NF90_NOWRITE,idf)
      IF (ierr.NE.0) STOP 'Error opening input restart file'

      ! Get variable ids
      IF (jpsto0d>0) ierr = NF90_INQ_VARID(idf,'sto0d',idv0d)
      IF (ierr.NE.0) STOP 'Variable sto0d not found in restart file'
      IF (jpsto2d>0) ierr = NF90_INQ_VARID(idf,'sto2d',idv2d)
      IF (ierr.NE.0) STOP 'Variable sto2d not found in restart file'
      IF (jpsto3d>0) ierr = NF90_INQ_VARID(idf,'sto3d',idv3d)
      IF (ierr.NE.0) STOP 'Variable sto3d not found in restart file'
      ierr = NF90_INQ_VARID(idf,'seed',idvseed)
      IF (ierr.NE.0) STOP 'Variable seed not found in restart file'

      ! Read stochastic arrays in restart file
      IF (jpsto0d>0) ierr = NF90_GET_VAR(idf,idv0d,sto0d)
      IF (ierr.NE.0) STOP 'Error reading sto0d in restart file'
      IF (jpsto2d>0) ierr = NF90_GET_VAR(idf,idv2d,sto2d)
      IF (ierr.NE.0) STOP 'Error reading sto2d in restart file'
      IF (jpsto3d>0) ierr = NF90_GET_VAR(idf,idv3d,sto3d)
      IF (ierr.NE.0) STOP 'Error reading sto3d in restart file'

      ! Get state of random number generator from restart file
      ! and reinitialize random number generator from there
      IF (ln_rstseed) THEN
        SELECT CASE (c_rngtype)
        CASE('kiss64')
          ierr = NF90_GET_VAR(idf,idvseed,zrseed8)
          IF (ierr.NE.0) STOP 'Error reading seed in restart file'
          ziseed8 = TRANSFER( zrseed8 , ziseed8 )
          CALL kiss_seed( ziseed8(1) , ziseed8(2) , ziseed8(3) , ziseed8(4) )
        CASE('kiss32')
          ierr = NF90_GET_VAR(idf,idvseed,ziseed4)
          IF (ierr.NE.0) STOP 'Error reading seed in restart file'
          CALL kiss_seed( ziseed4(1) , ziseed4(2) , ziseed4(3) , ziseed4(4) )
        CASE('shr3')
          ierr = NF90_GET_VAR(idf,idvseed,ziseed4(1))
          IF (ierr.NE.0) STOP 'Error reading seed in restart file'
          CALL shr3_seed( ziseed4(1) )
        CASE DEFAULT
          STOP 'Bad type of random number generator in storst'
        END SELECT
      ENDIF

      ! Close NetCDF file
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) STOP 'Error closing input restart file'

   END SUBROUTINE sto_rst_read


   SUBROUTINE sto_rst_write
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_rst_read  ***
      !!
      !! ** Purpose :   write stochastic parameters in restart file
      !!----------------------------------------------------------------------
      ! Variables to describe NetCDF file
      CHARACTER(len=lc) :: rstfile
      LOGICAL :: filexists
      INTEGER :: nf90_type
      INTEGER :: ierr, idf, idx, idy, idz
      INTEGER :: idt0d, idt2d, idt3d, idn0d, idn2d, idn3d, idv0d, idv2d, idv3d
      INTEGER :: idseed, idvseed
      CHARACTER(len=4) :: cn_area
      ! Variables to describe seed of random number generator
      INTEGER :: seed_type, seed_size
      INTEGER(KIND=8) :: ziseed8(4) ! RNG seeds in 64-bit integer type
      INTEGER(KIND=4) :: ziseed4(4) ! RNG seeds in 32-bit integer type
      REAL(KIND=8)    :: zrseed8(4) ! RNG seeds in 64-bit real type (with same bits to save in restart)
      REAL(KIND=4)    :: zrseed4(4) ! RNG seeds in 32-bit real type (with same bits to save in restart)

      ! Set name of restart file (function of processor index and ensemble member index)
      rstfile = cn_storst_out
      WRITE(cn_area,'(i4.4)') narea
      IF ( narea > 1   ) rstfile = rstfile(1:len(rstfile)) //  cn_area
      IF ( ln_ensemble ) rstfile = cn_mem // rstfile(1:len(rstfile))

      ! Get type and state of random number generator to save in restart file
      SELECT CASE (c_rngtype)
      CASE('kiss64')
         seed_type = NF90_DOUBLE ; seed_size = 4
         CALL kiss_state( ziseed8(1) , ziseed8(2) , ziseed8(3) , ziseed8(4) )
         zrseed8 = TRANSFER( ziseed8 , zrseed8 )
      CASE('kiss32')
         seed_type = NF90_INT    ; seed_size = 4
         CALL kiss_state( ziseed4(1) , ziseed4(2) , ziseed4(3) , ziseed4(4) )
      CASE('shr3')
         seed_type = NF90_INT    ; seed_size = 1
         CALL shr3_state( ziseed4(1) )
      CASE DEFAULT
         STOP 'Bad type of random number generator in storst'
      END SELECT

      ! Set type of output arrays
      IF (wp == SELECTED_REAL_KIND(12,307) ) THEN
        nf90_type = NF90_DOUBLE
      ELSE
        nf90_type = NF90_FLOAT
      ENDIF

      ! Check if output restart file already exists
      INQUIRE (FILE=rstfile,EXIST=filexists)
      IF (filexists) THEN
        ! Open output NetCDF file
        ierr = NF90_OPEN(rstfile,NF90_WRITE,idf)
        IF (ierr.NE.0) STOP 'Error opening existing output restart file'
        ! Get variable ids
        IF (jpsto0d>0) ierr = NF90_INQ_VARID(idf,'sto0d',idv0d)
        IF (jpsto2d>0) ierr = NF90_INQ_VARID(idf,'sto2d',idv2d)
        IF (jpsto3d>0) ierr = NF90_INQ_VARID(idf,'sto3d',idv3d)
      ELSE
        ! Create output NetCDF file
        ierr = NF90_CREATE(rstfile,NF90_CLOBBER,idf)
        IF (ierr.NE.0) STOP 'Error creating output restart file'
        IF ( (jpsto2d>0) .OR. (jpsto3d>0) ) ierr = NF90_DEF_DIM(idf,'x',jpi,idx)
        IF (ierr.NE.0) STOP 'Error in writing restart file'
        IF ( (jpsto2d>0) .OR. (jpsto3d>0) ) ierr = NF90_DEF_DIM(idf,'y',jpj,idy)
        IF (ierr.NE.0) STOP 'Error in writing restart file'
        ! Create dimensions for stochastic arrays
        IF (jpsto3d>0) ierr = NF90_DEF_DIM(idf,'z',jpk,idz)
        IF (ierr.NE.0) STOP 'Error in writing restart file'
        IF (jpsto0d>0) ierr = NF90_DEF_DIM(idf,'tslices_0d',jpidx0d+jpidxsup0d,idt0d)
        IF (ierr.NE.0) STOP 'Error in writing restart file'
        IF (jpsto2d>0) ierr = NF90_DEF_DIM(idf,'tslices_2d',jpidx2d+jpidxsup2d,idt2d)
        IF (ierr.NE.0) STOP 'Error in writing restart file'
        IF (jpsto3d>0) ierr = NF90_DEF_DIM(idf,'tslices_3d',jpidx3d+jpidxsup3d,idt3d)
        IF (ierr.NE.0) STOP 'Error in writing restart file'
        IF (jpsto0d>0) ierr = NF90_DEF_DIM(idf,'nsto_0d',jpsto0d,idn0d)
        IF (ierr.NE.0) STOP 'Error in writing restart file'
        IF (jpsto2d>0) ierr = NF90_DEF_DIM(idf,'nsto_2d',jpsto2d,idn2d)
        IF (ierr.NE.0) STOP 'Error in writing restart file'
        IF (jpsto3d>0) ierr = NF90_DEF_DIM(idf,'nsto_3d',jpsto3d,idn3d)
        IF (ierr.NE.0) STOP 'Error in writing restart file'
        ! Create dimension for state of random number generator
        ierr = NF90_DEF_DIM(idf,'nseed',seed_size,idseed)
        IF (ierr.NE.0) STOP 'Error in writing restart file'
        ! Create variables
        IF (jpsto0d>0) ierr = NF90_DEF_VAR(idf,'sto0d',nf90_type,(/idt0d,idn0d/),            idv0d)
        IF (ierr.NE.0) STOP 'Error in writing restart file'
        IF (jpsto2d>0) ierr = NF90_DEF_VAR(idf,'sto2d',nf90_type,(/idx,idy,idt2d,idn2d/),    idv2d)
        IF (ierr.NE.0) STOP 'Error in writing restart file'
        IF (jpsto3d>0) ierr = NF90_DEF_VAR(idf,'sto3d',nf90_type,(/idx,idy,idz,idt3d,idn3d/),idv3d)
        IF (ierr.NE.0) STOP 'Error in writing restart file'
        ! Create variable to store state of random number generator
        ierr = NF90_DEF_VAR(idf,'seed',seed_type,(/idseed/),idvseed)
        IF (ierr.NE.0) STOP 'Error in writing restart file'
        ierr = NF90_ENDDEF(idf)
        IF (ierr.NE.0) STOP 'Error in writing restart file'
      ENDIF

      ! Store stochastic arrays in restart file
      IF (jpsto0d>0) ierr = NF90_PUT_VAR(idf,idv0d,sto0d)
      IF (ierr.NE.0) STOP 'Error writing sto0d in restart file'
      IF (jpsto2d>0) ierr = NF90_PUT_VAR(idf,idv2d,sto2d)
      IF (ierr.NE.0) STOP 'Error writing sto2d in restart file'
      IF (jpsto3d>0) ierr = NF90_PUT_VAR(idf,idv3d,sto3d)
      IF (ierr.NE.0) STOP 'Error writing sto3d in restart file'

      ! Store state of random number generator
      SELECT CASE (c_rngtype)
      CASE('kiss64')
         ierr = NF90_PUT_VAR(idf,idvseed,zrseed8)
      CASE('kiss32')
         ierr = NF90_PUT_VAR(idf,idvseed,ziseed4)
      CASE('shr3')
         ierr = NF90_PUT_VAR(idf,idvseed,ziseed4(1))
      CASE DEFAULT
         STOP 'Bad type of random number generator in storst'
      END SELECT
      IF (ierr.NE.0) STOP 'Error writing rng state in restart file'

      ! Close NetCDF file
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) STOP 'Error closing output restart file'

   END SUBROUTINE sto_rst_write

END MODULE storst
