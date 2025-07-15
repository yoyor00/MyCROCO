!
! Helper header for latency hiding.
!
! It extends the iteration ranges set up by `compute_auxiliary_bounds.h`
! to include also halo layers used for latency hiding / communication avoiding /
! overlapping Schwarz.
!

#ifdef EW_PERIODIC

      imin = imin - MPI_LAT_HID_2D_ADD_LAYERS_AUX
      IstrR = IstrR - MPI_LAT_HID_2D_ADD_LAYERS_AUX
      ! IstrU was defined by precompiler, points to Istr
      !IstrU = IstrU - MPI_LAT_HID_2D_ADD_LAYERS_AUX
      !IstrV = IstrV - MPI_LAT_HID_2D_ADD_LAYERS_AUX

      imax = imax + MPI_LAT_HID_2D_ADD_LAYERS_AUX
      IendR = IendR + MPI_LAT_HID_2D_ADD_LAYERS_AUX
      !IendU = IendU + MPI_LAT_HID_2D_ADD_LAYERS_AUX
      !IendV = IendV + MPI_LAT_HID_2D_ADD_LAYERS_AUX

      Istr = Istr - MPI_LAT_HID_2D_ADD_LAYERS_AUX
      Iend = Iend + MPI_LAT_HID_2D_ADD_LAYERS_AUX

#else

      if (.not. (WESTERN_EDGE)) then
        imin = imin - MPI_LAT_HID_2D_ADD_LAYERS_AUX
        IstrR = IstrR - MPI_LAT_HID_2D_ADD_LAYERS_AUX
        IstrU = IstrU - MPI_LAT_HID_2D_ADD_LAYERS_AUX
        !IstrV = IstrV - MPI_LAT_HID_2D_ADD_LAYERS_AUX

        Istr = Istr - MPI_LAT_HID_2D_ADD_LAYERS_AUX
      endif

      if (.not. (EASTERN_EDGE)) then
        imax = imax + MPI_LAT_HID_2D_ADD_LAYERS_AUX
        IendR = IendR + MPI_LAT_HID_2D_ADD_LAYERS_AUX
        !IendU = IendU + MPI_LAT_HID_2D_ADD_LAYERS_AUX
        !IendV = IendV + MPI_LAT_HID_2D_ADD_LAYERS_AUX

        Iend = Iend + MPI_LAT_HID_2D_ADD_LAYERS_AUX
      endif

#endif


#ifdef NS_PERIODIC

      jmin = jmin - MPI_LAT_HID_2D_ADD_LAYERS_AUX
      JstrR = JstrR - MPI_LAT_HID_2D_ADD_LAYERS_AUX
      !JstrU = JstrU - MPI_LAT_HID_2D_ADD_LAYERS_AUX
      ! JstrV was defined by precompiler, points to Jstr
      !JstrV = JstrV - MPI_LAT_HID_2D_ADD_LAYERS_AUX

      jmax = jmax + MPI_LAT_HID_2D_ADD_LAYERS_AUX
      JendR = JendR + MPI_LAT_HID_2D_ADD_LAYERS_AUX
      !JendU = JendU + MPI_LAT_HID_2D_ADD_LAYERS_AUX
      !JendV = JendV + MPI_LAT_HID_2D_ADD_LAYERS_AUX

      Jstr = Jstr - MPI_LAT_HID_2D_ADD_LAYERS_AUX
      Jend = Jend + MPI_LAT_HID_2D_ADD_LAYERS_AUX

#else

      if (.not. (SOUTHERN_EDGE)) then
        jmin = jmin - MPI_LAT_HID_2D_ADD_LAYERS_AUX
        JstrR = JstrR - MPI_LAT_HID_2D_ADD_LAYERS_AUX
        JstrV = JstrV - MPI_LAT_HID_2D_ADD_LAYERS_AUX

        Jstr = Jstr - MPI_LAT_HID_2D_ADD_LAYERS_AUX
      endif

      if (.not. (NORTHERN_EDGE)) then
        jmax = jmax + MPI_LAT_HID_2D_ADD_LAYERS_AUX
        JendR = JendR + MPI_LAT_HID_2D_ADD_LAYERS_AUX
        !JendU = JendU + MPI_LAT_HID_2D_ADD_LAYERS_AUX
        !JendV = JendV + MPI_LAT_HID_2D_ADD_LAYERS_AUX

        Jend = Jend + MPI_LAT_HID_2D_ADD_LAYERS_AUX
      endif

#endif
