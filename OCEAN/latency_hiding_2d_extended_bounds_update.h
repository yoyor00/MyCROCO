!
! Helper header for latency hiding.
!
! It extends the iteration ranges set up by `compute_extended_bounds.h`
!

#ifdef EW_PERIODIC

      IstrR = IstrR - MPI_LAT_HID_2D_ADD_LAYERS
      IendR = IendR + MPI_LAT_HID_2D_ADD_LAYERS

      Istr = Istr - MPI_LAT_HID_2D_ADD_LAYERS
      Iend = Iend + MPI_LAT_HID_2D_ADD_LAYERS

#else

      if (.not. (WESTERN_EDGE)) then
        IstrR = IstrR - MPI_LAT_HID_2D_ADD_LAYERS
        Istr = Istr - MPI_LAT_HID_2D_ADD_LAYERS
      endif

      if (.not. (EASTERN_EDGE)) then
        IendR = IendR + MPI_LAT_HID_2D_ADD_LAYERS
        Iend = Iend + MPI_LAT_HID_2D_ADD_LAYERS
      endif

#endif


#ifdef NS_PERIODIC

      JstrR = JstrR - MPI_LAT_HID_2D_ADD_LAYERS
      JendR = JendR + MPI_LAT_HID_2D_ADD_LAYERS

      Jstr = Jstr - MPI_LAT_HID_2D_ADD_LAYERS
      Jend = Jend + MPI_LAT_HID_2D_ADD_LAYERS

#else

      if (.not. (SOUTHERN_EDGE)) then
        JstrR = JstrR - MPI_LAT_HID_2D_ADD_LAYERS
        Jstr = Jstr - MPI_LAT_HID_2D_ADD_LAYERS
      endif

      if (.not. (NORTHERN_EDGE)) then
        JendR = JendR + MPI_LAT_HID_2D_ADD_LAYERS
        Jend = Jend + MPI_LAT_HID_2D_ADD_LAYERS
      endif

#endif
