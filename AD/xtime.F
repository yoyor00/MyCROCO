c     timer
      function xtime()
      implicit none
      real xtime
#ifdef MPI
      include 'mpif.h'

      xtime = MPI_wtime()
#else
      call cpu_time(xtime)
#endif
      end function xtime
