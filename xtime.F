c     timer
      function xtime()
      implicit none
      real xtime
      include 'mpif.h'

      xtime = MPI_wtime()

      end function xtime
