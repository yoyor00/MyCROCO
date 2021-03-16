C     some dummy procedures for tapenade

      subroutine output
      end subroutine output

      subroutine diag (tile)
      integer tile
      end subroutine diag
      
      subroutine diag_tile(Istr,Iend,Jstr,Jend, ke2d,pe2d)
      integer Istr,Iend,Jstr,Jend
      real ke2d(1,1,1), pe2d(1,1,1)
      end subroutine diag_tile

      subroutine get_ssh
      end subroutine get_ssh

      subroutine set_ssh
      end subroutine set_ssh

c      subroutine get_vbc
c      end subroutine get_vbc

c      subroutine set_vbc
c      end subroutine set_vbc

      subroutine clm_tides(tile)
      integer tile
      
      end subroutine clm_tides

      subroutine pushpointer8(a)

      end subroutine

      subroutine cpu_time()

      end subroutine
      

      subroutine trim(str)
      character*30 str
      end subroutine

      subroutine read_inp()

      end subroutine
