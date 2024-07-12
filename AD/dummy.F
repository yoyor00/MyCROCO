C     some dummy procedures for tapenade

      subroutine output
      end subroutine output

      subroutine diag(tile)
      end subroutine diag

      subroutine diag_tile(Istr,Iend,Jstr,Jend, ke2d,pe2d)
      integer Istr,Iend,Jstr,Jend
      real ke2d(1,1,1), pe2d(1,1,1)
      end subroutine diag_tile

      subroutine get_ssh
      end subroutine get_ssh

      subroutine set_ssh
      end subroutine set_ssh

      subroutine get_vbc
      end subroutine get_vbc

      subroutine set_vbc
      end subroutine set_vbc

      subroutine get_smflux
      end subroutine get_smflux

      subroutine clm_tides(tile)
      integer tile
      
      end subroutine clm_tides
