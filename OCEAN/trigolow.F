#include "cppdefs.h"
#if defined RVTK_DEBUG             
      module my_functions_1
      contains     
      elemental function mylog(x)
!$acc routine(mylog) seq
      real,intent(in) :: x
      real :: mylog
      mylog = real(log(X),kind=4)
      end function mylog

      elemental function mytanh(x)
!$acc routine(mytanh) seq
      real,intent(in) :: x
      real :: mytanh
      mytanh = real(tanh(X),kind=4)
      end function mytanh 

      elemental function myexp(x)
!$acc routine(myexp) seq
      real,intent(in) :: x
      real :: myexp
      myexp = real(exp(X),kind=4)
      end function myexp

      elemental function mycos(x)
!$acc routine(mycos) seq
      real,intent(in) :: x
      real :: mycos
      mycos = real(cos(X),kind=4)
      end function mycos

      elemental function mysin(x)
!$acc routine(mysin) seq
      real,intent(in) :: x
      real :: mysin
      mysin = real(sin(X),kind=4)
      end function mysin

      elemental function mycosh(x)
!$acc routine(mycosh) seq
      real,intent(in) :: x
      real :: mycosh
      mycosh = real(cosh(X),kind=4)
      end function mycosh

      elemental function mysinh(x)
!$acc routine(mysinh) seq
      real,intent(in) :: x
      real :: mysinh
      mysinh = real(sinh(X),kind=4)
      end function mysinh

      elemental function myatan(x)
!$acc routine(myatan) seq
      real,intent(in) :: x
      real :: myatan
      myatan = real(atan(X),kind=4)
      end function myatan

      elemental function mysqrt(x)
!$acc routine(mysqrt) seq
      real,intent(in) :: x
      real :: mysqrt
      mysqrt = real(sqrt(X),kind=4)
      end function mysqrt


      end module my_functions_1 

      module my_functions
      USE my_functions_1, only : log=>mylog
      USE my_functions_1, only : tanh=>mytanh
      USE my_functions_1, only : exp=>myexp
      USE my_functions_1, only : cos=>mycos
      USE my_functions_1, only : sin=>mysin
      USE my_functions_1, only : cosh=>mycosh
      USE my_functions_1, only : sinh=>mysinh
      USE my_functions_1, only : atan=>myatan
      USE my_functions_1, only : sqrt=>mysqrt
      end module my_functions
#else
      module my_functions_empty
      end module my_functions_empty    
#endif
