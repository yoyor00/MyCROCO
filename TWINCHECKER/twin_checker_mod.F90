!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            PROJECT  : psyclone-roofline            !
!            VERSION  : 0.0.0                        !
!            DATE     : 07/2023                      !
!            AUTHOR   : Valat Sébastien              !
!            LICENSE  : CeCILL-C                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module twin_checker
    use, intrinsic :: iso_c_binding
    implicit none
    interface
        subroutine twin_check_float(value) bind(C, name="twin_check_float")
            use iso_c_binding, only: c_char, c_size_t, c_float
            real(c_float), value, intent(in) :: value
        end subroutine twin_check_float

        subroutine twin_check_float_fixable(value) bind(C, name="twin_check_float_fixable")
            use iso_c_binding, only: c_char, c_size_t, c_float
            real(c_float), intent(inout) :: value
        end subroutine twin_check_float_fixable

        subroutine twin_check_double(value) bind(C, name="twin_check_double")
            use iso_c_binding, only: c_char, c_size_t, c_double
            real(c_double), value, intent(in) :: value
        end subroutine twin_check_double

        subroutine twin_check_double_fixable(value) bind(C, name="twin_check_double_fixable")
            use iso_c_binding, only: c_char, c_size_t, c_double
            real(c_double), intent(inout) :: value
        end subroutine twin_check_double_fixable

        subroutine twin_check_integer(value) bind(C, name="twin_check_int")
            use iso_c_binding, only: c_char, c_int, c_size_t, c_double
            integer(c_int), value, intent(in) :: value
        end subroutine twin_check_integer

        subroutine twin_check_integer_fixable(value) bind(C, name="twin_check_int_fixable")
            use iso_c_binding, only: c_char, c_int, c_size_t, c_double
            integer(c_int), intent(inout) :: value
        end subroutine twin_check_integer_fixable
    end interface
end module
