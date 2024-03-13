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
        subroutine twin_check_float(value, equation, equation_size, location_id, source_line) &
            & bind(C, name="twin_check_float")
            use iso_c_binding, only: c_char, c_size_t, c_float, c_int, c_int64_t
            real(c_float), value, intent(in) :: value
            character(kind=c_char), dimension(*) :: equation
            integer(c_size_t), value :: equation_size
            integer (c_int64_t), value :: location_id
            integer (c_int), value :: source_line
        end subroutine twin_check_float

        subroutine twin_check_float_fixable(value, equation, equation_size, location_id, source_line) &
            & bind(C, name="twin_check_float_fixable")
            use iso_c_binding, only: c_char, c_size_t, c_float, c_int, c_int64_t
            real(c_float), intent(inout) :: value
            character(kind=c_char), dimension(*) :: equation
            integer(c_size_t), value :: equation_size
            integer (c_int64_t), value :: location_id
            integer (c_int), value :: source_line
        end subroutine twin_check_float_fixable

        subroutine twin_check_double(value, equation, equation_size, location_id, source_line) &
            & bind(C, name="twin_check_double")
            use iso_c_binding, only: c_char, c_size_t, c_double, c_int64_t, c_int
            real(c_double), value, intent(in) :: value
            character(kind=c_char), dimension(*) :: equation
            integer(c_size_t), value :: equation_size
            integer (c_int64_t), value :: location_id
            integer (c_int), value :: source_line
        end subroutine twin_check_double

        subroutine twin_check_double_fixable(value, equation, equation_size, location_id, source_line) &
            & bind(C, name="twin_check_double_fixable")
            use iso_c_binding, only: c_char, c_size_t, c_double, c_int, c_int64_t
            real(c_double), intent(inout) :: value
            character(kind=c_char), dimension(*) :: equation
            integer(c_size_t), value :: equation_size
            integer (c_int64_t), value :: location_id
            integer (c_int), value :: source_line
        end subroutine twin_check_double_fixable

        subroutine twin_check_integer(value, equation, equation_size, location_id, source_line) &
            & bind(C, name="twin_check_int")
            use iso_c_binding, only: c_char, c_int, c_size_t, c_double, c_int, c_int64_t
            integer(c_int), value, intent(in) :: value
            character(kind=c_char), dimension(*) :: equation
            integer(c_size_t), value :: equation_size
            integer (c_int64_t), value :: location_id
            integer (c_int), value :: source_line
        end subroutine twin_check_integer

        subroutine twin_check_integer_fixable(value, equation, equation_size, location_id, source_line) &
            & bind(C, name="twin_check_int_fixable")
            use iso_c_binding, only: c_char, c_int, c_size_t, c_double, c_int, c_int64_t
            integer(c_int), intent(inout) :: value
            character(kind=c_char), dimension(*) :: equation
            integer(c_size_t), value :: equation_size
            integer (c_int64_t), value :: location_id
            integer (c_int), value :: source_line
        end subroutine twin_check_integer_fixable

        subroutine twin_check_bool(value, equation, equation_size, location_id, source_line) &
            & bind(C, name="twin_check_bool")
            use iso_c_binding, only: c_char, c_size_t, c_bool, c_int, c_int64_t
            integer(c_int), value, intent(in) :: value
            character(kind=c_char), dimension(*) :: equation
            integer(c_size_t), value :: equation_size
            integer (c_int64_t), value :: location_id
            integer (c_int), value :: source_line
        end subroutine twin_check_bool

        subroutine twin_check_bool_fixable(value, equation, equation_size, location_id, source_line) &
            & bind(C, name="twin_check_bool_fixable")
            use iso_c_binding, only: c_char, c_size_t, c_bool, c_int, c_int64_t
            integer(c_int), intent(inout) :: value
            character(kind=c_char), dimension(*) :: equation
            integer(c_size_t), value :: equation_size
            integer (c_int64_t), value :: location_id
            integer (c_int), value :: source_line
        end subroutine twin_check_bool_fixable

        subroutine twin_register_site(id, source_file, source_file_size) &
            & bind(C, name="twin_register_site")
            use iso_c_binding, only: c_char, c_int, c_size_t, c_double, c_int, c_int64_t
            integer(c_int64_t), value :: id
            character(kind=c_char), dimension(*) :: source_file
            integer(c_size_t), value :: source_file_size
        end subroutine twin_register_site
    end interface
end module
