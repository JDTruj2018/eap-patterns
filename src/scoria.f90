!!  =======================================================================
!!  © (or copyright) 2022. Triad National Security, LLC. All rights
!!  reserved.  This program was produced under U.S. Government contract
!!  89233218CNA000001 for Los Alamos National Laboratory (LANL), which is
!!  operated by Triad National Security, LLC for the U.S.  Department of
!!  Energy/National Nuclear Security Administration. All rights in the
!!  program are reserved by Triad National Security, LLC, and the
!!  U.S. Department of Energy/National Nuclear Security
!!  Administration. The Government is granted for itself and others acting
!!  on its behalf a nonexclusive, paid-up, irrevocable worldwide license
!!  in this material to reproduce, prepare derivative works, distribute
!!  copies to the public, perform publicly and display publicly, and to
!!  permit others to do so.
!!
!!  See LICENSE file for details
!!  =======================================================================

module scoria
    use iso_c_binding
    use iso_fortran_env, only: REAL64, INT64, INT32, INT8
    implicit none

    interface
        subroutine print_tsum(tsum, ndim) BIND(C, name="print_tsum")
            use iso_c_binding
            implicit none
            real(c_double), value, intent(in) :: tsum
            integer(c_int32_t), value, intent(in) :: ndim
        end subroutine print_tsum

        subroutine print_face_num(fnum, len) BIND(C, name="print_face_num")
            use iso_c_binding
            implicit none
            integer(c_int32_t), value, intent(in) :: len
            integer(c_int32_t), intent(in) :: fnum(0:len)
        end subroutine print_face_num

        subroutine read_single_thread_0_wrapper(res, buffer, n, intrin) &
                   BIND(C, name="read_single_thread_0_wrapper")
            use iso_c_binding
            implicit none
            
            integer(c_size_t), value, intent(in) :: n
            real(c_double), intent(in) :: res(0:n)
            real(c_double), intent(in) :: buffer(0:n)
            integer(c_int), value, intent(in) :: intrin
        end subroutine read_single_thread_0_wrapper

        subroutine read_single_thread_1_wrapper(res, buffer, n, ind1, n1, intrin) &
                   BIND(C, name="read_single_thread_1_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            real(c_double), intent(in) :: res(0:n1)
            real(c_double), intent(in) :: buffer(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n)
            integer(c_int), value, intent(in) :: intrin
        end subroutine read_single_thread_1_wrapper

        subroutine read_single_thread_2_wrapper(res, buffer, n, ind1, n1, ind2, n2, intrin) &
                   BIND(C, name="read_single_thread_2_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            integer(c_size_t), value, intent(in) :: n2
            real(c_double), intent(in) :: res(0:n2)
            real(c_double), intent(in) :: buffer(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n)
            integer(c_size_t), intent(in) :: ind2(0:n1)
            integer(c_int), value, intent(in) :: intrin
        end subroutine read_single_thread_2_wrapper

        subroutine read_multi_thread_0_wrapper(res, buffer, n, nthreads, intrin) &
                   BIND(C, name="read_multi_thread_0_wrapper")
            use iso_c_binding
            implicit none
            
            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: nthreads
            real(c_double), intent(in) :: res(0:n)
            real(c_double), intent(in) :: buffer(0:n)
            integer(c_int), value, intent(in) :: intrin 
        end subroutine read_multi_thread_0_wrapper

        subroutine read_multi_thread_1_wrapper(res, buffer, n, ind1, n1, nthreads, intrin) &
                   BIND(C, name="read_multi_thread_1_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            integer(c_size_t), value, intent(in) :: nthreads
            real(c_double), intent(in) :: res(0:n1)
            real(c_double), intent(in) :: buffer(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n)
            integer(c_int), value, intent(in) :: intrin 
        end subroutine read_multi_thread_1_wrapper

        subroutine read_multi_thread_2_wrapper(res, buffer, n, ind1, n1, ind2, n2, nthreads, intrin) &
                   BIND(C, name="read_multi_thread_2_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            integer(c_size_t), value, intent(in) :: n2
            integer(c_size_t), value, intent(in) :: nthreads
            real(c_double), intent(in) :: res(0:n2)
            real(c_double), intent(in) :: buffer(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n)
            integer(c_size_t), intent(in) :: ind2(0:n1)
            integer(c_int), value, intent(in) :: intrin 
        end subroutine read_multi_thread_2_wrapper


        subroutine write_single_thread_0_wrapper(buffer, input, n, intrin) &
                   BIND(C, name="write_single_thread_0_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            real(c_double), intent(in) :: buffer(0:n)
            real(c_double), intent(in) :: input(0:n)
            integer(c_int), value, intent(in) :: intrin
        end subroutine write_single_thread_0_wrapper

        subroutine write_single_thread_1_wrapper(buffer, input, n, ind1, n1, intrin) &
                   BIND(C, name="write_single_thread_1_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            real(c_double), intent(in) :: buffer(0:n1)
            real(c_double), intent(in) :: input(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n)
            integer(c_int), value, intent(in) :: intrin 
        end subroutine write_single_thread_1_wrapper

        subroutine write_single_thread_2_wrapper(buffer, input, n, ind1, n1, ind2, n2, intrin) &
                   BIND(C, name="write_single_thread_2_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            integer(c_size_t), value, intent(in) :: n2
            real(c_double), intent(in) :: buffer(0:n2)
            real(c_double), intent(in) :: input(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n)
            integer(c_size_t), intent(in) :: ind2(0:n1)
            integer(c_int), value, intent(in) :: intrin
        end subroutine write_single_thread_2_wrapper

        subroutine write_multi_thread_0_wrapper(buffer, input, n, nthreads, intrin) &
                   BIND(C, name="write_multi_thread_0_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: nthreads
            real(c_double), intent(in) :: buffer(0:n)
            real(c_double), intent(in) :: input(0:n)
            integer(c_int), value, intent(in) :: intrin
        end subroutine write_multi_thread_0_wrapper

        subroutine write_multi_thread_1_wrapper(buffer, input, n, ind1, n1, nthreads, intrin) &
                   BIND(C, name="write_multi_thread_1_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            integer(c_size_t), value, intent(in) :: nthreads
            real(c_double), intent(in) :: buffer(0:n1)
            real(c_double), intent(in) :: input(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n)
            integer(c_int), value, intent(in) :: intrin  
        end subroutine write_multi_thread_1_wrapper

        subroutine write_multi_thread_2_wrapper(buffer, input, n, ind1, n1, ind2, n2, nthreads, intrin) &
                   BIND(C, name="write_multi_thread_2_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            integer(c_size_t), value, intent(in) :: n2
            integer(c_size_t), value, intent(in) :: nthreads
            real(c_double), intent(in) :: buffer(0:n2)
            real(c_double), intent(in) :: input(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n)
            integer(c_size_t), intent(in) :: ind2(0:n1)
            integer(c_int), value, intent(in) :: intrin
        end subroutine write_multi_thread_2_wrapper


        subroutine writeadd_single_thread_0_wrapper(buffer, input, n, intrin) &
                   BIND(C, name="writeadd_single_thread_0_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            real(c_double), intent(in) :: buffer(0:n)
            real(c_double), intent(in) :: input(0:n)
            integer(c_int), value, intent(in) :: intrin
        end subroutine writeadd_single_thread_0_wrapper

        subroutine writeadd_single_thread_1_wrapper(buffer, input, n, ind1, n1, intrin) &
                   BIND(C, name="writeadd_single_thread_1_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            real(c_double), intent(in) :: buffer(0:n1)
            real(c_double), intent(in) :: input(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n)
            integer(c_int), value, intent(in) :: intrin 
        end subroutine writeadd_single_thread_1_wrapper

        subroutine writeadd_single_thread_2_wrapper(buffer, input, n, ind1, n1, ind2, n2, intrin) &
                   BIND(C, name="writeadd_single_thread_2_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            integer(c_size_t), value, intent(in) :: n2
            real(c_double), intent(in) :: buffer(0:n2)
            real(c_double), intent(in) :: input(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n)
            integer(c_size_t), intent(in) :: ind2(0:n1)
            integer(c_int), value, intent(in) :: intrin
        end subroutine writeadd_single_thread_2_wrapper

        subroutine writeadd_multi_thread_0_wrapper(buffer, input, n, nthreads, intrin) &
                   BIND(C, name="writeadd_multi_thread_0_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: nthreads
            real(c_double), intent(in) :: buffer(0:n)
            real(c_double), intent(in) :: input(0:n)
            integer(c_int), value, intent(in) :: intrin 
        end subroutine writeadd_multi_thread_0_wrapper

        subroutine writeadd_multi_thread_1_wrapper(buffer, input, n, ind1, n1, nthreads, intrin) &
                   BIND(C, name="writeadd_multi_thread_1_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            integer(c_size_t), value, intent(in) :: nthreads
            real(c_double), intent(in) :: buffer(0:n1)
            real(c_double), intent(in) :: input(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n)
            integer(c_int), value, intent(in) :: intrin  
        end subroutine writeadd_multi_thread_1_wrapper

        subroutine writeadd_multi_thread_2_wrapper(buffer, input, n, ind1, n1, ind2, n2, nthreads, intrin) &
                   BIND(C, name="writeadd_multi_thread_2_wrapper")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            integer(c_size_t), value, intent(in) :: n2
            integer(c_size_t), value, intent(in) :: nthreads
            real(c_double), intent(in) :: buffer(0:n2)
            real(c_double), intent(in) :: input(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n)
            integer(c_size_t), intent(in) :: ind2(0:n1)
            integer(c_int), value, intent(in) :: intrin
        end subroutine writeadd_multi_thread_2_wrapper
    end interface
end module scoria
