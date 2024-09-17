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

        subroutine scoria_read_1(res, buffer, n, ind1, n1) &
                   BIND(C, name="scoria_read_1")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            real(c_double), intent(inout) :: res(0:n1)
            real(c_double), intent(in) :: buffer(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n1)
        end subroutine scoria_read_1

        subroutine scoria_write_1(buffer, input, n, ind1, n1) &
                   BIND(C, name="scoria_write_1")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            real(c_double), intent(inout) :: buffer(0:n)
            real(c_double), intent(in) :: input(0:n1)
            integer(c_size_t), intent(in) :: ind1(0:n1)
        end subroutine scoria_write_1

        subroutine scoria_read_1_many_2(resA, resB, bufferA, bufferB, n, ind1, n1) &
                   BIND(C, name="scoria_read_1_many_2")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            real(c_double), intent(inout) :: resA(0:n1)
            real(c_double), intent(inout) :: resB(0:n1)
            real(c_double), intent(in) :: bufferA(0:n)
            real(c_double), intent(in) :: bufferB(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n1)
        end subroutine scoria_read_1_many_2

        subroutine scoria_read_1_many_3(resA, resB, resC, bufferA, bufferB, bufferC, n, ind1, n1) &
                   BIND(C, name="scoria_read_1_many_3")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            real(c_double), intent(inout) :: resA(0:n1)
            real(c_double), intent(inout) :: resB(0:n1)
            real(c_double), intent(inout) :: resC(0:n1)
            real(c_double), intent(in) :: bufferA(0:n)
            real(c_double), intent(in) :: bufferB(0:n)
            real(c_double), intent(in) :: bufferC(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n1)
        end subroutine scoria_read_1_many_3

        subroutine scoria_read_1_many_4(resA, resB, resC, resD, bufferA, bufferB, bufferC, bufferD, n, ind1, n1) &
                   BIND(C, name="scoria_read_1_many_4")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            real(c_double), intent(inout) :: resA(0:n1)
            real(c_double), intent(inout) :: resB(0:n1)
            real(c_double), intent(inout) :: resC(0:n1)
            real(c_double), intent(inout) :: resD(0:n1)
            real(c_double), intent(in) :: bufferA(0:n)
            real(c_double), intent(in) :: bufferB(0:n)
            real(c_double), intent(in) :: bufferC(0:n)
            real(c_double), intent(in) :: bufferD(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n1)
        end subroutine scoria_read_1_many_4

        subroutine scoria_read_1_many_6(resA, resB, resC, resD, resE, resF, bufferA, bufferB, bufferC, bufferD, bufferE, bufferF, n, ind1, n1) &
                   BIND(C, name="scoria_read_1_many_6")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            real(c_double), intent(inout) :: resA(0:n1)
            real(c_double), intent(inout) :: resB(0:n1)
            real(c_double), intent(inout) :: resC(0:n1)
            real(c_double), intent(inout) :: resD(0:n1)
            real(c_double), intent(inout) :: resE(0:n1)
            real(c_double), intent(inout) :: resF(0:n1)
            real(c_double), intent(in) :: bufferA(0:n)
            real(c_double), intent(in) :: bufferB(0:n)
            real(c_double), intent(in) :: bufferC(0:n)
            real(c_double), intent(in) :: bufferD(0:n)
            real(c_double), intent(in) :: bufferE(0:n)
            real(c_double), intent(in) :: bufferF(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n1)
        end subroutine scoria_read_1_many_6

        subroutine scoria_read_1_many_7(resA, resB, resC, resD, resE, resF, resG, bufferA, bufferB, bufferC, bufferD, bufferE, bufferF, bufferG, n, ind1, n1) &
                   BIND(C, name="scoria_read_1_many_7")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            real(c_double), intent(inout) :: resA(0:n1)
            real(c_double), intent(inout) :: resB(0:n1)
            real(c_double), intent(inout) :: resC(0:n1)
            real(c_double), intent(inout) :: resD(0:n1)
            real(c_double), intent(inout) :: resE(0:n1)
            real(c_double), intent(inout) :: resF(0:n1)
            real(c_double), intent(inout) :: resG(0:n1)
            real(c_double), intent(in) :: bufferA(0:n)
            real(c_double), intent(in) :: bufferB(0:n)
            real(c_double), intent(in) :: bufferC(0:n)
            real(c_double), intent(in) :: bufferD(0:n)
            real(c_double), intent(in) :: bufferE(0:n)
            real(c_double), intent(in) :: bufferF(0:n)
            real(c_double), intent(in) :: bufferG(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n1)
        end subroutine scoria_read_1_many_7

        subroutine scoria_read_ranged_1_b(res, buffer, n, ind1, n1, R) &
                   BIND(C, name="scoria_read_ranged_1_b")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            integer(c_size_t), intent(in) :: R(0:n1)
            real(c_double), intent(inout) :: res(0:n1)
            real(c_double), intent(in) :: buffer(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n1)
        end subroutine scoria_read_ranged_1_b

        subroutine scoria_read_ranged_1_c(res, buffer, n, ind1, n1, R) &
                   BIND(C, name="scoria_read_ranged_1_c")
            use iso_c_binding
            implicit none

            integer(c_size_t), value, intent(in) :: n
            integer(c_size_t), value, intent(in) :: n1
            integer(c_size_t), value, intent(in) :: R
            real(c_double), intent(inout) :: res(0:n1)
            real(c_double), intent(in) :: buffer(0:n)
            integer(c_size_t), intent(in) :: ind1(0:n1)
        end subroutine scoria_read_ranged_1_c

    end interface
end module scoria
