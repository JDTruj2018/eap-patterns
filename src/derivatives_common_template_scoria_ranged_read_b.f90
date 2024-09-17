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

#if defined(__INTEL_COMPILER) || defined(_SEQUOIA) || defined(CRAY)
#define NOT_GCC
#endif
  
#if defined(_SEQUOIA) && defined(__GNUC__)
#undef NOT_GCC
#endif
  
#ifndef NOT_GCC
  /* The fortran preprocessor is in traditional mode, so token concatentation is like below:*/
#define stringify(s1) "s1"
#else
#define stringify(s1) #s1
#endif

#ifdef CRAY
#define stringify(s1) "s1"
#endif

#define stringjoin(s1,s2) (stringify(s1) // "_" // stringify(s2))
#define TIMER_NAME(key) stringjoin(derivatives_common,key)
#define TIMERSET(val,key)

  ! Define this as needed to get timer
#define TIMERSET_UNCOND(val,key) call timerset(val,TIMER_NAME(key))

#define get_value(l,nm) merge(value_cloned(l,nm),invalue(l,nm),l .gt. mesh%cells%numcell)

module my_scoria_ranged_read_b_derivatives
    use iso_fortran_env, only : REAL64
    use define_kind, only: ZERO, HALF, ONE, TWO, HI_SIDE, LO_SIDE, mype, iope
    use sim_types, only : sim_info_t
    use mesh_types, only : mesh_t, cells_t, faces_t
    use mesh_state_types, only : mesh_state_frac_core_t, mesh_state_core_t
    implicit none
    public
    integer, parameter :: NO_DERIV = 0
    integer, parameter :: INT_MM = 1
    integer, parameter :: INT_EMM = 2
    integer, parameter :: INT_LV = 3
    integer, parameter :: INT_NL = 4

    integer :: itrmax = 1
    integer :: method = 0
    integer :: nitr = 1
    logical :: limit_slope = .false.
  contains
    subroutine derivatives_common_split_scoria_ranged_read_b(sim, mesh, &
      frac_core, core, &
      gradp, intopt, &
      cell_dim, numitr, nvec, kode, &
      noslope_cell, deriv, do_fincom, faceval, &
      deriv_weight, invalue, value_cloned, do_special)
      use sim_types, only : sim_info_t
      use mesh_types, only : mesh_t
      use mesh_state_types, only : mesh_state_frac_core_t, mesh_state_core_t
      use gradient_types,         only : gradient_prop_t
      use interface_types,        only : interface_option_t
      class(sim_info_t), intent(in) :: sim
      type(mesh_t), intent(in) :: mesh
      type(mesh_state_frac_core_t), intent(in) :: frac_core
      type(mesh_state_core_t), intent(in) :: core
      type(gradient_prop_t), intent(in) :: gradp
      type(interface_option_t), intent(in) :: intopt
      integer,     intent(in) :: cell_dim
      integer,     intent(in) :: numitr
      integer,     intent(in) :: nvec
      integer,     intent(in) :: kode(:,:)
      logical,     intent(in), allocatable :: noslope_cell(:)
      logical,     intent(in) :: do_fincom
      logical,     intent(in), optional :: do_special

      real(REAL64), intent(out),  dimension(:,:,:) :: deriv

      ! for inflow bndy cond'ns
      real(REAL64), intent(in), optional            :: faceval(:,:)
      real(REAL64), intent(in), optional            :: deriv_weight(:)

      real(REAL64), intent(in),    dimension(:,:)   :: invalue
      real(REAL64), intent(in), optional, dimension(:,:)   :: value_cloned
      
      call derivatives_common_internal_split_scoria_ranged_read_b(sim, mesh, cell_dim, numitr, nvec, kode, &
        noslope_cell, deriv, do_fincom, &
        invalue, &
        faceval=faceval, deriv_weight=deriv_weight, &
        value_cloned=value_cloned, do_special=do_special, &
        frac_core=frac_core, core=core, &
        gradp=gradp, intopt=intopt)
    end subroutine derivatives_common_split_scoria_ranged_read_b

    subroutine derivatives_common_internal_split_scoria_ranged_read_b(sim, mesh, &
      cell_dim, numitr, nvec, kode, &
      noslope_cell, deriv, do_fincom, &
      invalue, &
      faceval, deriv_weight, &
      value_cloned, do_special, &
      frac_core, core, &
      gradp, intopt)
      use my_derivatives
      use sim_types, only : sim_info_t
      use mesh_types, only : mesh_t
      use mesh_state_types, only : mesh_state_frac_core_t, mesh_state_core_t
      use gradient_types,         only : gradient_prop_t
      use interface_types,        only : interface_option_t
      use clone_lib_module,       only : clone_get
      use timer_module  ,         only : timerset
      use util,                   only : global_error
      use var_wrapper_class,      only : var_wrapper, vw_set
#ifdef EAP_KOKKOS_GRADIENTS
      use flcl_mod, only: nd_array_t, to_nd_array
      use iso_c_binding
      use gradients_interfaces
#endif
      ! ------------------------------------------------------------------------------
      ! ----- Purpose
      !       This routine calculates the derivative for variable value for a
      !       variable with dimension (numcell_clone,nvec). When using OpenMP,
      !       it calculates the derivative of value, i.e. dvalue/ddirection,
      !       distributed into threads.

      !       In certain cases, this routine makes an approximation: It considers
      !       zero the derivatives that are very close to zero (TINY_DERIVATIVE).
      !       It will only save time if we have a large number of derivatives close to
      !       zero - abs(value(a)-value(b))/abs(value(a)+value(b)) .le. TINY_DERIVATIVE -
      !       it will make these derivatives = 0 without trying to calculate the exact
      !       value.
      !       The subroutine deriv_details_scoria specify how we are calculating the derivatives.
      ! ------------------------------------------------------------------------------

      ! ----- calling arguments
      class(sim_info_t), intent(in)                :: sim
      type(mesh_t), intent(in)                     :: mesh
      integer,     intent(in)                      :: cell_dim
      integer,     intent(in)                      :: numitr
      integer,     intent(in)                      :: nvec
      integer,     intent(in)                      :: kode(:,:)
      logical,     intent(in), allocatable         :: noslope_cell(:)
      logical,     intent(in)                      :: do_fincom
      logical,     intent(in),    optional         :: do_special
      type(mesh_state_frac_core_t), optional, intent(in) :: frac_core
      type(mesh_state_core_t), optional, intent(in) :: core
      type(gradient_prop_t), optional, intent(in)                 :: gradp
      type(interface_option_t), optional, intent(in)           :: intopt

      real(REAL64), intent(out),  dimension(:,:,:) :: deriv

      ! for inflow bndy cond'ns
      real(REAL64), intent(in), optional            :: faceval(:,:)
      real(REAL64), intent(in), optional            :: deriv_weight(:)

      real(REAL64), intent(in),    dimension(:,:)   :: invalue
      real(REAL64), intent(in), optional, dimension(:,:)   :: value_cloned


      ! ----- local variables

      logical      :: do_pressure

      integer      :: dir, itr, ivec, nm
      integer      :: lvofd
      integer      :: priv_suntop, suntop
      real(REAL64) :: lv_weight = ZERO
      type(var_wrapper) :: vw(1)

      real(REAL64) :: cell_val_flcl(mesh%cells%numcell_clone,2,nvec)
      real(REAL64) :: cell_value_hilo(mesh%cells%numcell_clone,2,nvec)
      real(REAL64) :: deriv_mm(mesh%cells%numcell_clone,nvec)
      real(REAL64) :: cell_val_mnmx_hilo(mesh%cells%numcell_clone,4)
      real(REAL64) :: cell_deriv_hilo(mesh%cells%numcell_clone,2)

      real(REAL64), allocatable :: cell_deriv_hilo_all_dir(:,:,:,:)    !JV left and right extimates of derivatives

#ifdef EAP_KOKKOS_GRADIENTS
      ! ----- interfacing variables

      logical(c_bool) :: ic_present_do_special
      logical(c_bool) :: ic_do_special
      logical(c_bool) :: ic_present_do_pressure
      logical(c_bool) :: ic_do_pressure
      logical(c_bool) :: ic_present_gradp
      logical(c_bool) :: ic_present_deriv_weight
      logical(c_bool) :: mype_is_iope
      logical(c_bool) :: ic_present_faceval
      logical(c_bool) :: ic_limit_slope
      type(nd_array_t) :: ic_faceval_wrapper
      type(nd_array_t) :: ic_deriv_weight_wrapper
      type(gradient_prop_t) :: ic_gradp

      ! ----- dummy variables
      real(REAL64) :: fake_faceval(2,2)
      real(REAL64) :: fake_deriv_weight(2)

#endif

      ! ----- code

      associate (cells => mesh%cells, &
           levs => mesh%levels, &
           faces => mesh%faces)

        TIMERSET_UNCOND(.true.,main)

        !JV allocate an array for left and right extimates of derivatives
        if( numitr.eq.7 .or. numitr.eq.8 ) allocate( cell_deriv_hilo_all_dir(sim%numdim,mesh%cells%numcell_clone,2,nvec) )

        call deriv_details_scoria_ranged_read_b(numitr, lv_weight, cell_dim, do_pressure,  &
             do_special=do_special)

#ifdef EAP_KOKKOS_GRADIENTS
        ! set up interfacing variables after deriv_details_scoria
        if( present(faceval) ) then
           ic_faceval_wrapper = to_nd_array(faceval)
           ic_present_faceval = .true.
        else
           ic_faceval_wrapper = to_nd_array(fake_faceval)
           ic_present_faceval = .false.
        end if

        ic_present_do_pressure = .true.
        if ( do_pressure ) then
           ic_do_pressure = .true.
        else
           ic_do_pressure = .false.
        end if

        if(present(do_special)) then
           ic_present_do_special = .true.
           ic_do_special = logical(do_special, c_bool)
        else
           ic_present_do_special = .false.
           ic_do_special = .false.
        end if

        if(present(gradp)) then
           ic_present_gradp = .true.
           ic_gradp = gradp
        else
           ic_present_gradp = .false.
           ! ic3_gradp will be uninitialized, don't use it
        end if

        if(present(deriv_weight)) then
           ic_present_deriv_weight = .true.
           ic_deriv_weight_wrapper = to_nd_array(deriv_weight)
        else
           ic_present_deriv_weight = .false.
           ic_deriv_weight_wrapper = to_nd_array(fake_deriv_weight)
        end if

        mype_is_iope = .false.
        if (mype.eq.iope) mype_is_iope = .true.

        ic_limit_slope = logical(limit_slope, c_bool)
#endif
        call inside_com1_split(sim, mesh, core, nvec, &
             & kode, deriv, cell_val_flcl, invalue=invalue, value_cloned=value_cloned)

        if (itrmax.gt.0) then
           !         calculate all, from now on:

#ifdef _DEBUG
            do ivec   = 1,nvec
              do dir = 1,sim%numdim
                 if(kode(dir,ivec).eq.-2)then
                    if(.not.present(faceval))then
                       call global_error(                                    &
                            'DERIVATIVES: optional argument missing for kode = -2')
                    endif

                 endif
              enddo
            enddo
#endif

           lvofd = 1
           if (intopt%interface_option .ge. 3) lvofd = mesh%cells%numcell_clone

           do dir    = 1,sim%numdim
              do itr  = 1,itrmax
                 TIMERSET_UNCOND(.true.,inner)
                 ! ......      Estimate values at edges of cells
#ifdef EAP_KOKKOS_GRADIENTS
                 call inside_com2_split_arrays(  itr, dir, nvec, &
                      & mesh%levels%allnumtop, &
                      & mesh%cells%numcell_clone, &
                      & mesh%cells%numcell, mype_is_iope, &
                      & to_nd_array(mesh%levels%alltop), &
                      & to_nd_array(deriv), &
                      & to_nd_array(deriv_mm), &
                      & to_nd_array(cell_value_hilo), &
                      & to_nd_array(invalue), &
                      & to_nd_array(value_cloned), &
                      & to_nd_array(cells%cell_half_lo), &
                      & to_nd_array(cells%cell_half_hi) )
#else
                 ! else ifdef EAP_KOKKOS_GRADIENTS
                 ! no EAP_KOKKOS_GRADIENTS - so we're running F only code
                 call inside_com2_split(mesh, itr, dir, nvec, deriv, deriv_mm, &
                      & cell_value_hilo, invalue=invalue, value_cloned=value_cloned)
#endif
                 !endif EAP_KOKKOS_GRADIENTS

                 !             Not necessary to communicate cell_value_hilo here because they have
                 !             been estimated in the ghost cells.
#ifdef EAP_KOKKOS_GRADIENTS
                 ! don't zero out cell_val_mnmx_hilo (initialized in inside_com3a)
                 ! don't zero out cell_deriv_hilo (zeroed out in inside_com3c)
                 do nm = 1,nvec
                    ! .....    calculate the area weighted average face values
                    if (mesh%levels%allnumtop.gt.0) then
                       call inside_com3a_split_arrays( mype_is_iope, dir, nm, cells%numcell,&
                            & cells%numcell_clone, huge(ONE), to_nd_array(faces%face_num), &
                            & to_nd_array(faces%face_id), to_nd_array(faces%face_lo), &
                            & to_nd_array(faces%face_hi), to_nd_array(faces%face_local), &
                            & to_nd_array(cell_val_mnmx_hilo), to_nd_array(invalue), &
                            & to_nd_array(value_cloned) )
                       call inside_com3b_arrays( mype_is_iope, ic_present_faceval, &
                            & ic_present_do_special, ic_do_special, ic_present_do_pressure, &
                            & ic_do_pressure, dir, nm, to_nd_array(faces%face_num), &
                            & to_nd_array(faces%face_id), to_nd_array(faces%face_lo), &
                            & to_nd_array(faces%face_hi), &
                            & to_nd_array(faces%face_local), to_nd_array(kode), &
                            & to_nd_array(cell_val_mnmx_hilo), to_nd_array(cell_value_hilo), &
                            & ic_faceval_wrapper, to_nd_array(core%rho), &
                            & to_nd_array(cells%cell_half_lo), &
                            & to_nd_array(cells%cell_half_hi) )
                       call inside_com3c_split_arrays( mype_is_iope, dir, nm, &
                            & levs%allnumtop, cells%numcell, to_nd_array(levs%alltop), &
                            & to_nd_array(cells%cell_half_lo), to_nd_array(cells%cell_half_hi), &
                            & to_nd_array(cell_val_mnmx_hilo), to_nd_array(cell_deriv_hilo), &
                            & to_nd_array(invalue), to_nd_array(value_cloned) )
                       call inside_com3d_split_arrays( mype_is_iope, dir, nm, cells%numcell, &
                            & to_nd_array(faces%face_num), to_nd_array(faces%face_id), &
                            & to_nd_array(faces%face_lo), to_nd_array(faces%face_hi), &
                            & to_nd_array(faces%face_local), to_nd_array(faces%face_flag), &
                            & to_nd_array(kode), to_nd_array(cells%cell_half_lo), &
                            & to_nd_array(cells%cell_half_hi), to_nd_array(cell_deriv_hilo), &
                            & to_nd_array(invalue), to_nd_array(value_cloned) )
                       call inside_com3e_arrays( mype_is_iope, ic_limit_slope, itr, nm, &
                            & levs%allnumtop, to_nd_array(levs%alltop), &
                            & to_nd_array(cell_deriv_hilo), to_nd_array(deriv_mm) )
                       call inside_com3f_arrays( mype_is_iope, dir, nm, numitr, &
                            & method, levs%allnumtop, lv_weight, &
                            & to_nd_array(levs%alltop), to_nd_array(deriv), &
                            & to_nd_array(cell_deriv_hilo), to_nd_array(cell_deriv_hilo_all_dir) )
                    end if
                 end do
                 call inside_com3g_arrays( mype_is_iope, ic_limit_slope, &
                      & ic_present_gradp, ic_present_deriv_weight, dir, nvec, &
                      & itr, itrmax, levs%allnumtop, ic_gradp, &
                      & to_nd_array(levs%alltop), to_nd_array(deriv), &
                      & to_nd_array(deriv_mm), ic_deriv_weight_wrapper )
                 call inside_com3h(sim, mesh, dir, numitr, nvec, &
                      & kode, deriv, frac_core=frac_core, intopt=intopt )
#else
                 ! else ifdef EAP_KOKKOS_GRADIENTS
                 ! no EAP_KOKKOS_GRADIENTS - so we're running only F code
                 ! don't zero out cell_val_mnmx_hilo (initialized in inside_com3a)
                 ! don't zero out cell_deriv_hilo (zeroed out in inside_com3c)
                 do nm = 1,nvec
                    ! .....    calculate the area weighted average face values
                    if (mesh%levels%allnumtop.gt.0) then
                       call inside_com3A_split(sim, mesh, dir, nm, cell_val_mnmx_hilo, &
                            & invalue, value_cloned=value_cloned )
                       call inside_com3b_scoria_ranged_read_b(sim, mesh, dir, nm, cell_val_mnmx_hilo, &
                            & kode, cell_value_hilo, do_special=do_special, &
                            & do_pressure=do_pressure, core=core, faceval=faceval )
                       call inside_com3C_split(sim, mesh, dir, nm, cell_deriv_hilo, &
                            & cell_val_mnmx_hilo, invalue, value_cloned=value_cloned )
                       call inside_com3D_split(sim, mesh, dir, nm, kode, cell_deriv_hilo, &
                            & invalue, value_cloned=value_cloned )
                       call inside_com3e(sim, mesh, itr, nm, limit_slope, &
                            & cell_deriv_hilo, deriv_mm )
                       call inside_com3f(sim, mesh, dir, nm, numitr, method, &
                            & lv_weight, deriv, cell_deriv_hilo, &
                            & cell_deriv_hilo_all_dir )
                    end if
                 end do
                 call inside_com3g(mesh, dir, nvec, itr, itrmax, limit_slope, &
                      & deriv, deriv_mm, gradp=gradp, deriv_weight=deriv_weight )
                 call inside_com3h(sim, mesh, dir, numitr, nvec, &
                      & kode, deriv, frac_core=frac_core, intopt=intopt )
#endif
                 ! endif EAP_KOKKOS_GRADIENTS
                 TIMERSET_UNCOND(.false.,inner)
                 if (itr.lt.itrmax) then
                    TIMERSET(.true., clone_get_1)
                    call clone_get(deriv(:cells%numcell_clone,dir,:nvec),cells%numcell_clone,nvec)
                    TIMERSET(.false., clone_get_1)
                 endif

              enddo ! itr
           enddo ! dir

           call inside_com4_split(sim, mesh, nvec, numitr, lv_weight, deriv, &
                & cell_val_flcl, cell_deriv_hilo_all_dir, invalue, &
                & value_cloned=value_cloned )

           ! ....    deallocate variables allocated before com1:

        endif ! itrmax

        !       Do the communications necessary before calling hydro_face_nocomm
        !       Communicate both value and deriv before calling hydro_face_nocomm

        if (do_fincom) then
           !           We do not need to communicate "value" here. It has been already
           !           communicated in the beginning of the routine.
           TIMERSET(.true., clone_get_2)
           call vw_set(vw(1), deriv, cells%numcell, sim%numdim, nvec)
           call clone_get(vw, 1)
           TIMERSET(.false., clone_get_2)
        endif

        if( numitr.eq.7 .or. numitr.eq.8 ) deallocate( cell_deriv_hilo_all_dir ) !JV

        !SS TIMERSET(.true., post_modify_slope)
        !SS if (allocated(noslope_cell)) &
        !SS      & call post_modify_slope(cell_dim,sim%numdim,nvec,noslope_cell,deriv)
        !SS TIMERSET(.false., post_modify_slope)

        TIMERSET_UNCOND(.false.,main)

      end associate

    contains
      ! ------------------------------------------------------------------------------
      subroutine deriv_details_scoria_ranged_read_b(numitr, lv_weight, cell_dim, &
           & do_pressure, do_special)

        use util       , only : global_error

        integer,      intent(in)     :: cell_dim
        integer,      intent(in)     :: numitr
        logical,      intent(out)    :: do_pressure
        real(REAL64), intent(inout)  :: lv_weight

        logical,      intent(in), optional :: do_special

        logical     :: dimerror = .false.

        !       numitr = 0 : no derivatives (1st order accurate)
        !              = 1 : MM
        !              = 2 : iterated MM
        !              = 3 : extended MM
        !              = 4 : LV
        !              = 5 : NL
        !              = 6 : modified LV

        do_pressure = .false.
        if (present(do_special)) then
           do_pressure = do_special
        endif

        nitr = numitr
        if (nitr.lt.0) nitr = -nitr

        ! ..... select method
        method = NO_DERIV
        itrmax = 0
        limit_slope = .false.

        select case (nitr)
        case (0)
           method = NO_DERIV
           itrmax = 0
           limit_slope = .false.
        case (1)
           method = INT_MM
           itrmax = 1
           limit_slope = .false.
        case (2)
           method = INT_MM
           itrmax = 2
           limit_slope = .true.
        case (3)
           method = INT_EMM
           itrmax = 1
           limit_slope = .true.
        case (4)
           method = INT_LV
           itrmax = 1
           limit_slope = .true.
           lv_weight = TWO
        case (5)
           method = INT_NL
           itrmax = 1
           limit_slope = .false.
        case (6,8)                
           method = INT_LV
           itrmax = 1
           limit_slope = .true.
           lv_weight = 1.5_REAL64
        case (7)                  
           method = INT_LV
           itrmax = 1
           limit_slope = .true.
           lv_weight = sqrt(3.0_REAL64)  !JV this value is reduced in inside_com4
        case default
           dimerror = .true.
        end select

        if (dimerror) then
           if (mype.eq.iope) write(*,*)'$$$ DERIVATIVES: nitr = ',nitr
           call global_error('DERIVATIVES: bad value for nitr')
        endif

        if (cell_dim .ne. mesh%cells%numcell_clone) then
           call global_error( &
                & 'derivatives_common_split: cell_dim must be equal to numcell_clone')
        endif

      end subroutine deriv_details_scoria_ranged_read_b
      ! ------------------------------------------------------------------------------
    end subroutine derivatives_common_internal_split_scoria_ranged_read_b

    subroutine inside_com3b_scoria_ranged_read_b(sim, mesh, dir, nm, cell_val_mnmx_hilo, kode, &
      & cell_value_hilo, do_special, do_pressure, core, faceval)
      use iso_c_binding
      use iso_fortran_env, only: REAL64, INT64, INT32, INT8
      use scoria, only: scoria_read_1, scoria_write_1, scoria_read_ranged_1_b

      ! non-optional scalars and derived types
      class(sim_info_t), intent(in) :: sim 
      type(mesh_t), intent(in) :: mesh
      integer, intent(in) :: dir, nm

      ! non-optional arrays
      real(REAL64), intent(inout), dimension(:,:) :: cell_val_mnmx_hilo
      integer, intent(in) :: kode(:,:)
      real(REAL64), intent(inout) :: cell_value_hilo(:,:,:)

      ! optional scalars and derived types
      logical, intent(in), optional :: do_special 
      logical, intent(in), optional :: do_pressure
      type(mesh_state_core_t), optional, intent(in) :: core 

      ! optional arrays
      ! for inflow bndy cond'ns
      real(REAL64), intent(in), optional :: faceval(:,:)

      ! ----- local variables
      integer :: loop, n, lhi, llo, lmi, lmo 
      real(REAL64) :: face_value
      real(REAL64) :: my_dt, start, end
      integer(INT64) :: c_now
      real(REAL64) :: c_rate

      ! ----- scoria variables
      integer(INT64) :: i, count, prev, range
      integer(c_size_t) :: scoria_n, transformed_hi_n, transformed_lo_n

      integer(c_size_t), allocatable, save :: face_local_hi_scoria(:), face_local_lo_scoria(:), &
        & face_local_hi_scoria_transformed(:), face_local_lo_scoria_transformed(:), &
        & range_face_local_hi_scoria(:), range_face_local_lo_scoria(:)

      real(REAL64), allocatable, save :: cell_half_lo_scoria(:), rho_scoria(:), &
        & cell_value_hilo_1_scoria(:), cell_value_hilo_2_scoria(:), cell_half_hi_scoria(:), &
        & cell_val_mnmx_hilo_hi_1_scoria(:), cell_val_mnmx_hilo_hi_2_scoria(:), &
        & cell_val_mnmx_hilo_lo_3_scoria(:), cell_val_mnmx_hilo_lo_4_scoria(:), &
        & cell_val_mnmx_hilo_mi_1_scoria(:), cell_val_mnmx_hilo_mi_2_scoria(:), &
        & cell_val_mnmx_hilo_mo_3_scoria(:), cell_val_mnmx_hilo_mo_4_scoria(:)

      real(REAL64), allocatable, save :: cell_half_lo_hi_packed(:), rho_hi_packed(:), &
        & cell_value_hilo_lo_packed(:), cell_half_hi_lo_packed(:), rho_lo_packed(:), &
        & cell_value_hilo_hi_packed(:), cell_val_mnmx_hilo_hi_1_packed(:), cell_val_mnmx_hilo_hi_2_packed(:), &
        & cell_val_mnmx_hilo_lo_3_packed(:), cell_val_mnmx_hilo_lo_4_packed(:), cell_val_mnmx_hilo_mi_1_packed(:), &
        & cell_val_mnmx_hilo_mi_2_packed(:), cell_val_mnmx_hilo_mo_3_packed(:), cell_val_mnmx_hilo_mo_4_packed(:)

      associate (cells => mesh%cells, &
        faces => mesh%faces)

      !open (unit=10,file="scoria-output.txt",action="write")

      TIMERSET(.true., inside_com3b)

      ! ..... calculate the area weighted average face values
      do loop = 1,faces%face_num(dir)
        ! Scoria Setup
        scoria_n = faces%face_hi(loop,dir) - faces%face_lo(loop,dir) + 1

        ! Scoria: Map to 1-D Arrays
        allocate(face_local_hi_scoria(scoria_n))
        allocate(face_local_lo_scoria(scoria_n))

        count = 1
        do i = faces%face_lo(loop,dir), faces%face_hi(loop,dir)
            face_local_hi_scoria(count) = faces%face_local(i, HI_SIDE, dir)
            face_local_lo_scoria(count) = faces%face_local(i, LO_SIDE, dir)
            count = count + 1
        enddo

        allocate(cell_half_lo_scoria(size(cells%cell_half_lo(:, dir))))
        allocate(cell_value_hilo_1_scoria(size(cell_value_hilo(:, 1, nm))))
        allocate(cell_value_hilo_2_scoria(size(cell_value_hilo(:, 2, nm))))
        allocate(cell_half_hi_scoria(size(cells%cell_half_hi(:, dir))))

        if (present(do_special)) then
            allocate(rho_scoria(size(core%rho)))
        endif

        count = 1
        do i = lbound(cells%cell_half_lo, 1), ubound(cells%cell_half_lo, 1)
            cell_half_lo_scoria(count) = cells%cell_half_lo(i, dir)
            count = count + 1
        enddo

        if (present(do_special)) then
          count = 1
          do i = lbound(core%rho, 1), ubound(core%rho, 1)
              rho_scoria(count) = core%rho(i)
              count = count + 1
          enddo
        endif

        count = 1
        do i = lbound(cell_value_hilo, 1), ubound(cell_value_hilo, 1)
            cell_value_hilo_1_scoria(count) = cell_value_hilo(i, 1, nm)
            count = count + 1
        enddo

        count = 1
        do i = lbound(cell_value_hilo, 1), ubound(cell_value_hilo, 1)
            cell_value_hilo_2_scoria(count) = cell_value_hilo(i, 2, nm)
            count = count + 1
        enddo

        count = 1
        do i = lbound(cells%cell_half_hi, 1), ubound(cells%cell_half_hi, 1)
            cell_half_hi_scoria(count) = cells%cell_half_hi(i, dir)
          count = count + 1
        enddo 

        allocate(cell_val_mnmx_hilo_hi_1_scoria(size(cell_val_mnmx_hilo(:, 1))))
        allocate(cell_val_mnmx_hilo_hi_2_scoria(size(cell_val_mnmx_hilo(:, 2))))
        allocate(cell_val_mnmx_hilo_lo_3_scoria(size(cell_val_mnmx_hilo(:, 3))))
        allocate(cell_val_mnmx_hilo_lo_4_scoria(size(cell_val_mnmx_hilo(:, 4))))
        allocate(cell_val_mnmx_hilo_mi_1_scoria(size(cell_val_mnmx_hilo(:, 1))))
        allocate(cell_val_mnmx_hilo_mi_2_scoria(size(cell_val_mnmx_hilo(:, 2))))
        allocate(cell_val_mnmx_hilo_mo_3_scoria(size(cell_val_mnmx_hilo(:, 3))))
        allocate(cell_val_mnmx_hilo_mo_4_scoria(size(cell_val_mnmx_hilo(:, 4))))

        ! Scoria: Allocate Packed Arrays
        allocate(cell_half_lo_hi_packed(scoria_n))
        allocate(cell_value_hilo_lo_packed(scoria_n))
        allocate(cell_half_hi_lo_packed(scoria_n))
        allocate(cell_value_hilo_hi_packed(scoria_n))

        if (present(do_special)) then
          allocate(rho_hi_packed(scoria_n))
          allocate(rho_lo_packed(scoria_n))
        endif

        allocate(cell_val_mnmx_hilo_hi_1_packed(scoria_n))
        allocate(cell_val_mnmx_hilo_hi_2_packed(scoria_n))
        allocate(cell_val_mnmx_hilo_lo_3_packed(scoria_n))
        allocate(cell_val_mnmx_hilo_lo_4_packed(scoria_n))
        allocate(cell_val_mnmx_hilo_mi_1_packed(scoria_n))
        allocate(cell_val_mnmx_hilo_mi_2_packed(scoria_n))
        allocate(cell_val_mnmx_hilo_mo_3_packed(scoria_n))
        allocate(cell_val_mnmx_hilo_mo_4_packed(scoria_n))

        transformed_hi_n = 0
        range = 0
        prev = 0
        do i = lbound(face_local_hi_scoria, 1), ubound(face_local_hi_scoria, 1)
          if (range .eq. 0) then
            transformed_hi_n = transformed_hi_n + 1
          else 
            if ((face_local_hi_scoria(i) - prev) .eq. 1) then
              range = range + 1
            else
              range = 0
              transformed_hi_n = transformed_hi_n + 1
            endif
          endif
          prev = face_local_hi_scoria(i)
        enddo

        transformed_lo_n = 0
        range = 0
        prev = 0
        do i = lbound(face_local_lo_scoria, 1), ubound(face_local_lo_scoria, 1)
          if (range .eq. 0) then
            transformed_lo_n = transformed_lo_n + 1
          else 
            if ((face_local_lo_scoria(i) - prev) .eq. 1) then
              range = range + 1
            else
              range = 0
              transformed_lo_n = transformed_lo_n + 1
            endif
          endif
          prev = face_local_lo_scoria(i)
        enddo

        allocate(face_local_hi_scoria_transformed(transformed_hi_n))
        allocate(face_local_lo_scoria_transformed(transformed_lo_n))

        allocate(range_face_local_hi_scoria(transformed_hi_n))
        allocate(range_face_local_lo_scoria(transformed_lo_n))

        count = 0
        range = 0
        prev = 0
        do i = lbound(face_local_hi_scoria, 1), ubound(face_local_hi_scoria, 1)
          if (range .eq. 0) then
            face_local_hi_scoria_transformed(count) = face_local_hi_scoria(i)
            count = count + 1
          else 
            if ((face_local_hi_scoria(i) - prev) .eq. 1) then
              range = range + 1
            else
              range_face_local_hi_scoria(count) = range 
              range = 0
              count = count + 1
            endif
          endif
          prev = face_local_hi_scoria(i)
        enddo

        count = 0
        range = 0
        prev = 0
        do i = lbound(face_local_lo_scoria, 1), ubound(face_local_lo_scoria, 1)
          if (range .eq. 0) then
            face_local_lo_scoria_transformed(count) = face_local_lo_scoria(i)
            count = count + 1
          else 
            if ((face_local_lo_scoria(i) - prev) .eq. 1) then
              range = range + 1
            else
              range_face_local_lo_scoria(count) = range 
              range = 0
              count = count + 1
            endif
          endif
          prev = face_local_lo_scoria(i)
        enddo

        ! End Scoria Setup

        ! BEGIN SCORIA TIMING
        call system_clock(c_now, c_rate)
        start = real(c_now, REAL64) / c_rate
        TIMERSET(.true., scoria_inside_com3b)

        ! Scoria Read (1-level Gather)
        call scoria_read_ranged_1_b(cell_half_lo_hi_packed, cell_half_lo_scoria, &
          & int(size(cell_half_lo_scoria), kind=c_size_t), face_local_hi_scoria_transformed, &
          & transformed_hi_n, range_face_local_hi_scoria)
        call scoria_read_ranged_1_b(cell_value_hilo_lo_packed, cell_value_hilo_2_scoria, & 
          & int(size(cell_value_hilo_2_scoria), kind=c_size_t), face_local_lo_scoria_transformed, &
          & transformed_lo_n, range_face_local_lo_scoria)

        call scoria_read_ranged_1_b(cell_half_hi_lo_packed, cell_half_hi_scoria, & 
          & int(size(cell_half_hi_scoria), kind=c_size_t), face_local_lo_scoria_transformed, &
          & transformed_lo_n, range_face_local_lo_scoria)
        call scoria_read_ranged_1_b(cell_value_hilo_hi_packed, cell_value_hilo_1_scoria, &
          & int(size(cell_value_hilo_1_scoria), kind=c_size_t), face_local_hi_scoria_transformed, &
          & transformed_hi_n, range_face_local_hi_scoria)

        if (present(do_special)) then
          call scoria_read_ranged_1_b(rho_hi_packed, rho_scoria, int(size(rho_scoria), kind=c_size_t), &
            & face_local_hi_scoria_transformed, transformed_hi_n, range_face_local_hi_scoria)
          call scoria_read_ranged_1_b(rho_lo_packed, rho_scoria, int(size(rho_scoria), kind=c_size_t), & 
          & face_local_lo_scoria_transformed, transformed_lo_n, range_face_local_lo_scoria)
        endif

        call scoria_read_ranged_1_b(cell_val_mnmx_hilo_hi_1_packed, cell_val_mnmx_hilo_hi_1_scoria, &
          & int(size(cell_val_mnmx_hilo_hi_1_scoria), kind=c_size_t), face_local_hi_scoria_transformed, &
          & transformed_hi_n, range_face_local_hi_scoria)
        call scoria_read_ranged_1_b(cell_val_mnmx_hilo_hi_2_packed, cell_val_mnmx_hilo_hi_2_scoria, &
          & int(size(cell_val_mnmx_hilo_hi_2_scoria), kind=c_size_t), face_local_hi_scoria_transformed, &
          & transformed_hi_n, range_face_local_hi_scoria)

        call scoria_read_ranged_1_b(cell_val_mnmx_hilo_lo_3_packed, cell_val_mnmx_hilo_lo_3_scoria, &
          & int(size(cell_val_mnmx_hilo_lo_3_scoria), kind=c_size_t), face_local_lo_scoria_transformed, & 
          & transformed_lo_n, range_face_local_lo_scoria)
        call scoria_read_ranged_1_b(cell_val_mnmx_hilo_lo_4_packed, cell_val_mnmx_hilo_lo_4_scoria, &
          & int(size(cell_val_mnmx_hilo_lo_4_scoria), kind=c_size_t), face_local_lo_scoria_transformed, &
          & transformed_lo_n, range_face_local_lo_scoria)

        call scoria_read_ranged_1_b(cell_val_mnmx_hilo_mi_1_packed, cell_val_mnmx_hilo_mi_1_scoria, &
          & int(size(cell_val_mnmx_hilo_mi_1_scoria), kind=c_size_t), face_local_hi_scoria_transformed, &
          & transformed_hi_n, range_face_local_hi_scoria)
        call scoria_read_ranged_1_b(cell_val_mnmx_hilo_mi_2_packed, cell_val_mnmx_hilo_mi_2_scoria, &
          & int(size(cell_val_mnmx_hilo_mi_2_scoria), kind=c_size_t), face_local_hi_scoria_transformed, &
          & transformed_hi_n, range_face_local_hi_scoria)

        call scoria_read_ranged_1_b(cell_val_mnmx_hilo_mo_3_packed, cell_val_mnmx_hilo_mo_3_scoria, &
          & int(size(cell_val_mnmx_hilo_mo_3_scoria), kind=c_size_t), face_local_lo_scoria_transformed, &
          & transformed_lo_n, range_face_local_lo_scoria)
        call scoria_read_ranged_1_b(cell_val_mnmx_hilo_mo_4_packed, cell_val_mnmx_hilo_mo_4_scoria, &
          & int(size(cell_val_mnmx_hilo_mo_4_scoria), kind=c_size_t), face_local_lo_scoria_transformed, &
          & transformed_lo_n, range_face_local_lo_scoria)

        if (faces%face_id(loop,dir) .gt. 2) then 
          if (present(do_special)) then
            ! Do the HI_SIDE part
            do i = 1, scoria_n
              face_value = ZERO

              if (rho_lo_packed(i) .gt. ZERO .or. rho_hi_packed(i) .gt. ZERO) then 
                if (do_pressure .and. cell_value_hilo_lo_packed(i) &
                  & * cell_value_hilo_hi_packed(i) .le. ZERO) then
                  ! ..... this coding addresses the hot spot problem
                  face_value = (cell_half_lo_hi_packed(i) * &
                    & rho_hi_packed(i) * &
                    & cell_value_hilo_lo_packed(i) + &
                    & cell_half_hi_lo_packed(i) * &
                    & rho_lo_packed(i) * &
                    & cell_value_hilo_hi_packed(i)) / &
                    & (cell_half_lo_hi_packed(i) * &
                    & rho_hi_packed(i) + &
                    & cell_half_hi_lo_packed(i) * &
                    & rho_lo_packed(i))
                else
                  face_value = (cell_half_lo_hi_packed(i) * &
                    & cell_value_hilo_lo_packed(i) + &
                    & cell_half_hi_lo_packed(i) * &
                    & cell_value_hilo_hi_packed(i)) / &
                    & (cell_half_lo_hi_packed(i) + &
                    & cell_half_hi_lo_packed(i))
                endif ! do_pressure ...
              endif ! core%rho

              !write(10, *) "1 Face Value: "
              !write(10, "(f8.2)") face_value

              cell_val_mnmx_hilo_hi_1_packed(i) = min(cell_val_mnmx_hilo_hi_1_packed(i), face_value)
              cell_val_mnmx_hilo_hi_2_packed(i) = min(cell_val_mnmx_hilo_hi_2_packed(i), face_value)
            enddo ! n

            ! Do the LO_SIDE part
            do i = 1, scoria_n
              face_value = ZERO 

              if (rho_lo_packed(i) .gt. ZERO .or. rho_hi_packed(i) .gt. ZERO) then
                if (do_pressure .and. cell_value_hilo_lo_packed(i) * cell_value_hilo_hi_packed(i) .le. ZERO) then
                  ! ..... this coding addresses the hot spot problem
                  face_value = (cell_half_lo_hi_packed(i) * &
                    & rho_hi_packed(i) * &
                    & cell_value_hilo_lo_packed(i) + &
                    & cell_half_hi_lo_packed(i) * &
                    & rho_lo_packed(i) * &
                    & cell_value_hilo_hi_packed(i)) / &
                    & (cell_half_lo_hi_packed(i) * &
                    & rho_hi_packed(i) + &
                    & cell_half_hi_lo_packed(i) * &
                    & rho_lo_packed(i))
                else
                  face_value = (cell_half_lo_hi_packed(i) * &
                    & cell_value_hilo_lo_packed(i) + &
                    & cell_half_hi_lo_packed(i) * &
                    & cell_value_hilo_hi_packed(i)) /&
                    & (cell_half_lo_hi_packed(i) + &
                    & cell_half_hi_lo_packed(i))
                endif ! do_pressure ...
              endif ! core%rho

              !write(10, *) "2 Face Value: "
              !write(10, "(f8.2)") face_value

              cell_val_mnmx_hilo_lo_3_packed(i) = min(cell_val_mnmx_hilo_lo_3_packed(i), face_value)
              cell_val_mnmx_hilo_lo_4_packed(i) = min(cell_val_mnmx_hilo_lo_4_packed(i), face_value)
            enddo ! n

          else ! not present(do_special)
            ! Do the HI_SIDE first
            do i = 1, scoria_n
              face_value = (cell_half_lo_hi_packed(i) * &
                & cell_value_hilo_lo_packed(i) + &
                & cell_half_hi_lo_packed(i) * &
                & cell_value_hilo_hi_packed(i)) / &
                & (cell_half_lo_hi_packed(i) + &
                & cell_half_hi_lo_packed(i))

              !write(10, *) "3 Face Value: "
              !write(10, "(f8.2)") face_value

              cell_val_mnmx_hilo_hi_1_packed(i) = min(cell_val_mnmx_hilo_hi_1_packed(i), face_value)
              cell_val_mnmx_hilo_hi_2_packed(i) = min(cell_val_mnmx_hilo_hi_2_packed(i), face_value)
            enddo ! n
            
            ! Do the LO_SIDE
            do i = 1, scoria_n
              face_value = (cell_half_lo_hi_packed(i) * &
                & cell_value_hilo_lo_packed(i) + &
                & cell_half_hi_lo_packed(i) * &
                & cell_value_hilo_hi_packed(i)) / &
                & (cell_half_lo_hi_packed(i) + &
                & cell_half_hi_lo_packed(i))

              !write(10, *) "4 Face Value: "
              !write(10, "(f8.2)") face_value

              cell_val_mnmx_hilo_lo_3_packed(i) = min(cell_val_mnmx_hilo_lo_3_packed(i), face_value)
              cell_val_mnmx_hilo_lo_4_packed(i) = min(cell_val_mnmx_hilo_lo_4_packed(i), face_value)
            enddo !n
          endif ! present(do_special)
        
        else if (faces%face_id(loop,dir) .eq. 2) then
          if (kode(dir,nm) .eq. -2) then 
            do i = 1, scoria_n
              face_value = faceval(i + faces%face_lo(loop,dir),dir)

              !write(10, *) "5 Face Value: "
              !write(10, "(f8.2)") face_value

              cell_val_mnmx_hilo_mo_3_packed(i) = min(cell_val_mnmx_hilo_mo_3_packed(i), face_value)
              cell_val_mnmx_hilo_mo_4_packed(i) = min(cell_val_mnmx_hilo_mo_4_packed(i), face_value)
            enddo ! n
          else
            do i = 1, scoria_n
              face_value = cell_value_hilo_lo_packed(i)

              !write(10, *) "6 Face Value: "
              !write(10, "(f8.2)") face_value

              cell_val_mnmx_hilo_lo_3_packed(i) = min(cell_val_mnmx_hilo_lo_3_packed(i), face_value)
              cell_val_mnmx_hilo_lo_4_packed(i) = min(cell_val_mnmx_hilo_lo_4_packed(i), face_value)
            enddo ! n
          endif

        else if (faces%face_id(loop,dir) .eq. 1) then
          if (kode(dir,nm) .eq. -2) then
            do i = 1, scoria_n
              face_value = faceval(i + faces%face_lo(loop,dir), dir)

              !write(10, *) "7 Face Value: "
              !write(10, "(f8.2)") face_value

              cell_val_mnmx_hilo_mi_1_packed(i) = min(cell_val_mnmx_hilo_mi_1_packed(i), face_value)
              cell_val_mnmx_hilo_mi_2_packed(i) = min(cell_val_mnmx_hilo_mi_2_packed(i), face_value)
            enddo ! n
          else
            do i = 1, scoria_n
              face_value = cell_value_hilo_hi_packed(i)

              !write(10, *) "8 Face Value: "
              !write(10, "(f8.2)") face_value

              cell_val_mnmx_hilo_hi_1_packed(i) = min(cell_val_mnmx_hilo_hi_1_packed(i), face_value)
              cell_val_mnmx_hilo_hi_2_packed(i) = min(cell_val_mnmx_hilo_hi_2_packed(i), face_value)
            enddo ! n
          endif
        endif ! face_id

        call scoria_write_1(cell_val_mnmx_hilo_hi_1_scoria, cell_val_mnmx_hilo_hi_1_packed, &
        & int(size(cell_val_mnmx_hilo_hi_1_scoria), 8), face_local_hi_scoria, scoria_n)
        call scoria_write_1(cell_val_mnmx_hilo_hi_2_scoria, cell_val_mnmx_hilo_hi_2_packed, &
        & int(size(cell_val_mnmx_hilo_hi_2_scoria), 8), face_local_hi_scoria, scoria_n)

        call scoria_write_1(cell_val_mnmx_hilo_lo_3_scoria, cell_val_mnmx_hilo_lo_3_packed, &
        & int(size(cell_val_mnmx_hilo_lo_3_scoria), 8), face_local_lo_scoria, scoria_n)
        call scoria_write_1(cell_val_mnmx_hilo_lo_4_scoria, cell_val_mnmx_hilo_lo_4_packed, &
        & int(size(cell_val_mnmx_hilo_lo_4_scoria), 8), face_local_lo_scoria, scoria_n)

        call scoria_write_1(cell_val_mnmx_hilo_mi_1_scoria, cell_val_mnmx_hilo_mi_1_packed, &
        & int(size(cell_val_mnmx_hilo_mi_1_scoria), 8), face_local_hi_scoria, scoria_n)
        call scoria_write_1(cell_val_mnmx_hilo_mi_2_scoria, cell_val_mnmx_hilo_mi_2_packed, &
        & int(size(cell_val_mnmx_hilo_mi_2_scoria), 8), face_local_hi_scoria, scoria_n)

        call scoria_write_1(cell_val_mnmx_hilo_mo_3_scoria, cell_val_mnmx_hilo_mo_3_packed, &
        & int(size(cell_val_mnmx_hilo_mo_3_scoria), 8), face_local_lo_scoria, scoria_n)
        call scoria_write_1(cell_val_mnmx_hilo_mo_4_scoria, cell_val_mnmx_hilo_mo_4_packed, &
        & int(size(cell_val_mnmx_hilo_mo_4_scoria), 8), face_local_lo_scoria, scoria_n)

        TIMERSET(.false., scoria_inside_com3b)
        call system_clock(c_now, c_rate)
        end = real(c_now, REAL64) / c_rate
        my_dt = my_dt + (end - start)
        ! END SCORIA TIMING

        deallocate(face_local_hi_scoria)
        deallocate(face_local_lo_scoria)

        deallocate(cell_half_lo_scoria)
        deallocate(cell_value_hilo_1_scoria)
        deallocate(cell_value_hilo_2_scoria)
        deallocate(cell_half_hi_scoria)

        if (present(do_special)) then
          deallocate(rho_scoria)
        endif 

        deallocate(cell_val_mnmx_hilo_hi_1_scoria)
        deallocate(cell_val_mnmx_hilo_hi_2_scoria)
        deallocate(cell_val_mnmx_hilo_lo_3_scoria)
        deallocate(cell_val_mnmx_hilo_lo_4_scoria)
        deallocate(cell_val_mnmx_hilo_mi_1_scoria)
        deallocate(cell_val_mnmx_hilo_mi_2_scoria)
        deallocate(cell_val_mnmx_hilo_mo_3_scoria)
        deallocate(cell_val_mnmx_hilo_mo_4_scoria)

        deallocate(cell_half_lo_hi_packed) 
        deallocate(cell_value_hilo_lo_packed)
        deallocate(cell_half_hi_lo_packed)
        deallocate(cell_value_hilo_hi_packed)

        if (present(do_special)) then
          deallocate(rho_hi_packed)
          deallocate(rho_lo_packed)
        endif

        deallocate(cell_val_mnmx_hilo_hi_1_packed)
        deallocate(cell_val_mnmx_hilo_hi_2_packed)
        deallocate(cell_val_mnmx_hilo_lo_3_packed)
        deallocate(cell_val_mnmx_hilo_lo_4_packed)
        deallocate(cell_val_mnmx_hilo_mi_1_packed)
        deallocate(cell_val_mnmx_hilo_mi_2_packed)
        deallocate(cell_val_mnmx_hilo_mo_3_packed)
        deallocate(cell_val_mnmx_hilo_mo_4_packed)

        deallocate(face_local_hi_scoria_transformed)
        deallocate(face_local_lo_scoria_transformed)

        deallocate(range_face_local_hi_scoria)
        deallocate(range_face_local_lo_scoria)
      enddo ! loop

      TIMERSET(.false., inside_com3b)
      write(*,*) '          ', my_dt, trim("scoria_inside_com3b")
      
      !close(10)

      end associate
    end subroutine inside_com3b_scoria_ranged_read_b
end module my_scoria_ranged_read_b_derivatives