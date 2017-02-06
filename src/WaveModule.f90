module WaveModule

  use Parameters

  implicit none

  contains
!===============================================================================
! Subroutine to set up the global mesh and local ot global numbering
!===============================================================================
  subroutine Setup_mesh(i_bool, rho, mu)
    integer i_spec, i_gll, i_glob
    real(dp), dimension(0:NUM_SPEC_EL) :: anchor_points
    real(dp), dimension(NUM_GLL) :: gll_points, gll_weights
    real(dp), dimension(NUM_GLOBAL_POINTS) :: global_points, rho, mu
    integer, dimension(NUM_GLL, NUM_SPEC_EL) :: i_bool

    ! GLL points and weights
    call zwgljd(gll_points,gll_weights,NUM_GLL,0.0_dp,0.0_dp)
    if(mod(NUM_GLL,2) /= 0) gll_points((NUM_GLL-1)/2+1) = 0.0_dp

    ! Setup anchor points
    do i_spec = 0, NUM_SPEC_EL
      anchor_points(i_spec) = dble(i_spec) * LENGTH / dble(NUM_SPEC_EL)
    enddo

    ! Local to global numbering
    i_glob = 1
    do i_spec = 1,NUM_SPEC_EL
      do i_gll = 1,NUM_GLL
        if(i_gll > 1) i_glob = i_glob+1
        i_bool(i_gll,i_spec) = i_glob
      enddo
    enddo

  ! Compute global gridpoint position
    do i_spec = 1,NUM_SPEC_EL
      do i_gll = 1,NUM_GLL
        i_glob = i_bool(i_gll,i_spec)
        global_points(i_glob) = 0.5*(1.-gll_points(i_gll))*anchor_points(i_spec)&
          +0.5*(1.+gll_points(i_gll))*anchor_points(i_spec - 1)
      enddo
    enddo

    ! Set mesh properties
    rho = density_fn(global_points)
    mu = rigidity_fn(global_points)

  end subroutine Setup_mesh

  subroutine Assemble_global_mass_matrix(rho, jacobian, i_bool, gll_weights)
    real(dp) rho(:), jacobian(:, :), gll_weights(:)
    integer i_bool(:, :)
    integer i_spec, i_gll, i_glob
    real(dp) mass_mat_local
    real(dp) mass_mat_glob(NUM_GLOBAL_POINTS)
    do i_spec = 1, size(i_bool, 2)
      do i_gll = 1, size(i_bool, 1)
        i_glob = i_bool(i_gll,i_spec)
        mass_mat_local = gll_weights(i_gll)*(rho(i_glob))*jacobian(i_gll,i_spec)
        mass_mat_glob(i_glob) = mass_mat_glob(i_glob) + mass_mat_local
      enddo
    enddo
  end subroutine Assemble_global_mass_matrix

  subroutine Increment_system()

  end subroutine Increment_system
end module WaveModule
