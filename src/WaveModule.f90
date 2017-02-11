module WaveModule

  use Parameters

  implicit none

  contains
!===============================================================================
! Subroutine to set up the global mesh and local ot global numbering
!===============================================================================
  subroutine Setup_mesh(i_bool, rho, mu, global_points, gll_weights, h_prime, jacobian_mat, jacobian)
    integer i_spec, i_gll, i_glob, i_h_prime, j_h_prime
    real(dp), dimension(0:NUM_SPEC_EL) :: anchor_points
    real(dp), dimension(NUM_GLL) :: gll_points, gll_weights
    real(dp), dimension(NUM_GLOBAL_POINTS) :: global_points, rho, mu
    integer, dimension(NUM_GLL, NUM_SPEC_EL) :: i_bool
    real(dp), dimension(NUM_GLL, NUM_SPEC_EL) :: jacobian_mat, jacobian
    real(dp), dimension(NUM_GLL, NUM_GLL) :: h_prime
    real(dp), external :: lagrange_deriv_GLL

    ! GLL points and weights
    call zwgljd(gll_points,gll_weights,NUM_GLL,0.0_dp,0.0_dp)
    if(mod(NUM_GLL,2) /= 0) gll_points((NUM_GLL-1)/2+1) = 0.0_dp

    ! get the derivatives of the Lagrange polynomials at 
    ! the GLL points; recall that  h_prime(i,j)=h'_{j}(xigll_{i}) 
    do i_h_prime=1,NUM_GLL
       do j_h_prime=1,NUM_GLL
          h_prime(i_h_prime,j_h_prime) = lagrange_deriv_GLL(j_h_prime-1,i_h_prime-1,gll_points,NUM_GLL)
       end do
    end do

    ! Setup anchor points
    do i_spec = 0, NUM_SPEC_EL
      anchor_points(i_spec) = dble(i_spec) * LENGTH / dble(NUM_SPEC_EL)
    enddo

    do i_spec = 1,NUM_SPEC_EL
      do i_gll = 1,NUM_GLL
        jacobian_mat(i_gll,i_spec) = 2. / (anchor_points(i_spec)-anchor_points(i_spec-1)) ! This is d(xi) / dx
        jacobian(i_gll,i_spec) = (anchor_points(i_spec)-anchor_points(i_spec-1)) / 2.
      enddo
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
        global_points(i_glob) = 0.5*(1.-gll_points(i_gll))*anchor_points(i_spec-1)&
          +0.5*(1.+gll_points(i_gll))*anchor_points(i_spec)
      enddo
    enddo

    ! Set mesh properties
    rho = density_fn(global_points)
    mu = rigidity_fn(global_points)

  end subroutine Setup_mesh

  function calculate_delta_t(rho, mu)
    real(dp) rho(NUM_GLOBAL_POINTS)
    real(dp) mu(NUM_GLOBAL_POINTS)
    real(dp) calculate_delta_t, delta_h, max_vel
    real(dp) :: courant_CFL = 0.4d0 !CFL stability condition

    delta_h = LENGTH / dble(NUM_GLOBAL_POINTS)
    max_vel = maxval(sqrt(mu / rho))
    calculate_delta_t = courant_CFL * delta_h / max_vel
  end function calculate_delta_t 

  function mass_mat_glob(rho, jacobian, i_bool, gll_weights)
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
  end function mass_mat_glob

  subroutine Increment_system(displ, vel, accel, rho, mu, mass_mat, i_bool, h_prime,&
   gll_weights, jacobian_mat, jacobian, delta_t)
    real(dp) displ(:), vel(:), accel(:), rho(:), mu(:), mass_mat(:), h_prime(:, :),&
     gll_weights(:), jacobian_mat(:,:), jacobian(:,:)
    integer i_bool(:, :)
    real(dp) delta_t, du_dxi, epsilon, sigma
    integer i_spec, i_gll, j_gll, k_gll, i_glob
    real(dp) temp(NUM_GLL)
    real(dp) templ
    real(dp) boundary_displ, boundary_vel, boundary_accel

    ! `Predictor' update displacement using explicit finite-difference time scheme (Newmark)
      displ(:) = displ(:) + delta_t * vel(:) + delta_t**2/2 * accel(:)
      vel(:) = vel(:) + delta_t/2  *accel(:)
      accel(:) = 0.
      do i_spec = 1,NUM_SPEC_EL
        do i_gll = 1,NUM_GLL
          ! Compute d(u) / d(xi)
          du_dxi = 0.
          do j_gll = 1,NUM_GLL
            i_glob = i_bool(j_gll,i_spec)
            du_dxi = du_dxi + displ(i_glob)*h_prime(i_gll,j_gll)
          enddo
          i_glob = i_bool(i_gll, i_spec)
          ! Strain (i.e., d(u) / dx in 1D)
          epsilon = du_dxi*jacobian_mat(i_gll,i_spec)
          ! stress
          sigma = mu(i_glob)*epsilon
          temp(i_gll) = jacobian(i_gll,i_spec)*sigma*jacobian_mat(i_gll,i_spec)
        enddo ! first loop over the GLL points
        do k_gll = 1,NUM_GLL
          templ = 0.
          do i_gll = 1,NUM_GLL
            templ = templ + temp(i_gll)*h_prime(i_gll,k_gll)*gll_weights(i_gll)
          enddo
          ! `Corrector' update of acceleration in the Newmark scheme
          ! The minus sign comes from the integration by part done in the weak formulation of the equations
          i_glob = i_bool(k_gll,i_spec)
          accel(i_glob) = accel(i_glob) - templ
        enddo ! Second loop over the GLL points
      enddo ! End loop over all spectral elements
      ! Periodic BCs
      boundary_displ = (displ(1) + displ(NUM_GLOBAL_POINTS)) / 2
      boundary_vel = (vel(1) + vel(NUM_GLOBAL_POINTS)) / 2
      boundary_accel = (accel(1) + accel(NUM_GLOBAL_POINTS)) / 2
      displ(1) = boundary_displ
      displ(NUM_GLOBAL_POINTS) = boundary_displ
      vel(1) = boundary_vel
      vel(NUM_GLOBAL_POINTS) = boundary_vel
      accel(1) = boundary_accel
      accel(NUM_GLOBAL_POINTS) = boundary_accel

      ! Fixed BCs
      ! accel(1) = 0
      ! accel(NUM_GLOBAL_POINTS) = 0

      ! Divide by the mass matrix, which is strictly (i.e. perfectly) diagonal
      accel(:) = accel(:)/mass_mat(:)
      ! `Corrector' update velocity
      vel(:) = vel(:) + delta_t/2 * accel(:)
  end subroutine Increment_system

  subroutine Output_snapshot(snapshot_dir, global_points, displ, timestep)
    character(len=*) snapshot_dir
    real(dp) global_points(:), displ(:)
    integer timestep, i_glob
    character(len=50) snapshot_file

    write(snapshot_file, "('snapshot',i5.5)") timestep

    open(unit=10, file=snapshot_dir//snapshot_file, action="write", status="unknown")
    do i_glob = 1, NUM_GLOBAL_POINTS
      write(10,*) sngl(global_points(i_glob)),sngl(displ(i_glob))
    enddo
  end subroutine Output_snapshot

end module WaveModule
