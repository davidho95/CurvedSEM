module WaveModule

  use MeshModule

  contains

  subroutine increment_system(a_mesh, displ, vel, accel, delta_t)

    implicit none

    type(Mesh) :: a_mesh
    real(dp) displ(:), vel(:), accel(:)
    real(dp) delta_t

    real(dp), allocatable :: force_vec(:), mass_mat(:)

    allocate(force_vec(size(accel)))
    allocate(mass_mat(size(accel)))

    mass_mat = mass_mat_glob(a_mesh)

    ! Use Newmark time-stepping method
    displ = displ + vel*delta_t + delta_t**2 / 2 * accel
    vel = vel + accel*delta_t / 2
    force_vec = force_vector(displ, a_mesh)
    accel = force_vec / mass_mat
    vel = vel + accel*delta_t / 2

    displ = periodic_bc(displ, a_mesh)
    vel = periodic_bc(vel, a_mesh)
    accel = periodic_bc(accel, a_mesh)
  end subroutine increment_system

  function mass_mat_glob(a_mesh) result(mass_mat)

    implicit none

    type(mesh) a_mesh
    integer total_num_spec
    integer i_spec, i_gll, j_gll, i_glob
    real(dp), allocatable :: mass_mat(:)

    total_num_spec = a_mesh%num_spec_el(1) * a_mesh%num_spec_el(2)

    allocate(mass_mat(a_mesh%total_num_glob))

    mass_mat(:) = 0d0

    do i_spec = 1,TOTAL_NUM_SPEC
       do i_gll = 1,a_mesh%NUM_GLL
          do j_gll = 1,a_mesh%NUM_GLL
             i_glob = a_mesh%i_bool(i_gll,j_gll,i_spec)
             mass_mat(i_glob) = mass_mat(i_glob) + a_mesh%gll_weights(i_gll)*a_mesh%gll_weights(j_gll) &
                                     *sqrt(a_mesh%metric_det(i_gll,j_gll,i_spec)) &
                                     *a_mesh%rho(i_gll,j_gll,i_spec)
          enddo
       enddo
    enddo
  end function mass_mat_glob

  function force_vector(displ, a_mesh)

    implicit none

    type(Mesh) a_mesh
    real(dp) displ(:)
    real(dp), allocatable, dimension(:, :, :) :: derivative_vec, temp
    real(dp), allocatable, dimension(:, :) :: displ_loc, temp2, gll_weights_mat, force_vec_loc
    integer total_num_spec

    integer i_spec, i_gll, j_gll, k_gll, i_glob
    real(dp), allocatable :: force_vector(:)

    total_num_spec = a_mesh%num_spec_el(1) * a_mesh%num_spec_el(2)

    allocate(derivative_vec(a_mesh%num_gll, a_mesh%num_gll, 2))
    allocate(temp(a_mesh%num_gll, a_mesh%num_gll, 2))
    allocate(displ_loc(a_mesh%num_gll, a_mesh%num_gll))
    allocate(temp2(a_mesh%num_gll, a_mesh%num_gll))
    allocate(gll_weights_mat(a_mesh%num_gll, a_mesh%num_gll))
    allocate(force_vec_loc(a_mesh%num_gll, a_mesh%num_gll))
    allocate(force_vector(a_mesh%total_num_glob))

    force_vector = 0d0

    do i_spec = 1, TOTAL_NUM_SPEC

      do i_gll = 1, a_mesh%NUM_GLL
        do j_gll = 1, a_mesh%NUM_GLL
          i_glob = a_mesh%i_bool(i_gll, j_gll, i_spec)
          displ_loc(i_gll, j_gll) = displ(i_glob)
        enddo
      enddo

      derivative_vec(:, :, 1) = matmul(a_mesh%hprime, displ_loc) + matmul(displ_loc, transpose(a_mesh%hprime))
      derivative_vec(:, :, 2) = matmul(a_mesh%hprime, displ_loc) + matmul(displ_loc, transpose(a_mesh%hprime))

      do i_gll = 1, a_mesh%NUM_GLL
        do j_gll = 1, a_mesh%NUM_GLL
          i_glob = a_mesh%i_bool(i_gll, j_gll, i_spec)
          temp(i_gll, j_gll, :) = matmul(a_mesh%inv_metric(i_gll,j_gll,i_spec, :, :), derivative_vec(i_gll, j_gll, :))&
            * sqrt(a_mesh%metric_det(i_gll, j_gll, i_spec)) * a_mesh%mu(i_gll, j_gll, i_spec)
        enddo
      enddo

      do i_gll = 1, a_mesh%NUM_GLL
        do j_gll = 1, a_mesh%NUM_GLL
          gll_weights_mat(i_gll, j_gll) = a_mesh%gll_weights(i_gll) * a_mesh%gll_weights(j_gll)
          temp2(i_gll, j_gll) = (temp(i_gll, j_gll, 1) + temp(i_gll, j_gll, 2))*a_mesh%hprime(i_gll, j_gll)
        enddo
      enddo

      force_vec_loc = -transpose(matmul(gll_weights_mat, temp2) / 2)

      do i_gll = 1, a_mesh%NUM_GLL
        do j_gll = 1, a_mesh%NUM_GLL
          i_glob = a_mesh%i_bool(i_gll, j_gll, i_spec)
          force_vector(i_glob) = force_vector(i_glob) + force_vec_loc(i_gll, j_gll)
        enddo
      enddo
    enddo
  end function force_vector

  function estimate_timestep(a_mesh) result(delta_t)

    implicit none

    type(Mesh) a_mesh
    integer num_global_points(2)
    real(dp) max_vel, delta_h, delta_t
    real(dp) :: courant_CFL = 0.4d0

    num_global_points = (a_mesh%num_gll - 1) * a_mesh%num_spec_el + (/1, 1/)
    delta_h = 1d0 / dble(maxval(num_global_points))
    max_vel = maxval(sqrt(a_mesh%mu / a_mesh%rho))

    delta_t = courant_CFL * delta_h / max_vel

    end function estimate_timestep

    function periodic_bc(vector, a_mesh) result(periodic_vector)

    implicit none

    type(Mesh) a_mesh
    real(dp) vector(:)
    real(dp) boundary_val
    integer i_spec_x, i_spec_y1, i_spec_y2, i_gll, i_glob_y1, i_glob_y2
    real(dp), allocatable :: periodic_vector(:)

    allocate(periodic_vector(a_mesh%total_num_glob))

    periodic_vector = vector

    do i_spec_x = 1, a_mesh%NUM_SPEC_EL(1)
       i_spec_y1 = (i_spec_x - 1) * a_mesh%NUM_SPEC_EL(2) + 1
       i_spec_y2 = i_spec_x * a_mesh%NUM_SPEC_EL(2)

    do i_gll = 1,a_mesh%NUM_GLL
        i_glob_y1 = a_mesh%i_bool(i_gll,1,i_spec_y1)
        i_glob_y2 = a_mesh%i_bool(i_gll,a_mesh%NUM_GLL,i_spec_y2)

        ! set new u to be average of us on the boundaries
        boundary_val = (vector(i_glob_y1) + vector(i_glob_y2)) / 2
        periodic_vector(i_glob_y1) = boundary_val
        periodic_vector(i_glob_y2) = boundary_val
      end do
    end do
  end function periodic_bc

  subroutine initial_conditions(a_mesh, displ, vel, accel)

    implicit none

    type(Mesh) a_mesh
    real(dp), allocatable :: displ(:), vel(:), accel(:)
    integer total_num_spec
    integer i_spec, i_gll, j_gll, i_glob

    total_num_spec = a_mesh%num_spec_el(1) * a_mesh%num_spec_el(2)

    allocate(displ(a_mesh%total_num_glob))
    allocate(vel(a_mesh%total_num_glob))
    allocate(accel(a_mesh%total_num_glob))

    vel = 1d0
    accel = 1d0

    do i_spec = 1, total_num_spec
      do i_gll = 1, a_mesh%num_gll
        do j_gll = 1, a_mesh%num_gll
          i_glob = a_mesh%i_bool(i_gll, j_gll, i_spec)
          displ(i_glob) = sin(2 * PI * a_mesh%nodes(i_gll, j_gll, i_spec, 1))
        enddo
      enddo
    enddo
  end subroutine initial_conditions

end module WaveModule