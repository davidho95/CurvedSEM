module WaveModule

  use MeshClass

  contains

  subroutine increment_system(a_mesh, delta_t)

    implicit none

    type(Mesh) :: a_mesh
    real(dp) delta_t

    real(dp), allocatable :: force_vec(:), mass_mat(:)

    allocate(force_vec(a_mesh%total_num_glob))
    allocate(mass_mat(a_mesh%total_num_glob))

    mass_mat = mass_mat_glob(a_mesh)

    ! Use Newmark time-stepping method
    a_mesh%displ = a_mesh%displ + a_mesh%vel*delta_t + delta_t**2 / 2d0 * a_mesh%accel
    a_mesh%vel = a_mesh%vel + a_mesh%accel*delta_t / 2d0
    force_vec = force_vector(a_mesh)
    a_mesh%accel = force_vec / mass_mat
    a_mesh%vel = a_mesh%vel + a_mesh%accel*delta_t / 2d0

    a_mesh%displ = periodic_bc(a_mesh%displ, a_mesh)
    a_mesh%vel = periodic_bc(a_mesh%vel, a_mesh)
    a_mesh%accel = periodic_bc(a_mesh%accel, a_mesh)
  end subroutine increment_system

  function mass_mat_glob(a_mesh) result(mass_mat)

    implicit none

    type(mesh) a_mesh
    integer i_spec, i_gll, j_gll, i_glob
    real(dp), allocatable :: mass_mat(:)

    allocate(mass_mat(a_mesh%total_num_glob))

    mass_mat(:) = 0d0

    do i_spec = 1,a_mesh%total_num_spec
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

  function force_vector(a_mesh)

    implicit none

    type(Mesh) a_mesh
    real(dp), allocatable, dimension(:, :, :) :: derivative_vec, temp
    real(dp), allocatable, dimension(:, :) :: displ_loc, temp2, gll_weights_mat, force_vec_loc

    integer i_spec, i_gll, j_gll, k_gll, i_glob
    real(dp), allocatable :: force_vector(:)

    allocate(derivative_vec(a_mesh%num_gll, a_mesh%num_gll, 2))
    allocate(temp(a_mesh%num_gll, a_mesh%num_gll, 2))
    allocate(displ_loc(a_mesh%num_gll, a_mesh%num_gll))
    allocate(temp2(a_mesh%num_gll, a_mesh%num_gll))
    allocate(gll_weights_mat(a_mesh%num_gll, a_mesh%num_gll))
    allocate(force_vec_loc(a_mesh%num_gll, a_mesh%num_gll))
    allocate(force_vector(a_mesh%total_num_glob))

    force_vector = 0d0

    do i_spec = 1, a_mesh%total_num_spec

      do i_gll = 1, a_mesh%NUM_GLL
        do j_gll = 1, a_mesh%NUM_GLL
          i_glob = a_mesh%i_bool(i_gll, j_gll, i_spec)
          displ_loc(i_gll, j_gll) = a_mesh%displ(i_glob)
        enddo
      enddo

      derivative_vec(:, :, 1) = matmul(a_mesh%hprime, displ_loc)
      derivative_vec(:, :, 2) =  matmul(displ_loc, transpose(a_mesh%hprime))

      do i_gll = 1, a_mesh%NUM_GLL
        do j_gll = 1, a_mesh%NUM_GLL
          i_glob = a_mesh%i_bool(i_gll, j_gll, i_spec)
          temp(i_gll, j_gll, :) = matmul(a_mesh%inv_metric(i_gll,j_gll,i_spec, :, :), derivative_vec(i_gll, j_gll, :))&
            * sqrt(a_mesh%metric_det(i_gll, j_gll, i_spec))
          temp(i_gll, j_gll, :) = matmul(a_mesh%mu(i_gll, j_gll, i_spec, :, :), temp(i_gll, j_gll, :))
          temp(i_gll, j_gll, :) = a_mesh%gll_weights(i_gll) * a_mesh%gll_weights(j_gll) * temp(i_gll, j_gll, :)
        enddo
      enddo

      force_vec_loc = -matmul(transpose(a_mesh%hprime), temp(:, :, 1)) - matmul(temp(:, :, 2), a_mesh%hprime)

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
    real(dp) :: courant_CFL = 0.1d0
    real(dp), allocatable :: trace_mu(:, :, :)

    allocate(trace_mu(a_mesh%num_gll, a_mesh%num_gll, a_mesh%total_num_spec))

    num_global_points = (a_mesh%num_gll - 1) * a_mesh%num_spec_el + (/1, 1/)
    delta_h = 1d0 / dble(maxval(num_global_points))
    trace_mu = a_mesh%mu(:, :, :, 1, 1) + a_mesh%mu(:, :, :, 2, 2)
    max_vel = maxval(sqrt(trace_mu / a_mesh%rho))

    delta_t = courant_CFL * delta_h / max_vel

    end function estimate_timestep

    function periodic_bc(vector, a_mesh) result(periodic_vector)

    implicit none

    type(Mesh) a_mesh
    real(dp) vector(:)
    real(dp) boundary_val
    integer i_spec_x, i_spec_y, i_spec_y1, i_spec_y2,&
     i_spec_x1, i_spec_x2, i_gll, i_glob_y1, i_glob_y2,&
     i_glob_x1, i_glob_x2, i_glob
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

    do i_spec_y = 1, a_mesh%NUM_SPEC_EL(1)
       i_spec_x1 = i_spec_y
       i_spec_x2 = a_mesh%NUM_SPEC_EL(1) * (a_mesh%num_spec_el(2) - 1) + i_spec_y

      do i_gll = 1,a_mesh%NUM_GLL
        i_glob_x1 = a_mesh%i_bool(1,i_gll,i_spec_x1)
        i_glob_x2 = a_mesh%i_bool(a_mesh%num_gll,i_gll,i_spec_x2)

        ! set new u to be average of us on the boundaries
        boundary_val = (vector(i_glob_x1) + vector(i_glob_x2)) / 2
        periodic_vector(i_glob_x1) = boundary_val
        periodic_vector(i_glob_x2) = boundary_val
      end do
    end do
  end function periodic_bc

  subroutine initial_conditions(a_mesh)

    implicit none

    type(Mesh) a_mesh
    integer i_spec, i_gll, j_gll, i_glob

    a_mesh%displ = 0d0

    do i_spec = 1, a_mesh%total_num_spec
      do i_gll = 1, a_mesh%num_gll
        do j_gll = 1, a_mesh%num_gll
          i_glob = a_mesh%i_bool(i_gll, j_gll, i_spec)
          a_mesh%displ(i_glob) = 5 * exp(-(20*(a_mesh%nodes(i_gll, j_gll, i_spec, 2)-0.5d0))**2)!&
            !-(20*(a_mesh%nodes(i_gll, j_gll, i_spec, 2)-0.5d0))**2)
        enddo
      enddo
    enddo

    a_mesh%displ = periodic_bc(a_mesh%displ, a_mesh)

    a_mesh%vel = 0d0
    a_mesh%accel = periodic_bc(force_vector(a_mesh) / mass_mat_glob(a_mesh), a_mesh)
  end subroutine initial_conditions

  function displ_on_torus(a_mesh, magnitude) result(coordinates)

    implicit none

    type(Mesh) a_mesh
    real(dp) magnitude
    real(dp), allocatable :: coordinates(:,:)
    integer i_glob

    allocate(coordinates(a_mesh%total_num_glob, 3))

    do i_glob = 1, a_mesh%total_num_glob
      coordinates(i_glob, :) = a_mesh%torus_points(i_glob, :) + magnitude&
       * a_mesh%displ(i_glob) * a_mesh%torus_normal(i_glob, :)
    enddo

  end function displ_on_torus

  subroutine write_mesh(this)
    type(Mesh) :: this
    character(len=43) :: output_path = "/home/davidho/CurvedSEM/CurvedSEM2d/output/"
    integer i_gll, j_gll, i_spec

    open(unit=10, file=output_path//"mesh.in", action="write", status="unknown")
    do i_spec = 1, this%total_num_spec
      do i_gll = 1, this%num_gll
        do j_gll = 1, this%num_gll
        write(10,*) this%nodes(i_gll, j_gll, i_spec, :),&
          this%metric_det(i_gll, j_gll, i_spec),&
          this%inv_metric(i_gll, j_gll, i_spec, :, :),&
          this%rho(i_gll, j_gll, i_spec),&
          this%mu(i_gll, j_gll, i_spec, :, :)
        enddo
      enddo
    enddo

  end subroutine write_mesh

end module WaveModule