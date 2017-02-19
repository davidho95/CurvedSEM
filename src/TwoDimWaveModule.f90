module TwoDimWaveModule

  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: PI = 3.141592653589793

  integer, dimension(2), parameter :: NUM_SPEC_EL = (/2,2/)
  integer, parameter :: TOTAL_NUM_SPEC = NUM_SPEC_EL(1) * NUM_SPEC_EL(2)
  integer, parameter :: NUM_GLL = 3
  integer, dimension(2), parameter :: NUM_GLOBAL_POINTS = NUM_SPEC_EL * (NUM_GLL - 1) + (/1, 1/)
  integer :: total_num_glob

  character(len = 43) :: output_path = "/home/davidho/WaveEqCurvedSpacetime/output/"
  contains

  subroutine Setup_mesh_2d()

    implicit none

    real(dp), dimension(NUM_GLL) :: gll_points, gll_weights
    real(dp), dimension(NUM_GLL, NUM_GLL), save :: hprime

    integer, parameter :: TOTAL_NUM_NODES = NUM_GLL * NUM_GLL * TOTAL_NUM_SPEC
    real(dp), dimension(NUM_GLL, NUM_GLL, TOTAL_NUM_SPEC, 2) :: nodes
    integer, dimension(NUM_GLL, NUM_GLL, TOTAL_NUM_SPEC) :: i_bool
    real(dp), dimension(TOTAL_NUM_NODES, 2) :: temp_points
    integer, dimension(TOTAL_NUM_NODES) :: locval
    logical, dimension(TOTAL_NUM_NODES) :: ifseg
    real(dp), dimension(NUM_GLL, NUM_GLL, TOTAL_NUM_SPEC, 2, 2) :: metric_tensor, inv_metric
    real(dp), dimension(NUM_GLL, NUM_GLL, TOTAL_NUM_SPEC) :: rho, mu, metric_det
    real(dp), dimension(:, :), allocatable :: global_points
    real(dp), dimension(:, :), allocatable :: torus_points
    real(dp), dimension(:), allocatable :: mass_mat

    real(dp), external :: lagrange_deriv_GLL

    integer i_gll, j_gll, i_spec, i_spec_x, i_spec_y, i_eoff, i_loc, i_glob

    ! get the GLL points and weights
    call zwgljd(gll_points,gll_weights,NUM_GLL,0.0_dp,0.0_dp)
    if(mod(NUM_GLL,2) /= 0) gll_points((NUM_GLL-1)/2+1) = 0.0_dp

    ! get the derivatives of the Lagrange polynomials at 
    ! the GLL points; recall that  hprime(i,j)=h'_{j}(xigll_{i}) 
    do j_gll = 1,NUM_GLL
      do i_gll=1,NUM_GLL
        hprime(i_gll, j_gll) = lagrange_deriv_GLL(j_gll - 1,i_gll - 1,gll_points,NUM_GLL)
      end do
    end do

    ! Setup nodes
    i_spec = 0
    do i_spec_x = 1, NUM_SPEC_EL(1)
      do i_spec_y = 1, NUM_SPEC_EL(2)
        i_spec = i_spec + 1
        do i_gll = 1, NUM_GLL
           do j_gll = 1, NUM_GLL
              nodes(i_gll,j_gll,i_spec, 1) =&
                (dble(i_spec_x) - 0.5_dp + gll_points(i_gll)*0.5_dp)/NUM_SPEC_EL(1)
              nodes(i_gll,j_gll,i_spec, 2) =&
                (dble(i_spec_y) - 0.5_dp + gll_points(j_gll)*0.5_dp)/NUM_SPEC_EL(2)
              rho(i_gll,j_gll,i_spec) = density_fn(nodes(i_gll, j_gll, i_spec, :))
              mu(i_gll,j_gll,i_spec) = rigidity_fn(nodes(i_gll, j_gll, i_spec, :))
           enddo
        enddo
      enddo
    enddo

    metric_tensor = compute_metric(nodes, hprime)
    metric_det(:, :, :) = metric_tensor(:, :, :, 1, 1) * metric_tensor(:, :, :, 2, 2)&
      - metric_tensor(:, :, :, 1, 2) * metric_tensor(:, :, :, 2, 1)

    do i_spec = 1, TOTAL_NUM_SPEC
       i_eoff = NUM_GLL*NUM_GLL*(i_spec-1)
       i_loc = 0
       do j_gll = 1,NUM_GLL
          do i_gll = 1,NUM_GLL
             i_loc = i_loc + 1
             temp_points(i_loc+i_eoff, :) = nodes(i_gll,j_gll,i_spec, :)
          enddo
       enddo
    enddo

    locval = 0.
    ifseg = .FALSE.

    call get_global(NUM_GLL,TOTAL_NUM_SPEC,temp_points(:, 1),temp_points(:, 2),&
      i_bool,locval,ifseg,total_num_glob,TOTAL_NUM_NODES)

    allocate(torus_points(total_num_glob, 3))

    torus_points = nodes_to_torus(nodes, i_bool, 2d0, 1d0)

    open(unit=10, file=output_path//"mesh.in", action="write", status="unknown")
    do i_spec = 1, TOTAL_NUM_SPEC
      do i_gll = 1, NUM_GLL
        do j_gll = 1, NUM_GLL
        write(10,*) nodes(i_gll, j_gll, i_spec, :),&
          metric_det(i_gll, j_gll, i_spec),&
          metric_tensor(i_gll, j_gll, i_spec, :, :),&
          rho(i_gll, j_gll, i_spec),&
          mu(i_gll, j_gll, i_spec)
        enddo
      enddo
    enddo

    allocate(mass_mat(total_num_glob))

    inv_metric(:,:,:,1,1) = metric_tensor(:,:,:,2,2) / metric_det
    inv_metric(:,:,:,1,2) = -metric_tensor(:,:,:,1,2) / metric_det
    inv_metric(:,:,:,2,1) = -metric_tensor(:,:,:,2,1) / metric_det
    inv_metric(:,:,:,2,2) = metric_tensor(:,:,:,1,1) / metric_det

    call increment_system(rho, mu, metric_det, inv_metric, i_bool, hprime, gll_weights)

  end subroutine Setup_mesh_2d

  subroutine increment_system(rho, mu, metric_det, inv_metric, i_bool, hprime, gll_weights)
    real(dp) rho(:, :, :), mu(:, :, :), inv_metric(:, :, :, :, :), metric_det(:, :, :)
    real(dp) hprime(:,:), gll_weights(:)
    integer i_bool(:, :, :)
    real(dp), dimension(total_num_glob) :: force_vec

    real(dp), dimension(25) :: displ = (/3.7266497013639212d-06,&
       3.7266497013639212d-06, 3.7266497013639212d-06, 3.7266497013639212d-06,&
       3.7266497013639212d-06, 6.1019324154531507d-13, 6.1019324154531507d-13,&
       6.1019324154531507d-13,  6.1019324154531507d-13, 6.1019324154531507d-13,&
       0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0/)

    integer i_spec, i_gll, j_gll

    force_vec = force_vector(displ, mu, metric_det, inv_metric, hprime, gll_weights, i_bool)

  end subroutine increment_system

  function mass_mat_glob(rho, metric_det, gll_weights, i_bool)
    real(dp) rho(:, :, :), metric_det(:, :, :), gll_weights(:)
    integer i_bool(:, :, :)
    integer i_spec, i_gll, j_gll, i_glob
    real(dp) mass_mat_glob(total_num_glob)

    mass_mat_glob(:) = 0d0

    do i_spec = 1,TOTAL_NUM_SPEC
       do i_gll = 1,NUM_GLL
          do j_gll = 1,NUM_GLL
             i_glob = i_bool(i_gll,j_gll,i_spec)
             mass_mat_glob(i_glob) = mass_mat_glob(i_glob) + gll_weights(i_gll)*gll_weights(j_gll) &
                                     *sqrt(metric_det(i_gll,j_gll,i_spec)) &
                                     *rho(i_gll,j_gll,i_spec)
          enddo
       enddo
    enddo
  end function mass_mat_glob

  function force_vector(u, mu, metric_det, inv_metric, hprime, gll_weights, i_bool)
    real(dp), dimension(:, :, :) :: mu, metric_det
    real(dp), dimension(:, :, :, :, :) :: inv_metric
    real(dp), dimension(:) :: gll_weights
    real(dp), dimension(:, :) :: hprime
    real(dp), dimension(NUM_GLL, NUM_GLL, 2) :: derivative_vec, temp
    real(dp), dimension(NUM_GLL, NUM_GLL) :: u_loc, temp2, gll_weights_mat, force_vec_loc
    integer i_spec, i_gll, j_gll, k_gll
    real(dp) u(total_num_glob)
    integer i_bool(:,:,:)
    real(dp) force_vector(total_num_glob)

    force_vector = 0d0

    do i_spec = 1, TOTAL_NUM_SPEC

      do i_gll = 1, NUM_GLL
        do j_gll = 1, NUM_GLL
          i_glob = i_bool(i_gll, j_gll, i_spec)
          u_loc(i_gll, j_gll) = u(i_glob)
        enddo
      enddo

      derivative_vec(:, :, 1) = matmul(hprime, u_loc) + matmul(u_loc, transpose(hprime))
      derivative_vec(:, :, 2) = matmul(hprime, u_loc) + matmul(u_loc, transpose(hprime))

      do i_gll = 1, NUM_GLL
        do j_gll = 1, NUM_GLL
          i_glob = i_bool(i_gll, j_gll, i_spec)
          temp(i_gll, j_gll, :) = matmul(inv_metric(i_gll,j_gll,i_spec, :, :), derivative_vec(i_gll, j_gll, :))&
            * sqrt(metric_det(i_gll, j_gll, i_spec)) * mu(i_gll, j_gll, i_spec)
        enddo
      enddo

      do i_gll = 1, NUM_GLL
        do j_gll = 1, NUM_GLL
          gll_weights_mat(i_gll, j_gll) = gll_weights(i_gll) * gll_weights(j_gll)
          temp2(i_gll, j_gll) = (temp(i_gll, j_gll, 1) + temp(i_gll, j_gll, 2))*hprime(i_gll, j_gll)
        enddo
      enddo

      force_vec_loc = -transpose(matmul(gll_weights_mat, temp2) / 2)

      do i_gll = 1, NUM_GLL
        do j_gll = 1, NUM_GLL
          i_glob = i_bool(i_gll, j_gll, i_spec)
          force_vector(i_glob) = force_vector(i_glob) + force_vec_loc(i_gll, j_gll)
        enddo
      enddo
    enddo
    print *, force_vector
    end function force_vector

  function density_fn(position_vec)
    real(dp) position_vec(2)
    real(dp) density_fn
    integer i_pos

    density_fn = 1d0
  end function density_fn

  function rigidity_fn(position_vec)
    real(dp) position_vec(2)
    real(dp) rigidity_fn
    integer i_pos

    rigidity_fn = 1d0
  end function rigidity_fn

  function compute_metric(nodes, hprime)

    implicit none

    integer i_gll, j_gll, k_gll, i_spec
    real(dp) nodes(NUM_GLL, NUM_GLL, TOTAL_NUM_SPEC, 2)
    real(dp) hprime(NUM_GLL, NUM_GLL)
    real(dp) jacobian(NUM_GLL, NUM_GLL, TOTAL_NUM_SPEC, 2, 2)
    real(dp) compute_metric(NUM_GLL, NUM_GLL, TOTAL_NUM_SPEC, 2, 2)

    jacobian = 0d0

    do i_spec = 1, TOTAL_NUM_SPEC
      jacobian(:, :, i_spec, 1, 1) = matmul(hprime, nodes(:, :, i_spec, 1))
      jacobian(:, :, i_spec, 1, 2) = matmul(nodes(:, :, i_spec, 1), transpose(hprime))
      jacobian(:, :, i_spec, 2, 1) = matmul(hprime, nodes(:, :, i_spec, 2))
      jacobian(:, :, i_spec, 2, 2) = matmul(nodes(:, :, i_spec, 2), transpose(hprime))
    enddo

    do i_spec = 1, TOTAL_NUM_SPEC
      do i_gll = 1, NUM_GLL
        do j_gll = 1, NUM_GLL
          compute_metric(i_gll, j_gll, i_spec, :, :) =&
            matmul(transpose(jacobian(i_gll, j_gll, i_spec, :, :)),&
            jacobian(i_gll, j_gll, i_spec, :, :))
        enddo
      enddo
    enddo
  end function compute_metric

  function nodes_to_torus(nodes, i_bool, r1, r2)

    implicit none

    real(dp) nodes(NUM_GLL, NUM_GLL, TOTAL_NUM_SPEC, 2)
    real(dp) r1, r2
    integer i_bool(NUM_GLL, NUM_GLL, TOTAL_NUM_SPEC) 
    integer i_spec, i_gll, j_gll, i_glob
    real(dp) nodes_to_torus(total_num_glob, 3)

    do i_spec = 1, TOTAL_NUM_SPEC
      do j_gll = 1, NUM_GLL
        do i_gll = 1, NUM_GLL
          i_glob = i_bool(i_gll, j_gll, i_spec)
            nodes_to_torus(i_glob, 1) = (r1 + r2 * cos(nodes(i_gll, j_gll, i_spec, 1)))&
             * cos(nodes(i_gll, j_gll, i_spec, 2))
            nodes_to_torus(i_glob, 2) = (r1 + r2 * cos(nodes(i_gll, j_gll, i_spec, 1)))&
             * sin(nodes(i_gll, j_gll, i_spec, 2))
            nodes_to_torus(i_glob, 3) = r2 * sin(nodes(i_gll, j_gll, i_spec, 1))
        enddo
      enddo
    enddo
  end function nodes_to_torus

end module TwoDimWaveModule