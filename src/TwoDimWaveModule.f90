module TwoDimWaveModule

  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: PI = 3.141592653589793

  integer, dimension(2), parameter :: NUM_SPEC_EL = (/2, 2/)
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
    real(dp), dimension(NUM_GLL, NUM_GLL, TOTAL_NUM_SPEC, 2, 2) :: metric_tensor
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

    mass_mat = mass_mat_glob(rho, metric_det, i_bool, gll_weights)
    print *, mass_mat

  end subroutine Setup_mesh_2d

  function mass_mat_glob(rho, metric_det, i_bool, gll_weights)
    real(dp) rho(:, :, :), metric_det(:, :, :), gll_weights(:)
    integer i_bool(:, :, :)
    integer i_spec, i_gll, j_gll, i_glob
    real(dp) mass_mat_glob(total_num_glob)

    mass_mat_glob(:) = 0d0

    do i_spec = 1,TOTAL_NUM_SPEC
       do i_gll = 1,NUM_GLL
          do j_ll = 1,NUM_GLL
             i_glob = i_bool(i_gll,j_ll,i_spec)
             mass_mat_glob(i_glob) = mass_mat_glob(i_glob) + gll_weights(i_gll)*gll_weights(j_ll) &
                                     *sqrt(metric_det(i_gll,j_ll,i_spec)) &
                                     *rho(i_gll,j_ll,i_spec)
          enddo
       enddo
    enddo
  end function mass_mat_glob

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

    ! do i_spec = 1, TOTAL_NUM_SPEC
    !   do j_gll = 1, NUM_GLL
    !     do i_gll = 1, NUM_GLL
    !       do k_gll = 1, NUM_GLL
    !         jacobian(i_gll, j_gll, i_spec, 1, 1) = jacobian(i_gll, j_gll, i_spec, 1, 1)+&
    !          nodes(k_gll,j_gll,i_spec, 1)*hprime(i_gll,k_gll)
    !         jacobian(i_gll, j_gll, i_spec, 1, 2) = jacobian(i_gll, j_gll, i_spec, 1, 2)+&
    !          nodes(i_gll,k_gll,i_spec, 1)*hprime(j_gll,k_gll)
    !         jacobian(i_gll, j_gll, i_spec, 2, 1) = jacobian(i_gll, j_gll, i_spec, 2, 1)+&
    !          nodes(k_gll,j_gll,i_spec, 2)*hprime(i_gll,k_gll)
    !         jacobian(i_gll, j_gll, i_spec, 2, 2) = jacobian(i_gll, j_gll, i_spec, 2, 2)+&
    !          nodes(i_gll,k_gll,i_spec, 2)*hprime(j_gll,k_gll)
    !       enddo
    !     enddo
    !   enddo
    ! enddo

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