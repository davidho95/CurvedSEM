module MeshModule

  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: PI = 3.141592653589793

  type Mesh
    integer, dimension(2) :: num_spec_el
    integer :: num_gll
    integer :: total_num_glob
    real(dp), allocatable :: nodes(:, :, :, :)
    real(dp), allocatable :: inv_metric(:, :, :, :, :)
    real(dp), allocatable :: metric_det(:, :, :)
    real(dp), allocatable :: rho(:, :, :)
    real(dp), allocatable :: mu(:, :, :)
    integer, allocatable :: i_bool(:, :, :)
    real(dp), allocatable :: gll_weights(:)
    real(dp), allocatable :: hprime(:,:)
    real(dp), allocatable :: torus_points(:, :)
    real(dp), allocatable :: torus_normal(:, :)
    real(dp) :: r1
    real(dp) :: r2
  end type Mesh
  contains

  subroutine initialise_mesh(this, num_spec_el, num_gll, r1, r2)

    implicit none

    type(Mesh) :: this
    integer, dimension(2) :: num_spec_el
    integer :: num_gll
    real(dp) r1, r2

    real(dp), dimension(num_gll) :: gll_points, gll_weights
    real(dp), dimension(num_gll, num_gll) :: hprime

    integer :: total_num_spec, total_num_nodes
    real(dp), allocatable :: temp_points(:,:)
    integer, allocatable :: locval(:)
    logical, allocatable :: ifseg(:)
    real(dp), allocatable :: metric_tensor(:,:,:,:,:)
    real(dp), allocatable :: torus_points(:, :)

    real(dp), external :: lagrange_deriv_GLL

    integer i_gll, j_gll, i_spec, i_spec_x, i_spec_y, i_eoff, i_loc, i_glob

    total_num_spec = num_spec_el(1) * num_spec_el(2)
    total_num_nodes = num_gll * num_gll * total_num_spec

    this%num_spec_el = num_spec_el
    this%num_gll = num_gll

    this%r1 = r1
    this%r2 = r2

    allocate(this%nodes(num_gll, num_gll, total_num_spec, 2))
    allocate(this%inv_metric(num_gll, num_gll, total_num_spec, 2, 2))
    allocate(this%metric_det(num_gll, num_gll, total_num_spec))
    allocate(this%rho(num_gll, num_gll, total_num_spec))
    allocate(this%mu(num_gll, num_gll, total_num_spec))
    allocate(this%i_bool(num_gll, num_gll, total_num_spec))
    allocate(this%gll_weights(num_gll))
    allocate(this%hprime(num_gll, num_gll))

    allocate(temp_points(total_num_nodes, 2))
    allocate(locval(total_num_nodes))
    allocate(ifseg(total_num_nodes))
    allocate(metric_tensor(num_gll, num_gll, total_num_spec, 2, 2))

    ! get the GLL points and weights
    call zwgljd(gll_points,this%gll_weights,NUM_GLL,0.0_dp,0.0_dp)
    if(mod(NUM_GLL,2) /= 0) gll_points((NUM_GLL-1)/2+1) = 0.0_dp

    ! get the derivatives of the Lagrange polynomials at 
    ! the GLL points; recall that  hprime(i,j)=h'_{j}(xigll_{i}) 
    do j_gll = 1,NUM_GLL
      do i_gll=1,NUM_GLL
        this%hprime(i_gll, j_gll) = lagrange_deriv_GLL(j_gll - 1,i_gll - 1,gll_points,NUM_GLL)
      end do
    end do

    ! Setup nodes
    i_spec = 0
    do i_spec_x = 1, NUM_SPEC_EL(1)
      do i_spec_y = 1, NUM_SPEC_EL(2)
        i_spec = i_spec + 1
        do i_gll = 1, NUM_GLL
           do j_gll = 1, NUM_GLL
              this%nodes(i_gll,j_gll,i_spec, 1) =&
                (dble(i_spec_x) - 0.5_dp + gll_points(i_gll)*0.5_dp)/NUM_SPEC_EL(1)
              this%nodes(i_gll,j_gll,i_spec, 2) =&
                (dble(i_spec_y) - 0.5_dp + gll_points(j_gll)*0.5_dp)/NUM_SPEC_EL(2)
              this%rho(i_gll,j_gll,i_spec) = density_fn(this%nodes(i_gll, j_gll, i_spec, :))
              this%mu(i_gll,j_gll,i_spec) = rigidity_fn(this%nodes(i_gll, j_gll, i_spec, :))
           enddo
        enddo
      enddo
    enddo

    metric_tensor = compute_metric(this%nodes, this%hprime)

    this%metric_det(:, :, :) = metric_tensor(:, :, :, 1, 1) * metric_tensor(:, :, :, 2, 2)&
      - metric_tensor(:, :, :, 1, 2) * metric_tensor(:, :, :, 2, 1)

    this%inv_metric(:,:,:,1,1) = metric_tensor(:,:,:,2,2) / this%metric_det
    this%inv_metric(:,:,:,1,2) = -metric_tensor(:,:,:,1,2) / this%metric_det
    this%inv_metric(:,:,:,2,1) = -metric_tensor(:,:,:,2,1) / this%metric_det
    this%inv_metric(:,:,:,2,2) = metric_tensor(:,:,:,1,1) / this%metric_det

    do i_spec = 1, TOTAL_NUM_SPEC
       i_eoff = NUM_GLL*NUM_GLL*(i_spec-1)
       i_loc = 0
       do j_gll = 1,NUM_GLL
          do i_gll = 1,NUM_GLL
             i_loc = i_loc + 1
             temp_points(i_loc+i_eoff, :) = this%nodes(i_gll,j_gll,i_spec, :)
          enddo
       enddo
    enddo

    locval = 0.
    ifseg = .FALSE.

    call get_global(NUM_GLL,TOTAL_NUM_SPEC,temp_points(:, 1),temp_points(:, 2),&
      this%i_bool,locval,ifseg,this%total_num_glob,total_num_nodes)

    allocate(this%torus_points(this%total_num_glob, 3))
    allocate(this%torus_normal(this%total_num_glob, 3))

    this%torus_points = mesh_to_torus(this)
    this%torus_normal = calculate_torus_normal(this)

  end subroutine initialise_mesh

  function density_fn(position_vec)

    implicit none

    real(dp) position_vec(2)
    real(dp) density_fn
    integer i_pos

    density_fn = 1d0
  end function density_fn

  function rigidity_fn(position_vec)

    implicit none

    real(dp) position_vec(2)
    real(dp) rigidity_fn
    integer i_pos

    rigidity_fn = 1d0
  end function rigidity_fn

  function compute_metric(nodes, hprime) result(metric_tensor)

    implicit none

    integer i_gll, j_gll, k_gll, i_spec
    real(dp) nodes(:, :, :, :)
    real(dp) hprime(:, :)
    real(dp), allocatable :: jacobian(:, :, :, :, :)
    real(dp), allocatable :: metric_tensor(:, :, :, :, :)
    integer, dimension(4) :: nodes_shape 
    integer :: total_num_spec
    integer :: num_gll

    nodes_shape = shape(nodes)
    total_num_spec = nodes_shape(3)
    num_gll = nodes_shape(2)

    allocate(jacobian(num_gll, num_gll, total_num_spec, 2, 2))
    allocate(metric_tensor(num_gll, num_gll, total_num_spec, 2, 2))

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
          metric_tensor(i_gll, j_gll, i_spec, :, :) =&
            matmul(transpose(jacobian(i_gll, j_gll, i_spec, :, :)),&
            jacobian(i_gll, j_gll, i_spec, :, :))
        enddo
      enddo
    enddo
  end function compute_metric

  subroutine write_mesh(this)
    type(Mesh) :: this
    character(len=43) :: output_path = "/home/davidho/WaveEqCurvedSpacetime/output/"
    integer :: total_num_spec

    total_num_spec = this%num_spec_el(1) * this%num_spec_el(2)

    open(unit=10, file=output_path//"mesh.in", action="write", status="unknown")
    do i_spec = 1, total_num_spec
      do i_gll = 1, this%num_gll
        do j_gll = 1, this%num_gll
        write(10,*) this%nodes(i_gll, j_gll, i_spec, :),&
          this%metric_det(i_gll, j_gll, i_spec),&
          this%inv_metric(i_gll, j_gll, i_spec, :, :),&
          this%rho(i_gll, j_gll, i_spec),&
          this%mu(i_gll, j_gll, i_spec)
        enddo
      enddo
    enddo

  end subroutine write_mesh

  function mesh_to_torus(this) result(torus_points)

    implicit none

    type(Mesh) this
    integer total_num_spec
    integer i_spec, i_gll, j_gll, i_glob
    real(dp) :: tangent1(3), tangent2(3), normal(3)

    real(dp), allocatable :: torus_points(:, :)

    allocate(torus_points(this%total_num_glob, 3))
    total_num_spec = this%num_spec_el(1) * this%num_spec_el(2)

    do i_spec = 1, total_num_spec
      do i_gll = 1, this%num_gll
        do j_gll = 1, this%num_gll
          i_glob = this%i_bool(i_gll, j_gll, i_spec)
            torus_points(i_glob, 1) = (this%r1 + this%r2 * cos(2 * PI * this%nodes(i_gll, j_gll, i_spec, 1)))&
             * cos(2 * PI * this%nodes(i_gll, j_gll, i_spec, 2))
            torus_points(i_glob, 2) = (this%r1 + this%r2 * cos(2 * PI * this%nodes(i_gll, j_gll, i_spec, 1)))&
             * sin(2 * PI * this%nodes(i_gll, j_gll, i_spec, 2))
            torus_points(i_glob, 3) = this%r2 * sin(2 * PI * this%nodes(i_gll, j_gll, i_spec, 1))
        enddo
      enddo
    enddo
  end function mesh_to_torus

  function calculate_torus_normal(this) result(normal)

    implicit none

    type(Mesh) this
    real(dp), allocatable :: normal(:, :)
    real(dp) :: tangent1(3), tangent2(3)
    integer total_num_spec
    integer i_spec, i_gll, j_gll, i_glob

    allocate(normal(this%total_num_glob, 3))

    total_num_spec = this%num_spec_el(1)*this%num_spec_el(2)

    do i_spec = 1, total_num_spec
      do i_gll = 1, this%num_gll
        do j_gll = 1, this%num_gll
          i_glob = this%i_bool(i_gll, j_gll, i_spec)
          tangent1 = (/-sin(2 * PI * this%nodes(i_gll, j_gll, i_spec, 2)),&
           cos(2 * PI * this%nodes(i_gll, j_gll, i_spec, 2)), 0d0/)
          tangent2 = (/-cos(2 * PI * this%nodes(i_gll, j_gll, i_spec, 2))&
           * sin(2 * PI * this%nodes(i_gll, j_gll, i_spec, 1)),&
          -sin(2 * PI * this%nodes(i_gll, j_gll, i_spec, 2))&
           * sin(2 * PI * this%nodes(i_gll, j_gll, i_spec, 1)),&
           cos(2 * PI * this%nodes(i_gll, j_gll, i_spec, 1))/)

          normal(i_glob, :) = (/tangent1(2)*tangent2(3) - tangent1(3)*tangent2(2),&
                                tangent1(3)*tangent2(1) - tangent1(1)*tangent2(3),&
                                tangent1(1)*tangent2(2) - tangent1(2)*tangent2(1)/)
        enddo
      enddo
    enddo
  end function calculate_torus_normal

end module MeshModule