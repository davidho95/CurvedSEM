module OneDimWaveModule

  implicit none

  integer, parameter :: dp = kind(1.0d0)

! System parameters
  character(len = 43) :: output_path = "/home/davidho/WaveEqCurvedSpacetime/output/"

! Mathematical constants
  real(dp), parameter :: PI = 3.141592653589793d0


  type Mesh
    integer :: num_spec_el
    integer :: num_gll
    integer :: num_global_points
    real(dp), allocatable :: global_points(:)
    real(dp), allocatable :: inv_metric(:, :)
    real(dp), allocatable :: metric_det(:, :)
    real(dp), allocatable :: rho(:)
    real(dp), allocatable :: mu(:)
    real(dp), allocatable :: mass_mat(:)
    integer, allocatable :: i_bool(:, :)
    real(dp), allocatable :: gll_weights(:)
    real(dp), allocatable :: h_prime(:, :)
    real(dp), allocatable :: circle_points(:, :)
    real(dp) :: r

    real(dp), allocatable :: displ(:), vel(:), accel(:)
  end type Mesh

  contains
!===============================================================================
! Subroutine to set up the global mesh and local ot global numbering
!===============================================================================
  subroutine Setup_mesh(this, num_spec_el, num_gll, r)
    type(Mesh) :: this
    integer num_spec_el, num_gll
    real(dp) r

    integer i_spec, i_gll, j_gll, i_glob
    real(dp), dimension(0:num_spec_el) :: anchor_points
    real(dp), dimension(num_gll) :: gll_points
    real(dp), external :: lagrange_deriv_GLL

    this%num_spec_el = num_spec_el
    this%num_gll = num_gll
    this%r = r
    this%num_global_points = num_spec_el * (num_gll - 1) + 1

    allocate(this%global_points(this%num_global_points))
    allocate(this%inv_metric(num_gll, num_spec_el))
    allocate(this%metric_det(num_gll, num_spec_el))
    allocate(this%rho(this%num_global_points))
    allocate(this%mu(this%num_global_points))
    allocate(this%mass_mat(this%num_global_points))
    allocate(this%i_bool(num_gll, num_spec_el))
    allocate(this%gll_weights(num_gll))
    allocate(this%h_prime(num_gll, num_gll))
    allocate(this%circle_points(this%num_global_points, 2))
    allocate(this%displ(this%num_global_points))
    allocate(this%vel(this%num_global_points))
    allocate(this%accel(this%num_global_points))

    ! GLL points and weights
    call zwgljd(gll_points,this%gll_weights,num_gll,0.0_dp,0.0_dp)
    if(mod(num_gll,2) /= 0) gll_points((num_gll-1)/2+1) = 0.0_dp

    ! get the derivatives of the Lagrange polynomials at 
    ! the GLL points; recall that  h_prime(i,j)=h'_{j}(xigll_{i}) 
    do i_gll=1,num_gll
       do j_gll=1,num_gll
          this%h_prime(i_gll,j_gll) = lagrange_deriv_GLL(j_gll-1,i_gll-1,gll_points,num_gll)
       end do
    end do

    ! Setup anchor points
    do i_spec = 0, num_spec_el
      anchor_points(i_spec) = 2 * PI * dble(i_spec) / dble(num_spec_el)
    enddo

    do i_spec = 1,num_spec_el
      this%metric_det(:, i_spec) = (2. * r / (anchor_points(i_spec)-anchor_points(i_spec-1)))**2
      this%inv_metric(:, i_spec) = ((anchor_points(i_spec)-anchor_points(i_spec-1)) / 2. * r)**2
    enddo

    ! Local to global numbering
    i_glob = 1
    do i_spec = 1,num_spec_el
      do i_gll = 1,num_gll
        if(i_gll > 1) i_glob = i_glob+1
        this%i_bool(i_gll,i_spec) = i_glob
      enddo
    enddo

    ! Compute global gridpoint position
    do i_spec = 1,num_spec_el
      do i_gll = 1,num_gll
        i_glob = this%i_bool(i_gll,i_spec)
        this%global_points(i_glob) = 0.5*(1.-gll_points(i_gll))*anchor_points(i_spec-1)&
          +0.5*(1.+gll_points(i_gll))*anchor_points(i_spec)
      enddo
    enddo

    ! Set mesh properties
    this%rho = density_fn(this%global_points)
    this%mu = rigidity_fn(this%global_points)

    this%mass_mat = mass_mat_glob(this)

    do i_glob = 1, this%num_global_points
      this%circle_points(i_glob, 1) = this%r * cos(this%global_points(i_glob))
      this%circle_points(i_glob, 2) = this%r * sin(this%global_points(i_glob))
    enddo

  end subroutine Setup_mesh

  function calculate_delta_t(this) result(delta_t)
    type(Mesh) this
    real(dp) delta_t, delta_h, max_vel
    real(dp) :: courant_CFL = 0.4d0 !CFL stability condition

    delta_h = 2 * PI * this%r / dble(this%num_global_points)
    max_vel = maxval(sqrt(this%mu / this%rho))
    delta_t = courant_CFL * delta_h / max_vel
  end function calculate_delta_t

  function mass_mat_glob(this) result(mass_mat)
    type(Mesh) this
    integer i_spec, i_gll, i_glob
    real(dp) mass_mat_local
    real(dp) mass_mat(this%num_global_points)
    do i_spec = 1, this%num_spec_el
      do i_gll = 1, this%num_gll
        i_glob = this%i_bool(i_gll,i_spec)
        mass_mat_local = this%gll_weights(i_gll)*(this%rho(i_glob))*sqrt(this%metric_det(i_gll,i_spec))
        mass_mat(i_glob) = mass_mat(i_glob) + mass_mat_local
      enddo
    enddo
  end function mass_mat_glob

  subroutine Increment_system(this, delta_t)
    type(Mesh) this
    real(dp) delta_t, du_dxi, epsilon, sigma
    integer i_spec, i_gll, j_gll, k_gll, i_glob
    real(dp) temp(this%num_gll)
    real(dp) templ
    real(dp) boundary_displ, boundary_vel, boundary_accel

    ! `Predictor' update displacement using explicit finite-difference time scheme (Newmark)
      this%displ(:) = this%displ(:) + delta_t * this%vel(:) + delta_t**2/2 * this%accel(:)
      this%vel(:) = this%vel(:) + delta_t/2  *this%accel(:)
      this%accel(:) = 0.
      do i_spec = 1,this%num_spec_el
        do i_gll = 1,this%num_gll
          ! Compute d(u) / d(xi)
          du_dxi = 0.
          do j_gll = 1,this%num_gll
            i_glob = this%i_bool(j_gll,i_spec)
            du_dxi = du_dxi + this%displ(i_glob)*this%h_prime(i_gll,j_gll)
          enddo
          i_glob = this%i_bool(i_gll, i_spec)
          temp(i_gll) = du_dxi*sqrt(this%metric_det(i_gll, i_spec))*this%mu(i_glob)*this%inv_metric(i_gll,i_spec)
        enddo ! first loop over the GLL points
        do k_gll = 1,this%num_gll
          templ = 0.
          do i_gll = 1,this%num_gll
            templ = templ + temp(i_gll)*this%h_prime(i_gll,k_gll)*this%gll_weights(i_gll)
          enddo
          ! `Corrector' update of acceleration in the Newmark scheme
          ! The minus sign comes from the integration by part done in the weak formulation of the equations
          i_glob = this%i_bool(k_gll,i_spec)
          this%accel(i_glob) = this%accel(i_glob) - templ
        enddo ! Second loop over the GLL points
      enddo ! End loop over all spectral elements
      ! Periodic BCs
      boundary_displ = (this%displ(1) + this%displ(this%num_global_points)) / 2
      boundary_vel = (this%vel(1) + this%vel(this%num_global_points)) / 2
      boundary_accel = (this%accel(1) + this%accel(this%num_global_points)) / 2
      this%displ(1) = boundary_displ
      this%displ(this%num_global_points) = boundary_displ
      this%vel(1) = boundary_vel
      this%vel(this%num_global_points) = boundary_vel
      this%accel(1) = boundary_accel
      this%accel(this%num_global_points) = boundary_accel

      ! Divide by the mass matrix, which is strictly (i.e. perfectly) diagonal
      this%accel(:) = this%accel(:)/this%mass_mat(:)
      ! `Corrector' update velocity
      this%vel(:) = this%vel(:) + delta_t/2 * this%accel(:)
  end subroutine Increment_system

  subroutine Output_snapshot(snapshot_dir, this, timestep)
    character(len=*) snapshot_dir
    type(Mesh) :: this
    integer timestep, i_glob
    character(len=50) snapshot_file

    write(snapshot_file, "('snapshot',i5.5)") timestep

    open(unit=10, file=snapshot_dir//snapshot_file, action="write", status="unknown")
    do i_glob = 1, this%num_global_points
      write(10,*) sngl(this%global_points(i_glob)),sngl(this%displ(i_glob))
    enddo
  end subroutine Output_snapshot

  subroutine Output_3d_snapshot(snapshot_dir, this, timestep)
    character(len=*) snapshot_dir
    type(Mesh) :: this
    integer timestep, i_glob
    character(len=50) snapshot_file

    write(snapshot_file, "('snapshot',i5.5)") timestep

    open(unit=10, file=snapshot_dir//snapshot_file, action="write", status="unknown")
    do i_glob = 1, this%num_global_points
      write(10,*) sngl(this%circle_points(i_glob, 1)), sngl(this%circle_points(i_glob, 2)), sngl(this%displ(i_glob))
    enddo
  end subroutine Output_3d_snapshot

!===============================================================================
! Function to determine the density on gridpoints
!===============================================================================
  function density_fn(mesh)
    real(dp) mesh(:)
    real(dp) density_fn(size(mesh))
    integer i_mesh
    do i_mesh = 1, size(mesh)
      density_fn(i_mesh) = 1d0
    enddo
  end function density_fn

!===============================================================================
! Function to determine the rigidity on gridpoints
!===============================================================================
  function rigidity_fn(mesh)
    real(dp) mesh(:)
    real(dp) rigidity_fn(size(mesh))
    integer i_mesh
    do i_mesh = 1, size(mesh)
      rigidity_fn(i_mesh) = 1d0
    enddo
  end function rigidity_fn

  subroutine Initial_conditions(this)
    type(Mesh) :: this
    integer i_displ, i_vel

    do i_displ = 1, this%num_global_points
      this%displ(i_displ) = exp(-10. * (this%global_points(i_displ) - this%global_points(this%num_global_points / 2))**2)
    enddo

    do i_vel = 1, this%num_global_points
      this%vel(i_vel) = 0.
    enddo
  end subroutine Initial_conditions

end module OneDimWaveModule
