
!===============================================================================
!
! MeshClass.f90
!
! Pseudo-class for a 1d spectral element mesh
!
! David Ho (2017)
! 
! Adapted from Specfem1d v1.0, Komatisch & Tromp (2007)
!
!===============================================================================
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!===============================================================================

module MeshClass

  implicit none

  integer, parameter :: dp = kind(1.0d0)

! Mathematical constants
  real(dp), parameter :: PI = 3.141592653589793d0

  type Mesh
    integer :: num_spec_el
    integer :: num_gll
    integer :: num_global_points
    real(dp), allocatable :: global_points(:)
    !Inverse metric is used more often than metric
    real(dp), allocatable :: inv_metric(:, :)
    real(dp), allocatable :: metric_det(:, :)
    real(dp), allocatable :: rho(:) !Density
    real(dp), allocatable :: mu(:) !rigidity
    real(dp), allocatable :: mass_mat(:)
    integer, allocatable :: i_bool(:, :) !For local to global numbering
    real(dp), allocatable :: gll_weights(:)
    !Derivatives of Lagrange polynomials at GLL points
    real(dp), allocatable :: h_prime(:, :)
    real(dp), allocatable :: circle_points(:, :) 
    real(dp) :: r

    real(dp), allocatable :: displ(:), vel(:), accel(:)
  end type Mesh

  contains
!===============================================================================
! Subroutine to set up the global mesh and local to global numbering
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
          this%h_prime(i_gll,j_gll)&
          = lagrange_deriv_GLL(j_gll-1,i_gll-1,gll_points,num_gll)
       end do
    end do

    ! Set Metric
    do i_spec = 1,num_spec_el
      this%metric_det(:, i_spec) = (2. * r * num_spec_el / (2 * PI))**2
      this%inv_metric(:, i_spec) = 1 / this%metric_det(:, i_spec)
    enddo

    ! Local to global numbering
    i_glob = 1
    do i_spec = 1,num_spec_el
      do i_gll = 1,num_gll
        if(i_gll > 1) i_glob = i_glob+1
        this%i_bool(i_gll,i_spec) = i_glob
      enddo
    enddo

    ! Setup anchor points
    do i_spec = 0, num_spec_el
      anchor_points(i_spec) = 2 * PI * dble(i_spec) / dble(num_spec_el)
    enddo

    ! Compute global gridpoint position
    do i_spec = 1,num_spec_el
      do i_gll = 1,num_gll
        i_glob = this%i_bool(i_gll,i_spec)
        this%global_points(i_glob)&
          = 0.5*(1.-gll_points(i_gll))*anchor_points(i_spec-1)&
          +0.5*(1.+gll_points(i_gll))*anchor_points(i_spec)
      enddo
    enddo

    ! Set mesh properties
    call set_density(this)
    call set_rigidity(this)
    this%mass_mat = set_mass_matrix(this)

    do i_glob = 1, this%num_global_points
      this%circle_points(i_glob, 1) = this%r * cos(this%global_points(i_glob))
      this%circle_points(i_glob, 2) = this%r * sin(this%global_points(i_glob))
    enddo

  end subroutine Setup_mesh

!==============================================================================
! Function to calculate global mass matrix
!==============================================================================
  function set_mass_matrix(this) result(mass_mat)
    type(Mesh) this
    integer i_spec, i_gll, i_glob
    real(dp) mass_mat_local
    real(dp) mass_mat(this%num_global_points)
    do i_spec = 1, this%num_spec_el
      do i_gll = 1, this%num_gll
        i_glob = this%i_bool(i_gll,i_spec)
        mass_mat_local = this%gll_weights(i_gll)*(this%rho(i_glob))&
        *sqrt(this%metric_det(i_gll,i_spec))
        mass_mat(i_glob) = mass_mat(i_glob) + mass_mat_local
      enddo
    enddo
  end function set_mass_matrix

!===============================================================================
! Function to determine the density on gridpoints
!===============================================================================
  subroutine set_density(this)
    type(Mesh) this
    integer i_glob

    this%rho = 1d0
  end subroutine set_density

!===============================================================================
! Function to determine the rigidity on gridpoints
!===============================================================================
  subroutine set_rigidity(this)
    type(Mesh) this
    integer i_glob

    this%mu = 1d0
  end subroutine set_rigidity

end module MeshClass
