module Parameters

  implicit none

  integer, parameter :: dp = kind(1.0d0)

! System parameters
  character(len = 43) :: output_path = "/home/davidho/WaveEqCurvedSpacetime/output/"

! Mathematical constants
  real(dp), parameter :: PI = 3.141592653589793d0

! Simulation parameters
  integer, parameter :: NUM_SPEC_EL = 200 ! No. of spectral elements
  integer, parameter :: NUM_GLL = 4 ! No. of Gauss-Lobatto-Legendre points
  integer, parameter :: NUM_TIMESTEPS = 100000 ! No. of time steps to calculate
  integer, parameter :: SNAPSHOT_TIMESTEP = 500 ! No. of Timesteps between snapshots

  integer, parameter :: NUM_GLOBAL_POINTS = NUM_SPEC_EL * (NUM_GLL - 1) + 1

! Model parameters
  real(dp), parameter :: LENGTH = 10000 ! 

  contains
!===============================================================================
! Function to determine the density on gridpoints
!===============================================================================
  function density_fn(mesh)
    real(dp) mesh(:)
    real(dp) density_fn(size(mesh))
    integer i_mesh
    do i_mesh = 1, size(mesh)
      density_fn(i_mesh) = 1.5d3
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
      rigidity_fn(i_mesh) = 3.0d13
    enddo
  end function rigidity_fn

  subroutine Initial_conditions(global_points, displ, vel)
    real(dp) global_points(:), displ(:), vel(:)
    integer i_displ, i_vel

    do i_displ = 1, size(displ)
      displ(i_displ) = cos(8 * PI * global_points(i_displ) / LENGTH)
    enddo

    do i_vel = 1, size(vel)
      vel(i_vel) = 0.d0
    enddo
  end subroutine Initial_conditions

end module Parameters
