module OneDimParameters

  implicit none

  integer, parameter :: dp = kind(1.0d0)

! System parameters
  character(len = 43) :: output_path = "/home/davidho/WaveEqCurvedSpacetime/output/"

! Mathematical constants
  real(dp), parameter :: PI = 3.141592653589793d0

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

end module OneDimParameters