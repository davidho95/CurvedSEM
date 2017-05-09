
!===============================================================================
!
! OneDimWaveEqn.f90
!
! Solves the wave equation on a (1 + 1) dimensional pseudo-Riemannian manifold
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

!===============================================================================
! Main program
!===============================================================================
program OneDimWaveEqn

  use MeshClass
  use WaveModule

  implicit none

  type(Mesh) my_mesh
  real(dp) delta_t
  integer :: num_spec_el = 50 ! No. of spectral elements
  integer :: num_gll = 4 ! No. of Gauss-Lobatto-Legendre points
  integer :: num_timesteps = 4000 ! No. of time steps to calculate
  integer :: snapshot_timestep = 40 ! No. of Timesteps between snapshots
  character(len = 43) :: output_path = "/home/davidho/CurvedSEM/CurvedSEM1d/output/"

  integer i_spec, i_gll, i_glob

! Model parameters
  real(dp) :: radius = 1d0 ! 


! Index variables
  integer timestep

  call Setup_mesh(my_mesh, num_spec_el, num_gll, radius)

  delta_t = calculate_delta_t(my_mesh)

  call Initial_conditions(my_mesh)

  do timestep = 1, NUM_TIMESTEPS
    call Increment_system(my_mesh, 4.9d0)
    if (mod(timestep, SNAPSHOT_TIMESTEP) == 0) then
      call Output_3d_snapshot(OUTPUT_PATH, my_mesh, timestep)
      !print *, calculate_energy(my_mesh)
    endif
  enddo

end program OneDimWaveEqn