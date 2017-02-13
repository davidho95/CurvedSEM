
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

  ! use Parameters
  ! use OneDimWaveModule
  use TwoDimWaveModule

  implicit none

  integer, dimension(NUM_GLL, NUM_SPEC_EL) :: i_bool
  real(dp), dimension(NUM_GLL, NUM_SPEC_EL) :: jacobian_mat, jacobian
  real(dp), dimension(NUM_GLOBAL_POINTS) :: global_points, rho, mu
  real(dp), dimension(NUM_GLOBAL_POINTS) :: displ, vel, accel, mass_mat
  real(dp), dimension(NUM_GLL, NUM_GLL) :: h_prime
  real(dp), dimension(NUM_GLL) :: gll_weights
  real(dp) delta_t

! Index variables
  integer timestep

  call Setup_mesh(i_bool, rho, mu, global_points, gll_weights, h_prime, jacobian_mat, jacobian)

  delta_t = calculate_delta_t(rho, mu)

  mass_mat = mass_mat_glob(rho, jacobian, i_bool, gll_weights)

  accel(:) = 0

  call Initial_conditions(global_points, displ, vel)

  do timestep = 1, NUM_TIMESTEPS
    call Increment_system(displ, vel, accel, rho, mu, mass_mat, i_bool, h_prime,&
      gll_weights, jacobian_mat, jacobian, delta_t)
    if (mod(timestep, SNAPSHOT_TIMESTEP) == 0) call Output_3d_snapshot(OUTPUT_PATH, global_points, displ, timestep)
  enddo

end program OneDimWaveEqn