
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

  use Parameters

  implicit none

! Index variables
  integer timestep

  call Setup_mesh()

  ! call Inital_conditions()


  ! do timestep = 1, NUM_TIMESTEPS
  !   call Increment_system()
  !   if (mod(timestep, SNAPSHOT_TIMESTEP) == 0) call Output_snapshot()
  ! enddo

end program OneDimWaveEqn

!===============================================================================
! Subroutine to set up the global mesh and local ot global numbering
!===============================================================================
subroutine Setup_mesh()

  use Parameters

  implicit none

  integer, parameter :: NUM_GLOBAL_POINTS = NUM_SPEC_EL * (NUM_GLL - 1) + 1
  integer i_spec, i_gll, i_glob
  real(dp), dimension(0:NUM_SPEC_EL) :: anchor_points
  real(dp), dimension(NUM_GLL) :: gll_points, gll_weights
  real(dp), dimension(NUM_GLOBAL_POINTS) :: global_points, rho, mu
  integer, dimension(NUM_GLL, NUM_SPEC_EL) :: i_bool

  ! GLL points and weights
  call zwgljd(gll_points,gll_weights,NUM_GLL,0.0_dp,0.0_dp)
  if(mod(NUM_GLL,2) /= 0) gll_points((NUM_GLL-1)/2+1) = 0.0_dp

  ! Setup anchor points
  do i_spec = 0, NUM_SPEC_EL
    anchor_points(i_spec) = dble(i_spec) * LENGTH / dble(NUM_SPEC_EL)
  enddo

  ! Local to global numbering
  i_glob = 1
  do i_spec = 1,NUM_SPEC_EL
    do i_gll = 1,NUM_GLL
      if(i_gll > 1) i_glob = i_glob+1
      i_bool(i_gll,i_spec) = i_glob
    enddo
  enddo

! Compute global gridpoint position
  do i_spec = 1,NUM_SPEC_EL
    do i_gll = 1,NUM_GLL
      i_glob = i_bool(i_gll,i_spec)
      global_points(i_glob) = 0.5*(1.-gll_points(i_gll))*anchor_points(i_spec)&
        +0.5*(1.+gll_points(i_gll))*anchor_points(i_spec - 1)
    enddo
  enddo

  ! Set mesh properties
  rho = density_fn(global_points)
  mu = rigidity_fn(global_points)

end subroutine Setup_mesh